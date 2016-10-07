/*
 
 Copyright (c) 2012-2015, Michael (Mikhail) Yudelson
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of the Michael (Mikhail) Yudelson nor the
 names of other contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS AND CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */

#include "HMMProblemComp.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemComp::HMMProblemComp(struct param *param) {
    for(NPAR i=0; i<3; i++) this->sizes[i] = param->nK;
//    this->sizes = {param->nK, param->nK, param->nK};
    this->n_params = param->nK * 4;
  
    init(param);
}

void HMMProblemComp::init(struct param *param) {
	this->p = param;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, (size_t)this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    
    NUMBER *a_PI, ** a_A, ** a_B;
    init3Params(a_PI, a_A, a_B, nS, nO);
    
    //
    // setup params
    //
	NPAR i, j, idx, offset;
	NUMBER sumPI = 0;
	NUMBER sumA[nS];
	NUMBER sumB[nS];
	for(i=0; i<nS; i++) {
		sumA[i] = 0;
		sumB[i] = 0;
	}
	// populate PI
	for(i=0; i<((nS)-1); i++) {
		a_PI[i] = this->p->init_params[i];
		sumPI  += this->p->init_params[i];
	}
	a_PI[nS-1] = 1 - sumPI;
	// populate A
	offset = (NPAR)(nS-1);
	for(i=0; i<nS; i++) {
		for(j=0; j<((nS)-1); j++) {
			idx = (NPAR)(offset + i*((nS)-1) + j);
			a_A[i][j] = this->p->init_params[idx];
			sumA[i]  += this->p->init_params[idx];
		}
		a_A[i][((this->p->nS)-1)]  = 1 - sumA[i];
	}
	// polupale B
	offset = (NPAR)((nS-1) + nS*(nS-1));
	for(i=0; i<nS; i++) {
		for(j=0; j<((nO)-1); j++) {
			idx = (NPAR)(offset + i*((nO)-1) + j);
			a_B[i][j] = this->p->init_params[idx];
			sumB[i] += this->p->init_params[idx];
		}
		a_B[i][((nO)-1)]  = 1 - sumB[i];
	}
    
    // mass produce PI's/PIg's, A's, B's
    if( this->p->do_not_check_constraints==0 && !checkPIABConstraints(a_PI, a_A, a_B)) {
        fprintf(stderr,"params do not meet constraints.\n");
        exit(1);
    }
    this->pi  = init2D<NUMBER>(nK, nS);
    this->A   = init3D<NUMBER>(nK, nS, nS);
    this->B   = init3D<NUMBER>(nK, nS, nO);
    NCAT x;
    for(x=0; x<nK; x++) {
        cpy1D<NUMBER>(a_PI, this->pi[x], nS);
        cpy2D<NUMBER>(a_A,  this->A[x],  nS, nS);
        cpy2D<NUMBER>(a_B,  this->B[x],  nS, nO);
    }
    // PIg start with same params
    // destroy setup params
	free(a_PI);
	free2D<NUMBER>(a_A, nS);
	free2D<NUMBER>(a_B, nS);
	
    // if needs be -- read in init params from a file
    if(param->initfile[0]!=0)
        this->readModel(param->initfile, false /* read and upload but not overwrite*/);
    
    // populate boundaries
	// populate lb*/ub*
	// *PI
    init3Params(this->lbPI, this->lbA, this->lbB, nS, nO);
    init3Params(this->ubPI, this->ubA, this->ubB, nS, nO);
	for(i=0; i<nS; i++) {
		lbPI[i] = this->p->param_lo[i];
		ubPI[i] = this->p->param_hi[i];
	}
	// *A
	offset = nS;
	for(i=0; i<nS; i++)
		for(j=0; j<nS; j++) {
			idx = (NPAR)(offset + i*nS + j);
			lbA[i][j] = this->p->param_lo[idx];
			ubA[i][j] = this->p->param_hi[idx];
		}
	// *B
	offset = (NPAR)(nS + nS*nS);
	for(i=0; i<nS; i++)
		for(j=0; j<nO; j++) {
			idx = (NPAR)(offset + i*nS + j);
			lbB[i][j] = this->p->param_lo[idx];
			ubB[i][j] = this->p->param_hi[idx];
		}
    // is_multi
    this->is_multi = Calloc(NPAR, (size_t)this->p->nK);
    // assign is_multi labels
    NCAT k, n;
    for(NDAT t=0; t<this->p->N; t++) {
        n = this->p->dat_skill_rcount[t];
        if(n>1) {
            for(NCAT l=0; l<n; l++) {
                k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[t] + l ];
                is_multi[k] = 1;
            }
        }
    }
}

HMMProblemComp::~HMMProblemComp() {
    destroy();
}

void HMMProblemComp::destroy() {
	// destroy additional model data
	free(this->is_multi);
}// ~HMMProblemComp

//NUMBER** HMMProblemComp::getPI() {
//	return this->pi;
//}
//
//NUMBER* HMMProblemComp::getPI(NCAT x) { // same as getPIk(x)
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->pi[x];
//}
//
//NUMBER** HMMProblemComp::getA(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->A[x];
//}
//
//NUMBER** HMMProblemComp::getB(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->B[x];
//}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getPI(struct data* dt, NPAR i) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if dt->cnt==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question
	if( dt->cnt==0 && this->is_multi[dt->k] == 0) { // this skill is never multiskill
		r = this->pi[dt->k][i];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->pi[dt->k][i];
		} else { // collect multiple \pi[i]
			NUMBER* p = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				p[l] = this->pi[k][i];
				r = squishing(p, n_skills);
			}
		} // definitely multiskill
	} // end sometimes multiskill
    return r;
}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getA(struct data* dt, NPAR i, NPAR j) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if dt->cnt==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question
	if( dt->cnt==0 && this->is_multi[dt->k] == 0) { // this skill is never multiskill
		r = this->A[dt->k][i][j];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->A[dt->k][i][j];
		} else { // collect multiple \pi[i]
			NUMBER* p = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				p[l] = this->A[k][i][j];
				r = squishing(p, n_skills);
			}
		} // definitely multiskill
	} // end sometimes multiskill
	return r;
}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getB(struct data* dt, NPAR i, NPAR m) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if dt->cnt==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question

	// special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
	// in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
	if(m<0) {
		r = 1;
	} else if( dt->cnt==0 && this->is_multi[dt->k] == 0) { // this skill is never multiskill
		r = this->B[dt->k][i][m];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->B[dt->k][i][m];
		} else { // collect multiple \pi[i]
			NUMBER* p = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				p[l] = this->B[k][i][m];
				r = squishing(p, n_skills);
			}
		} // definitely multiskill
	} // end sometimes multiskill
	return r;
}


void HMMProblemComp::setGradPI(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0, ndat = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
//    o = dt->obs[t];
        o = this->p->dat_obs[ dt->ix[t] ];//->get( dt->ix[t] );
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( fb->pi[i] * (1-fb->pi[i]) );
            fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_PI, this->p, (NUMBER)ndat);
}

void HMMProblemComp::toFile(const char *filename) {
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    
    // write solved id
    writeSolverInfo(fid, this->p);
    
	fprintf(fid,"Null skill ratios\t");
	for(NPAR m=0; m<this->p->nO; m++)
		fprintf(fid," %12.10f%s",this->null_obs_ratio[m],(m==(this->p->nO-1))?"\n":"\t");
	NCAT k, g;
    NPAR i,j,m;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<this->p->nG; g++) {
		it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		fprintf(fid,"PIg\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->PIg[g][i],(i==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->pi[k][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%12.10f%s",this->A[k][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%12.10f%s",this->B[k][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblemComp::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure == STRUCTURE_PIgk) {
        loglik_rmse[0] += GradientDescent();
    } else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemComp::GradientDescent() {
	NCAT k, g, /*ki, gi, nX, */x;
    NCAT nK = this->p->nK, nG = this->p->nG;

	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill>0) {
        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
        fb->link( HMMProblem::getPI(0), HMMProblem::getA(0), HMMProblem::getB(0), this->p->nSeq, this->p->k_data);// link skill 0 (we'll copy fit parameters to others
        if(this->p->block_fitting[0]!=0) fb->pi = NULL;
        if(this->p->block_fitting[1]!=0) fb->A  = NULL;
        if(this->p->block_fitting[2]!=0) fb->B  = NULL;
        
        fb->init(FBS_PARm1);
        fb->init(FBS_PARm2);
        fb->init(FBS_GRAD);
        if(this->p->solver==METHOD_CGD) {
            fb->init(FBS_GRADm1);
            fb->init(FBS_DIRm1);
        }
        
        NCAT* original_ks = Calloc(NCAT, (size_t)this->p->nSeq);
        for(x=0; x<this->p->nSeq; x++) { original_ks[x] = this->p->all_data[x].k; this->p->all_data[x].k = 0; } // save progonal k's
        FitResult fr = GradientDescentBit(fb);
        for(x=0; x<this->p->nSeq; x++) { this->p->all_data[x].k = original_ks[x]; } // restore original k's
        free(original_ks);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
        delete fb;
    }
	//
	// Main fit
	//
    if( this->p->single_skill!=2 ) { // if not "force single skill"
        NCAT first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
        NCAT iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
        NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK);
        NCAT* iter_qual_group = Calloc(NCAT, (size_t)nG);
        NCAT* iter_fail_skill = Calloc(NCAT, (size_t)nK); // counting concecutive number of failures to converge for skills
        NCAT* iter_fail_group = Calloc(NCAT, (size_t)nG); // counting concecutive number of failures to converge for groups
        int skip_k = 0, skip_g = 0;

        int i = 0; // count runs
//        int parallel_now = this->p->parallel==1; //PAR
//        #pragma omp parallel if(parallel_now) shared(iter_qual_group,iter_qual_skill)//PAR
//        {//PAR
        while(skip_k<nK || skip_g<nG) {
            //
            // Skills first
            //
            
            if(skip_k<nK) {
//                #pragma omp for schedule(dynamic) //PAR
                for(k=0; k<nK; k++) { // for all A,B-by-skill
                    if(iter_qual_skill[k]==iterations_to_qualify || iter_fail_skill[k]==iterations_to_qualify)
                        continue;
                    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                    // link
                    fb->link( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), this->p->k_numg[k], this->p->k_g_data[k]);// link skill 0 (we'll copy fit parameters to others
                    if(this->p->block_fitting[0]!=0) fb->pi = NULL;
                    if(this->p->block_fitting[1]!=0) fb->A  = NULL;
                    if(this->p->block_fitting[2]!=0) fb->B  = NULL;

                    fb->Cslice = 0;
                    fb->init(FBS_PARm1);
                    fb->init(FBS_PARm2);
                    fb->init(FBS_GRAD);
                    if(this->p->solver==METHOD_CGD) {
                        fb->init(FBS_GRADm1);
                        fb->init(FBS_DIRm1);
                    }
                    FitResult fr = GradientDescentBit(fb);
                    // decide on convergence
                    if(i>=first_iteration_qualify || fb->xndat==0) {
                        if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==nG || fb->xndat==0) { // converged quick, or don't care (others all converged
                            iter_qual_skill[k]++;
                            if(fb->xndat==0) {
                                iter_qual_skill[k]=iterations_to_qualify;
                                fr.conv = 1;
                            }
                            if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
                                if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
//                                #pragma omp critical(update_skip_k)//PAR
//                                {//PAR
                                    skip_k++;
//                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                    printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        }
                        else
                            iter_qual_skill[k]=0;
                    } // decide on convergence
                    delete fb;
                } // for all skills
            }
            //
            // PIg second
            //
            if(skip_g<nG){
//                #pragma omp for schedule(dynamic)//PAR
                for(g=0; g<nG; g++) { // for all PI-by-user
                    if(iter_qual_group[g]==iterations_to_qualify)
                        continue;
                    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                    // link
                    fb->link(this->getPIg(g), NULL, NULL, this->p->g_numk[g], this->p->g_k_data[g]);
                    if(this->p->block_fitting[0]!=0) fb->pi = NULL;
                    if(this->p->block_fitting[1]!=0) fb->A  = NULL;
                    if(this->p->block_fitting[2]!=0) fb->B  = NULL;

                    fb->Cslice = 1;
                    fb->init(FBS_PARm1);
                    fb->init(FBS_PARm2);
                    fb->init(FBS_GRAD);
                    if(this->p->solver==METHOD_CGD) {
                        fb->init(FBS_GRADm1);
                        fb->init(FBS_DIRm1);
                    }
                    FitResult fr = GradientDescentBit(fb);
                    // decide on convergence
                    if(i>=first_iteration_qualify || fb->xndat==0) { //can qualify or  student had no skill labelled rows
                        if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==nK || fb->xndat==0) { // converged quick, or don't care (others all converged), or
                            iter_qual_group[g]++;
                            if(fb->xndat==0) {
                                iter_qual_group[g]=iterations_to_qualify;
                                fr.conv = 1;
                            }
                            if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                                if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
//                                #pragma omp critical(update_skip_g)//PAR
//                                {//PAR
                                    skip_g++;
//                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                    printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        }
                        else
                            iter_qual_group[g]=0;
                    } // decide on convergence
                    delete fb;
                } // for all groups
            }
//            fr.conv = checkConvergenceBig(fbs, nK+nG, this->p->tol, &crit);
//            
//            computeAlphaAndPOParam(this->p->nSeq, this->p->k_data);
//            ll = HMMProblem::getSumLogPOPara(this->p->nSeq, this->p->k_data);
//            printf("*%i ll=%15.7f, crit=%12.10f\n",i,ll,crit);
//            #pragma omp single//PAR
//            {//PAR
            i++;
//            }//PAR
        }
//        }//PAR
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_qual_group != NULL) free(iter_qual_group);
        if( iter_fail_skill != NULL ) free(iter_fail_skill);
        if( iter_fail_group != NULL) free(iter_fail_group);

    } // if not "force single skill
        
//    for(x=0;x<nX;x++) delete fbs[x];
//    if(fbs!=NULL) free(fbs);
        
    if(false){
        NCAT q, x;
        NCAT nQ = this->p->nK + this->p->nG;
        FitBit **fbs = Calloc(FitBit*, (size_t)nQ);
        for(q=0;q<nQ;q++) {
            NCAT xndat;
            struct data** x_data;
            if(q<this->p->nK) { // skills
                x = q;
                xndat = this->p->k_numg[x];
                x_data = this->p->k_g_data[x];
                fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                fbs[q]->link( HMMProblem::getPI(x), HMMProblem::getA(x), HMMProblem::getB(x), xndat, x_data);
            } else { // groups
                x = q - this->p->nK;
                xndat = this->p->g_numk[x];
                x_data = this->p->g_k_data[x];
                fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                fbs[q]->link( this->getPIg(x), NULL, NULL, xndat, x_data);
            }
            fbs[q]->init(FBS_PARm1);
            fbs[q]->init(FBS_GRAD);
            if(this->p->solver==METHOD_CGD) {
                fbs[q]->init(FBS_GRADm1);
                fbs[q]->init(FBS_DIRm1);
            }
        }
        /*FitResult fr=*/GradientDescentBitBig(fbs, nQ);
        
        for(q=0;q<nQ;q++) delete fbs[q];
        if(fbs!=NULL) free(fbs);
    }
    
    // compute loglik
    return getSumLogPOPara(this->p->nSeq, this->p->k_data);
}


//struct data***  HMMProblemComp::getExdendedData(NCAT xndat, struct data** x_data, NPAR kg_flag, NCAT* xxndat) {
//    NCAT q, r, g, k;
//    struct data*** xx_data;
//    *xxndat = 0;
//    NCAT yndat;
//    struct data** z_data;
//    struct data** y_data;
//    if(kg_flag==0) {
//        for(q=0; q<xndat; q++)
//            *xxndat += this->p->g_numk[ x_data[q]->g ];
//        
//    } else if (kg_flag==1) {
//        for(q=0; q<xndat; q++)
//            *xxndat += this->p->k_numg[ x_data[q]->k ];
//        
//    } else {
//        fprintf(stderr,"Unexpected value of kg_flag\n");
//        return NULL;
//    }
//    z_data = (Malloc(struct data *, (size_t)*xxndat));
//    xx_data = &z_data;
//    NDAT i = 0;
//    if(kg_flag==0) {
//        for(q=0; q<xndat; q++) {
//            g = x_data[q]->g;
//            yndat = this->p->g_numk[g];
//            y_data = this->p->g_k_data[g];
//            for(r=0; r<yndat; r++)
//                z_data[i++] = y_data[r];
//        }
//        
//    } else if (kg_flag==1) {
//        for(q=0; q<xndat; q++) {
//            k = x_data[q]->k;
//            yndat = this->p->k_numg[k];
//            y_data = this->p->k_g_data[k];
//            for(r=0; r<yndat; r++)
//                z_data[i++] = y_data[r];
//        }
//    }
//    return xx_data;
//}

//NUMBER HMMProblemComp::GradientDescentX() {
//	NCAT k, ki, g, gi;
//    /*NPAR nS = this->p->nS, nO = this->p->nO;*/ NCAT nK = this->p->nK, nG = this->p->nG;
//    NUMBER loglik = 0;
//    FitResult fr;
//    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
//    fb->init(FBS_PARm1);
//    fb->init(FBS_GRAD);
//    if(this->p->solver==METHOD_CGD) {
//        fb->init(FBS_GRADm1);
//        fb->init(FBS_DIRm1);
//    }
//	//
//	// fit all as 1 skill first, set group gradients to 0, and do not fit them
//	//
//	if(this->p->single_skill>0) {
//        fb->linkPar( HMMProblem::getPI(0), HMMProblem::getA(0), HMMProblem::getB(0), this->p->nSeq, this->p->k_data);// link skill 0 (we'll copy fit parameters to others
//        fr = GradientDescentBit(fb, true /*is1SkillForAll*/);
//        if( !this->p->quiet )
//            printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
//    }
//	//
//	// Main fit
//	//
//    if( this->p->single_skill!=2 ) { // if not "force single skill"
//        NCAT first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
//        NCAT iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
//        NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK);
//        NCAT* iter_qual_group = Calloc(NCAT, (size_t)nG);
//        int skip_k = 0, skip_g = 0;
//        
//        //        NUMBER **cpyPI;
//        ////        NUMBER **cpyPIg;
//        //        NUMBER ***cpyA;
//        //        NUMBER ***cpyB;
//        //        if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
//        //            cpyPI  = init2D<NUMBER>(this->p->nK, this->p->nS);
//        //            cpyA   = init3D<NUMBER>(this->p->nK, this->p->nS, this->p->nS);
//        //            cpyB   = init3D<NUMBER>(this->p->nK, this->p->nS, this->p->nO);
//        ////            cpyPIg = init2D<NUMBER>(this->p->nG, this->p->nS);
//        //            NCAT x;
//        //            for(x=0; x<nK; x++) {
//        //                cpy1D<NUMBER>(this->getPI(x), cpyPI[x], this->p->nS);
//        //                cpy2D<NUMBER>(this->getA(x),  cpyA[x],  this->p->nS, this->p->nS);
//        //                cpy2D<NUMBER>(this->getB(x),  cpyB[x],  this->p->nS, this->p->nO);
//        //            }
//        ////            // PIg start with same params
//        ////            for(x=0; x<nG; x++)
//        ////                cpy1D<NUMBER>(cpyPIg[x], this->PIg[x], this->p->nS);
//        //        }
//        
//        // utilize fitting larger data first
//        NDAT newnK=0, newnG=0;
//        this->reorderSequences(&newnK, &newnG);
//        
//        
//        int i = 0; // count runs
//        while(skip_k<newnK || skip_g<newnG) {
//            //
//            // Skills first
//            //
//            for(ki=0; ki<newnK && skip_k<newnK; ki++) { // for all A,B-by-skill
//                k = sortstrip_k[ki].id; // grab reordered k
//                
//                if(iter_qual_skill[k]==iterations_to_qualify)
//                    continue;
//                NCAT xndat = this->p->k_numg[k];
//                struct data** x_data = this->p->k_g_data[k];
//                
//                NCAT yndat;
//                struct data** y_data;
//                y_data = *getExdendedData(xndat, x_data, 0/*skill*/, &yndat);
//                xndat = yndat;
//                x_data = y_data;
//                
//                // link and fit
//                fb->linkPar( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k));// link skill 0 (we'll copy fit parameters to others
//                fr = GradientDescentBit(k/*use skill x*/, xndat, x_data, 0/*skill*/, fb, false /*is1SkillForAll*/);
//                // decide on convergence
//                if(i>=first_iteration_qualify) {
//                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==newnG) { // converged quick, or don't care (others all converged
//                        iter_qual_skill[k]++;
//                        if(iter_qual_skill[k]==iterations_to_qualify || skip_g==newnG) {// criterion met, or don't care (others all converged)
//                            if(skip_g==newnG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
//                            skip_k++;
//                            if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
//                                computeAlphaAndPOParam(xndat, x_data);
//                                printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
//                            }
//                        }
//                    }
//                    else
//                        iter_qual_skill[k]=0;
//                } // decide on convergence
//                free(y_data);
//                
//                //                //
//                //                // make copies of parameters to do gradient descend, link to copies, but grads computed from actual params
//                //                //
//                //                if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
//                //                    swap1D<NUMBER>(cpyPI[k], this->getPI(k), this->p->nS);
//                //                    swap2D<NUMBER>(cpyA[k],   this->getA(k), this->p->nS, this->p->nS);
//                //                    swap2D<NUMBER>(cpyB[k],   this->getB(k), this->p->nS, this->p->nO);
//                //                }
//            } // for all skills
//            //
//            // PIg second
//            //
//            //            int z = 0;
//            for(gi=0; gi<newnG && skip_g<newnG; gi++) { // for all PI-by-user
//                g = sortstrip_g[gi].id; // grab reordered g
//                
//                if(iter_qual_group[g]==iterations_to_qualify)
//                    continue;
//                NCAT xndat = this->p->g_numk[g];
//                struct data** x_data = this->p->g_k_data[g];
//                
//                NCAT yndat;
//                struct data** y_data;
//                y_data = *getExdendedData(xndat, x_data, 1/*group*/, &yndat);
//                xndat = yndat;
//                x_data = y_data;
//                
//                // vvvvvvvvvvvvvvvvvvvvv ONLY PART THAT IS DIFFERENT FROM others
//                fb->linkPar(this->getPIg(g), NULL, NULL);
//                // ^^^^^^^^^^^^^^^^^^^^^
//                // decide on convergence
//                fr = GradientDescentBit(g/*use group x*/, xndat, x_data, 1/*group*/, fb, false /*is1SkillForAll*/);
//                if(i>=first_iteration_qualify || xndat==0) { //can qualify or  student had no skill labelled rows
//                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==newnK || xndat==0) { // converged quick, or don't care (others all converged), or  student had no skill labelled rows
//                        iter_qual_group[g]++;
//                        if(xndat==0) {
//                            iter_qual_group[g]=iterations_to_qualify;
//                            fr.conv = 1;
//                        }
//                        if(iter_qual_group[g]==iterations_to_qualify || skip_k==newnK) {// criterion met, or don't care (others all converged)
//                            if(skip_k==newnK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
//                            skip_g++;
//                            if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
//                                computeAlphaAndPOParam(xndat, x_data);
//                                printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
//                            }
//                        }
//                    }
//                    else
//                        iter_qual_group[g]=0;
//                } // decide on convergence
//                free(y_data);
//                //                //
//                //                // make copies of parameters to do gradient descend, link to copies, but grads computed from actual params
//                //                //
//                //                if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
//                //                    swap1D<NUMBER>(this->PIg[g], cpyPIg[g], this->p->nS);
//                //                }
//            } // for all groups
//            loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
//            i++;
//        }
//        // recycle qualifications
//        if( iter_qual_skill != NULL ) free(iter_qual_skill);
//        if( iter_qual_group != NULL) free(iter_qual_group);
//        
//        //        if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
//        //            free2D<NUMBER>(cpyPI, nK);
//        //            free3D<NUMBER>(cpyA,  nK, this->p->nS);
//        //            free3D<NUMBER>(cpyB,  nK, this->p->nS);
//        ////            free2D<NUMBER>(cpyPIg, this->p->nG);
//        //        }
//    } // if not "force single skill
//    
//    delete fb;
//    // compute loglik
//    loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
//    return loglik;
//}
//
void HMMProblemComp::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
	NPAR i,j,m;
	NCAT k = 0, g = 0, idxk = 0, idxg = 0;
	string s;
    std::map<std::string,NCAT>::iterator it;
    char col[2048];
    //
    readNullObsRatio(fid, param, line_no);
    //
    // init param
    //
    if(overwrite) {
        this->p->map_group_fwd = new map<string,NCAT>();
        this->p->map_group_bwd = new map<NCAT,string>();
        this->p->map_skill_fwd = new map<string,NCAT>();
        this->p->map_skill_bwd = new map<NCAT,string>();
    }
	//
	// read grouped PIg
	//
    for(g=0; g<this->p->nG; g++) {
		// read group label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;        
        if(overwrite) {
            this->p->map_group_fwd->insert(pair<string,NCAT>(s, (NCAT)this->p->map_group_fwd->size()));
            this->p->map_group_bwd->insert(pair<NCAT,string>((NCAT)this->p->map_group_bwd->size(), s));
            idxg = g;
        } else {
            it = this->p->map_group_fwd->find(s);
            if( it==this->p->map_group_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                continue; // skip this iteration
            }
            else
                idxg =it->second;
        }

        // read PIg
        fscanf(fid,"PIg\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->PIg[idxg][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->PIg[idxg][i] = atof(col);
        (*line_no)++;
    }
    //
    // read skills
    //
	for(k=0; k<this->p->nK; k++) {
		// read skill label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;
        if(overwrite) {
            this->p->map_skill_fwd->insert(pair<string,NCAT>(s, (NCAT)this->p->map_skill_fwd->size()));
            this->p->map_skill_bwd->insert(pair<NCAT,string>((NCAT)this->p->map_skill_bwd->size(), s));
            idxk = k;
        } else {
            it = this->p->map_skill_fwd->find(s);
            if( it==this->p->map_skill_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                continue; // skip this iteration
            }
            else
                idxk =it->second;
        }
        
        // read PI
        fscanf(fid,"PIk\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->pi[idxk][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->pi[idxk][i] = atof(col);
        (*line_no)++;
		// read A
        fscanf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++) {
                if(i==(this->p->nS-1) && j==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->A[idxk][i][j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->A[idxk][i][j] = atof(col);
                }
			}
        (*line_no)++;
		// read B
        fscanf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nO-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->B[idxk][i][m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->B[idxk][i][m] = atof(col);
                }
			}
        (*line_no)++;
	} // for all k
}
