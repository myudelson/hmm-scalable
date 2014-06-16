//
//  HMMProblemPiGKPloGK.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 8/31/12.
//
//

#include "HMMProblemPiGK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemPiGK::HMMProblemPiGK(struct param *param) {
    for(NPAR i=0; i<3; i++) this->sizes[i] = param->nK;
//    this->sizes = {param->nK, param->nK, param->nK};
    this->n_params = param->nK * 4 + param->nG;
    init(param);
}

void HMMProblemPiGK::init(struct param *param) {
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
	if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
		this->pi  = init2D<NUMBER>(nK, nS);
		this->A   = init3D<NUMBER>(nK, nS, nS);
		this->B   = init3D<NUMBER>(nK, nS, nO);
		this->PIg = init2D<NUMBER>(nG, nS);
        NCAT x;
		for(x=0; x<this->p->nK; x++) {
			cpy1D<NUMBER>(a_PI, this->pi[x], nS);
			cpy2D<NUMBER>(a_A,  this->A[x],  nS, nS);
			cpy2D<NUMBER>(a_B,  this->B[x],  nS, nO);
        }
        // PIg start with same params
		for(x=0; x<this->p->nG; x++)
			cpy1D<NUMBER>(a_PI, this->PIg[x], this->p->nS);
	} else {
		fprintf(stderr,"params do not meet constraints.\n");
		exit(1);
	}
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
	for(i=0; i<this->p->nS; i++) {
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
}

HMMProblemPiGK::~HMMProblemPiGK() {
    destroy();
}

void HMMProblemPiGK::destroy() {
	// destroy additional model data
	free2D<NUMBER>(this->PIg, this->p->nG);
}// ~HMMProblemPiGK

//NUMBER** HMMProblemPiGK::getPI() { // same as getPIk
//	return this->PI;
//}

NUMBER** HMMProblemPiGK::getPIk() {
	return this->pi;
}

NUMBER** HMMProblemPiGK::getPIg() {
	return this->PIg;
}

//NUMBER* HMMProblemPiGK::getPI(NCAT x) { // same as getPIk(x)
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->PI[x];
//}
//
//NUMBER** HMMProblemPiGK::getA(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->A[x];
//}
//
//NUMBER** HMMProblemPiGK::getB(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->B[x];
//}

NUMBER* HMMProblemPiGK::getPIk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->pi[x];
}

NUMBER* HMMProblemPiGK::getPIg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing PI_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->PIg[x];
}

NUMBER HMMProblemPiGK::getPI(struct data* dt, NPAR i) {
    NUMBER p = this->pi[dt->k][i], q = this->PIg[dt->g][i];
    return 1/( 1 + (1-p)*(1-q)/(p*q) );

//    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
//    NUMBER item = this->p->item_complexity[ this->p->dat_item->get( dt->ix[0] ) ];
//    NUMBER v = 1/( 1 + (1-p)*(1-q)*(1-item)/(p*q*item) );
//    return v;
    
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiGK::getA(struct data* dt, NPAR i, NPAR j) {
    return this->A[dt->k][i][j];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiGK::getB(struct data* dt, NPAR i, NPAR m) {
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;
    return this->B[dt->k][i][m];
}

void HMMProblemPiGK::setGradPI(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
//    o = dt->obs[t];
        o = this->p->dat_obs[ dt->ix[t] ];//->get( dt->ix[t] );
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( fb->pi[i] * (1-fb->pi[i]) );
			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
        }
        // penalty
        for(i=0; i<fb->nS && this->p->C!=0; i++)
            fb->gradPI[i] += L2penalty(this->p,fb->pi[i], 0.5);
    }
}

void HMMProblemPiGK::toFile(const char *filename) {
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    
    // write solved id
    writeSolverInfo(fid, this->p);
    
	fprintf(fid,"Null skill ratios\t");
	for(NPAR m=0; m<this->p->nO; m++)
		fprintf(fid," %10.7f%s",this->null_obs_ratio[m],(m==(this->p->nO-1))?"\n":"\t");
	NCAT k, g;
    NPAR i,j,m;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<this->p->nG; g++) {
		it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		fprintf(fid,"PIg\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PIg[g][i],(i==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->pi[k][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%10.8f%s",this->A[k][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%10.8f%s",this->B[k][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblemPiGK::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure == STRUCTURE_PIgk) {
        loglik_rmse[0] += GradientDescent();
    } else {
        fprintf(stderr,"Solver specified is not supported.\n");
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemPiGK::GradientDescent() {
	NCAT k, g, /*ki, gi, nX, */x;
    NCAT nK = this->p->nK, nG = this->p->nG;
//    NUMBER ll;
    
//    NUMBER crit;
//    NCAT q;
//    NCAT nQ = this->p->nK + this->p->nG;
//    FitBit **fbs = Calloc(FitBit*, nQ);
//    for(q=0;q<nQ;q++) {
//        NCAT xndat;
//        struct data** x_data;
//        if(q<this->p->nK) { // skills
//            x = q;
//            xndat = this->p->k_numg[x];
//            x_data = this->p->k_g_data[x];
//            fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
//            fbs[q]->link( HMMProblem::getPI(x), HMMProblem::getA(x), HMMProblem::getB(x), xndat, x_data);
//        } else { // groups
//            x = q - this->p->nK;
//            xndat = this->p->g_numk[x];
//            x_data = this->p->g_k_data[x];
//            fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
//            fbs[q]->link( this->getPIg(x), NULL, NULL, xndat, x_data);
//        }
//        fbs[q]->init(FBS_PARm1);
//        fbs[q]->init(FBS_GRAD);
//        if(this->p->solver==METHOD_CGD) {
//            fbs[q]->init(FBS_GRADm1);
//            fbs[q]->init(FBS_DIRm1);
//        }
//    }
    
    
        
	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill>0) {
        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
        fb->init(FBS_PARm1);
        fb->init(FBS_GRAD);
        if(this->p->solver==METHOD_CGD) {
            fb->init(FBS_GRADm1);
            fb->init(FBS_DIRm1);
        }
        fb->link( HMMProblem::getPI(0), HMMProblem::getA(0), HMMProblem::getB(0), this->p->nSeq, this->p->k_data);// link skill 0 (we'll copy fit parameters to others
        NCAT* original_ks = Calloc(NCAT, (size_t)this->p->nSeq);
        for(x=0; x<this->p->nSeq; x++) { original_ks[x] = this->p->all_data[x].k; this->p->all_data[x].k = 0; } // save progonal k's
        FitResult fr = GradientDescentBit(fb);
        for(x=0; x<this->p->nSeq; x++) { this->p->all_data[x].k = original_ks[x]; } // restore original k's
        free(original_ks);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
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
        int skip_k = 0, skip_g = 0;
        // utilize fitting larger data first
//        NDAT newnK=nK, newnG=nG;
//        this->reorderSequences(&newnK, &newnG, false);

        int i = 0; // count runs
//        while(skip_k<newnK || skip_g<newnG) {
        int parallel_now = this->p->parallel==1; //PAR
        #pragma omp parallel if(parallel_now) shared(iter_qual_group,iter_qual_skill)//PAR
        {//PAR
        while(skip_k<nK || skip_g<nG) {
            //
            // Skills first
            //
//            for(ki=0; ki<newnK && skip_k<newnK; ki++) { // for all A,B-by-skill

            
            if(skip_k<nK) {
                #pragma omp for schedule(dynamic) //PAR
                for(k=0; k<nK; k++) { // for all A,B-by-skill
                    if(iter_qual_skill[k]==iterations_to_qualify)
                        continue;
                    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
                    fb->init(FBS_PARm1);
                    fb->init(FBS_GRAD);
                    if(this->p->solver==METHOD_CGD) {
                        fb->init(FBS_GRADm1);
                        fb->init(FBS_DIRm1);
                    }
    //                k = sortstrip_k[ki].id; // grab reordered k
                    
    //                NCAT xndat = this->p->k_numg[k];
    //                struct data** x_data = this->p->k_g_data[k];
                    // link and fit
                    fb->link( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), this->p->k_numg[k], this->p->k_g_data[k]);// link skill 0 (we'll copy fit parameters to others
                    FitResult fr = GradientDescentBit(fb);
    //                fbs[k]->link( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), this->p->k_numg[k], this->p->k_g_data[k]);// link skill 0 (we'll copy fit parameters to others
    //                fr = GradientDescentBit(fbs[k], false /*is1SkillForAll*/);
                    // decide on convergence
                    if(i>=first_iteration_qualify) {
    //                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==newnG) { // converged quick, or don't care (others all converged
                        if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==nG) { // converged quick, or don't care (others all converged
                            iter_qual_skill[k]++;
    //                        if(iter_qual_skill[k]==iterations_to_qualify || skip_g==newnG) {// criterion met, or don't care (others all converged)
                            if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
    //                            if(skip_g==newnG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                                if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                                #pragma omp critical(update_skip_k)//PAR
                                {//PAR
                                    skip_k++;
                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
    //                                computeAlphaAndPOParam(fb->xndat, fb->x_data);
    //                                computeAlphaAndPOParam(fbs[k]->xndat, fbs[k]->x_data);
                                    printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        }
                        else
                            iter_qual_skill[k]=0;
                    } // decide on convergence
    //                //
    //                // make copies of parameters to do gradient descend, link to copies, but grads computed from actual params
    //                //
    //                if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
    //                    swap1D<NUMBER>(cpyPI[k], this->getPI(k), this->p->nS);
    //                    swap2D<NUMBER>(cpyA[k],   this->getA(k), this->p->nS, this->p->nS);
    //                    swap2D<NUMBER>(cpyB[k],   this->getB(k), this->p->nS, this->p->nO);
    //                }
                    delete fb;
                } // for all skills
            }
            //
            // PIg second
            //
//            int z = 0;
//            for(gi=0; gi<newnG && skip_g<newnG; gi++) { // for all PI-by-user
            
            if(skip_g<nG){
                #pragma omp for schedule(dynamic)//PAR
                for(g=0; g<nG; g++) { // for all PI-by-user
                    if(iter_qual_group[g]==iterations_to_qualify)
                        continue;
                    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
                    fb->init(FBS_PARm1);
                    fb->init(FBS_GRAD);
                    if(this->p->solver==METHOD_CGD) {
                        fb->init(FBS_GRADm1);
                        fb->init(FBS_DIRm1);
                    }
    //                g = sortstrip_g[gi].id; // grab reordered g

    //                NCAT xndat = this->p->g_numk[g];
    //                struct data** x_data = this->p->g_k_data[g];
                    fb->link(this->getPIg(g), NULL, NULL, this->p->g_numk[g], this->p->g_k_data[g]);
                    FitResult fr = GradientDescentBit(fb);
    //                fbs[nK+g]->link(this->getPIg(g), NULL, NULL, this->p->g_numk[g], this->p->g_k_data[g]);
    //                fr = GradientDescentBit(fbs[nK+g], false /*is1SkillForAll*/);
    //                if(i>=first_iteration_qualify || fbs[nK+g]->xndat==0) { //can qualify or  student had no skill labelled rows
    //                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==newnK || fbs[nK+g]->xndat==0) { // converged quick, or don't care (others all converged), or  student had no skill labelled rows
                    if(i>=first_iteration_qualify || fb->xndat==0) { //can qualify or  student had no skill labelled rows
    //                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==newnK || fb->xndat==0) { // converged quick, or don't care (others all converged), or                          iter_qual_group[g]++;
                        if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==nK || fb->xndat==0) { // converged quick, or don't care (others all converged), or
                            if(fb->xndat==0) {
    //                        if(fbs[nK+g]->xndat==0) {
                                iter_qual_group[g]=iterations_to_qualify;
                                fr.conv = 1;
                            }
    //                        if(iter_qual_group[g]==iterations_to_qualify || skip_k==newnK) {// criterion met, or don't care (others all converged)
    //                            if(skip_k==newnK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                        if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                            if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                                #pragma omp critical(update_skip_g)//PAR
                                {//PAR
                                    skip_g++;
                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
    //                                computeAlphaAndPOParam(fb->xndat, fb->x_data);
    //                                computeAlphaAndPOParam(fbs[nK+g]->xndat, fbs[nK+g]->x_data);
                                    printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        }
                        else
                            iter_qual_group[g]=0;
                    } // decide on convergence
    //                //
    //                // make copies of parameters to do gradient descend, link to copies, but grads computed from actual params
    //                //
    //                if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
    //                    swap1D<NUMBER>(this->PIg[g], cpyPIg[g], this->p->nS);
    //                }
                    delete fb;
                } // for all groups
            }
//            fr.conv = checkConvergenceBig(fbs, nK+nG, this->p->tol, &crit);
//            
//            computeAlphaAndPOParam(this->p->nSeq, this->p->k_data);
//            ll = HMMProblem::getSumLogPOPara(this->p->nSeq, this->p->k_data);
//            printf("*%i ll=%15.7f, crit=%10.7f\n",i,ll,crit);
            #pragma omp single//PAR
            {//PAR
            i++;
            }//PAR
        }
        }//PAR
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_qual_group != NULL) free(iter_qual_group);

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
                fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
                fbs[q]->link( HMMProblem::getPI(x), HMMProblem::getA(x), HMMProblem::getB(x), xndat, x_data);
            } else { // groups
                x = q - this->p->nK;
                xndat = this->p->g_numk[x];
                x_data = this->p->g_k_data[x];
                fbs[q] = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
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


struct data***  HMMProblemPiGK::getExdendedData(NCAT xndat, struct data** x_data, NPAR kg_flag, NCAT* xxndat) {
    NCAT q, r, g, k;
    struct data*** xx_data;
    *xxndat = 0;
    NCAT yndat;
    struct data** z_data;
    struct data** y_data;
    if(kg_flag==0) {
        for(q=0; q<xndat; q++)
            *xxndat += this->p->g_numk[ x_data[q]->g ];
        
    } else if (kg_flag==1) {
        for(q=0; q<xndat; q++)
            *xxndat += this->p->k_numg[ x_data[q]->k ];
        
    } else {
        fprintf(stderr,"Unexpected value of kg_flag\n");
        return NULL;
    }
    z_data = (Malloc(struct data *, (size_t)*xxndat));
    xx_data = &z_data;
    NDAT i = 0;
    if(kg_flag==0) {
        for(q=0; q<xndat; q++) {
            g = x_data[q]->g;
            yndat = this->p->g_numk[g];
            y_data = this->p->g_k_data[g];
            for(r=0; r<yndat; r++)
                z_data[i++] = y_data[r];
        }
        
    } else if (kg_flag==1) {
        for(q=0; q<xndat; q++) {
            k = x_data[q]->k;
            yndat = this->p->k_numg[k];
            y_data = this->p->k_g_data[k];
            for(r=0; r<yndat; r++)
                z_data[i++] = y_data[r];
        }
    }
    return xx_data;
}

//NUMBER HMMProblemPiGK::GradientDescentX() {
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
//            printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
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
//                                printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
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
//                                printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
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
void HMMProblemPiGK::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
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
			for(m=0; m<this->p->nS; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nS-1)) {
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
