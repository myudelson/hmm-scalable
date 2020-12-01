/*
 
 Copyright (c) 2012-2017, Michael (Mikhail) Yudelson
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

#include "HMMProblemPiABGK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemPiABGK::HMMProblemPiABGK(struct param *param) {
    //    this->sizes = {param->nK, param->nK, param->nK};
    this->sizes[0] = param->nK;
    this->sizes[1] = param->nK;
    this->sizes[2] = param->nK;
    this->n_params = param->nK * 4 + param->nG * 4;
    init(param);
}

void HMMProblemPiABGK::init(struct param *param) {
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, (size_t)this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    
    NUMBER *a_PI, ** a_A, ** a_B;
    init3Params(a_PI, a_A, a_B, nS, nO);
    
    //
    // setup params
    //
	NPAR i, j, m, idx, offset;
	NUMBER sumPI = 0;
	NUMBER sumA[nS];
	NUMBER sumB[nS];
	for(i=0; i<nS; i++) {
		sumA[i] = 0;
		sumB[i] = 0;
	}
    
    // write default parameters first
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
		a_A[i][((nS)-1)]  = 1 - sumA[i];
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
    this->PIg = init2D<NUMBER>(nG, nS);
    this->Ag  = init3D<NUMBER>(nG, nS, nS);
    this->Bg  = init3D<NUMBER>(nG, nS, nO);
    NCAT x;
    for(x=0; x<nK; x++) {
        cpy1D<NUMBER>(a_PI, this->pi[x], nS);
        cpy2D<NUMBER>(a_A,  this->A[x],  nS, nS);
        cpy2D<NUMBER>(a_B,  this->B[x],  nS, nO);
    }
    // PIg start with "no-effect" params,PI[i] = 1/nS
    for(x=0; x<nG; x++) {
        for(i=0; i<nS; i++) {
            this->PIg[x][i] = (NUMBER)1/nS;
            for(j=0; j<nS; j++)
                this->Ag[x][i][j] = (NUMBER)1/nS;
            for(m=0; m<nS; m++)
                this->Bg[x][i][m] = (NUMBER)1/nO;
        }
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
}

HMMProblemPiABGK::~HMMProblemPiABGK() {
    destroy();
}

void HMMProblemPiABGK::destroy() {
	// destroy model data - only additional stuff
	free2D<NUMBER>(this->PIg, this->p->nG);
	free3D<NUMBER>(this->Ag, this->p->nG, this->p->nS);
	free3D<NUMBER>(this->Bg, this->p->nG, this->p->nS);
}// ~HMMProblemPiABGK

NUMBER** HMMProblemPiABGK::getPI() { // same as getPIk
	return this->pi;
}

NUMBER** HMMProblemPiABGK::getPIk() {
	return this->pi;
}

NUMBER** HMMProblemPiABGK::getPIg() {
	return this->PIg;
}

NUMBER*** HMMProblemPiABGK::getA() {
	return this->A;
}

NUMBER*** HMMProblemPiABGK::getAk() {
	return this->A;
}

NUMBER*** HMMProblemPiABGK::getAg() {
	return this->Ag;
}

NUMBER*** HMMProblemPiABGK::getBk() {
	return this->B;
}

NUMBER*** HMMProblemPiABGK::getBg() {
	return this->Bg;
}

NUMBER* HMMProblemPiABGK::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->pi[x];
}

NUMBER** HMMProblemPiABGK::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemPiABGK::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER* HMMProblemPiABGK::getPIk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->pi[x];
}

NUMBER* HMMProblemPiABGK::getPIg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing PI_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->PIg[x];
}

NUMBER** HMMProblemPiABGK::getAk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemPiABGK::getAg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing A_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->Ag[x];
}

NUMBER** HMMProblemPiABGK::getBk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER** HMMProblemPiABGK::getBg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing B_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->Bg[x];
}

NUMBER HMMProblemPiABGK::getPI(struct data* dt, NPAR i) {
    NUMBER p = this->pi[dt->k][i], q = this->PIg[dt->g][i];
    return pairing(p,q);
    //    return sigmoid( logit( this->pi[dt->k][i] ) + logit( this->PIg[dt->g][i] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiABGK::getA(struct data* dt, NPAR i, NPAR j) {
    NUMBER p = this->A[dt->k][i][j], q = this->Ag[dt->g][i][j];
    return pairing(p,q);
    //    return sigmoid( logit( this->A[dt->k][i][j] ) + logit( this->Ag[dt->k][i][j] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiABGK::getB(struct data* dt, NPAR i, NPAR m) {
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;
    NUMBER p = this->B[dt->k][i][m], q = this->Bg[dt->g][i][m];
    return pairing(p,q);
}

void HMMProblemPiABGK::setGradPI(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0, ndat = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    //    o = dt->obs[t];
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        o = this->p->dat_obs[ dt->ix[t] ];
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->pi[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( fb->pi[i] * (1-fb->pi[i]) );
            fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_PI, this->p, (NUMBER)ndat);
}

void HMMProblemPiABGK::setGradA (FitBit *fb){
    if(this->p->block_fitting[1]>0) return;
    NDAT t, ndat = 0;
    NPAR o, i, j;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        for(t=1; t<dt->n; t++) {
            o = this->p->dat_obs[ dt->ix[t] ];
            for(i=0; i<fb->nS; i++)
                for(j=0; j<fb->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = 1 / safe0num( fb->A[i][j] * (1-fb->A[i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * ((o<0)?1:fb->B[j][o]) * dt->alpha[t-1][i] / safe0num(dt->p_O_param) ;
                }
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_A, this->p, (NUMBER)ndat);
}

void HMMProblemPiABGK::setGradB (FitBit *fb){
    if(this->p->block_fitting[2]>0) return;
    NDAT t, ndat = 0;
    NPAR o=0, i;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        for(t=0; t<dt->n; t++) {
            o = this->p->dat_obs[ dt->ix[t] ];
            if(o<0)
                continue;
            for(i=0; i<fb->nS; i++) {
                combined = getB(dt,i,o);
                deriv_logit = 1 / safe0num( fb->B[i][o] * (1-fb->B[i][o]) );
                fb->gradB[i][o] -= combined * (1-combined) * deriv_logit * dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * ((o<0)?1:fb->B[i][o]));
            }
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_B, this->p, (NUMBER)ndat);
}


void HMMProblemPiABGK::toFile(const char *filename) {
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
		fprintf(fid,"Ag\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%12.10f%s",this->Ag[g][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"Bg\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%12.10f%s",this->Bg[g][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->pi[k][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"Ak\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%12.10f%s",this->A[k][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"Bk\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%12.10f%s",this->B[k][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblemPiABGK::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure==STRUCTURE_PIABgk)
        loglik_rmse[0] += GradientDescent();
    else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemPiABGK::GradientDescent() {
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
        NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK); // counting concecutive number of convergences for skills
        NCAT* iter_qual_group = Calloc(NCAT, (size_t)nG); // counting concecutive number of convergences for groups
        NCAT* iter_fail_skill = Calloc(NCAT, (size_t)nK); // counting concecutive number of failures to converge for skills
        NCAT* iter_fail_group = Calloc(NCAT, (size_t)nG); // counting concecutive number of failures to converge for groups
        int skip_k = 0, skip_g = 0;
        // utilize fitting larger data first
        
        int i = 0; // count runs
//        int parallel_now = this->p->parallel==1; //PAR
//        #pragma omp parallel if(parallel_now) shared(iter_qual_group,iter_qual_skill,iter_fail_skill,iter_fail_group)//PAR
//        {//PAR
            while(skip_k<nK || skip_g<nG) {
                //
                // Skills first
                //
                if(skip_k<nK) {
//                    #pragma omp for schedule(dynamic) //PAR
                    for(k=0; k<nK; k++) { // for all A,B-by-skill
                        if(iter_qual_skill[k]==iterations_to_qualify || iter_fail_skill[k]==iterations_to_qualify)
                            continue;
                        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                        // link and fit
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
                        
                        // count number of concecutive failures
                        if(fr.iter==this->p->maxiter) {
                            iter_fail_skill[k]++;
                        } else {
                            iter_fail_skill[k]=0;
                        }

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
//                                    #pragma omp critical(update_skip_k)//PAR
//                                    {//PAR
                                        skip_k++;
//                                    }//PAR
                                    if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                        printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                    }
                                }
                            } else
                                iter_qual_skill[k]=0;
                            // failed too many times in a row
                            if( iter_fail_skill[k]==iterations_to_qualify ) {
//                                #pragma omp critical(update_skip_k)//PAR
//                                {//PAR
                                skip_k++;
//                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                    printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        } // decide on convergence
                        delete fb;
                    } // for all skills
                }
                //
                // PIg second
                //
                if(skip_g<nG){
//                    #pragma omp for schedule(dynamic)//PAR
                    for(g=0; g<nG; g++) { // for all PI-by-user
                        if(iter_qual_group[g]==iterations_to_qualify || iter_fail_group[g]==iterations_to_qualify)
                            continue;
                        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                        // link
                        // vvvvvvvvvvvvvvvvvvvvv ONLY PART THAT IS DIFFERENT FROM HMMProblemPiGK
                        fb->link(this->getPIg(g), this->getAg(g), this->getBg(g), this->p->g_numk[g], this->p->g_k_data[g]);
                        if(this->p->block_fitting[0]!=0) fb->pi = NULL;
                        if(this->p->block_fitting[1]!=0) fb->A  = NULL;
                        if(this->p->block_fitting[2]!=0) fb->B  = NULL;
                        // ^^^^^^^^^^^^^^^^^^^^^

                        fb->Cslice = 1;
                        fb->init(FBS_PARm1);
                        fb->init(FBS_PARm2);
                        fb->init(FBS_GRAD);
                        if(this->p->solver==METHOD_CGD) {
                            fb->init(FBS_GRADm1);
                            fb->init(FBS_DIRm1);
                        }
                        FitResult fr = GradientDescentBit(fb);

                        // count number of concecutive failures
                        if(fr.iter==this->p->maxiter) {
                            iter_fail_group[g]++;
                        } else {
                            iter_fail_group[g]=0;
                        }

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
//                                    #pragma omp critical(update_skip_g)//PAR
//                                    {//PAR
                                        skip_g++;
//                                    }//PAR
                                    if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                        printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                                    }
                                }
                            }
                            else
                                iter_qual_group[g]=0;
                            // failed too many times in a row
                            if( iter_fail_group[g]==iterations_to_qualify ) {
//                                #pragma omp critical(update_skip_k)//PAR
//                                {//PAR
                                skip_g++;
//                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                    printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        } // decide on convergence
                        delete fb;
                    } // for all groups
                }
//                #pragma omp single//PAR
//                {//PAR
                    i++;
//                }//PAR
            }
//        }//PAR
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_qual_group != NULL) free(iter_qual_group);
        if( iter_fail_skill != NULL ) free(iter_fail_skill);
        if( iter_fail_group != NULL) free(iter_fail_group);
        
    } // if not "force single skill
    
    // compute loglik
    return getSumLogPOPara(this->p->nSeq, this->p->k_data);
}

void HMMProblemPiABGK::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
	NPAR i,j,m;
	NCAT k = 0, g = 0, idxk = 0, idxg = 0;
	string s;
    char col[2048];
    std::map<std::string,NCAT>::iterator it;
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
	// read grouped PIg, Ag, and Bg
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
		// read Ag
        fscanf(fid,"Ag\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++) {
                if(i==(this->p->nS-1) && j==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->Ag[idxg][i][j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->Ag[idxg][i][j] = atof(col);
                }
			}
		// read Bg
        fscanf(fid,"Bg\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nO-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->Bg[idxg][i][m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->Bg[idxg][i][m] = atof(col);
                }
			}
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
        fscanf(fid,"Ak\t");
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
