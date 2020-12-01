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
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK;//, nG = this->p->nG;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, (size_t)this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    
    NUMBER *a_PI, ** a_A, ** a_B;
    init3Params(a_PI, a_A, a_B, nS, nO);
    
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
//		if(this->is_multi[x]==0) {
			cpy1D<NUMBER>(a_PI, this->pi[x], nS);
			cpy2D<NUMBER>(a_A,  this->A[x],  nS, nS);
			cpy2D<NUMBER>(a_B,  this->B[x],  nS, nO);
//		} else {
//			for(i=0; i<nS; i++) {
//				this->pi[x][i] = 0.5;
//				for(j=0; j<nS; j++)
//					this->A[x][i][j] = 0.5;
//				for(m=0; m<nO; m++)
//					this->B[x][i][m] = 0.5;
//			}
//		}
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
	// set up parameter union
//	this->pu = new PUCorbettianAdditive();
	this->pu = new PULogistic();
}

HMMProblemComp::~HMMProblemComp() {
    destroy();
}

void HMMProblemComp::destroy() {
	// destroy additional model data
	free(this->is_multi);
	// delete parameter union
	free(this->pu);
}// ~HMMProblemComp

//NUMBER** HMMProblemComp::getPI() {
//	return this->pi;
//}
//
NUMBER* HMMProblemComp::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->pi[x];
}

NUMBER** HMMProblemComp::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemComp::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getPI(struct data* dt, NPAR i) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if param->tag1==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question
	if( this->p->tag1==0 || this->is_multi[dt->k] == 0) { // we ask for single skill specifically or this skill is never multiskill
		r = this->pi[dt->k][i];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->pi[dt->k][i];
		} else { // collect multiple \pi[i]
			NUMBER* q = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				q[l] = this->pi[k][i];
			}
			r = this->pu->unite(q, n_skills, NULL, 0); //squishing(q, n_skills);
			free(q);
		} // definitely multiskill
	} // end sometimes multiskill
    return r;
}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getA(struct data* dt, NPAR i, NPAR j) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if param->tag1==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question
	if( this->p->tag1==0 || this->is_multi[dt->k] == 0) { // we ask for single skill specifically or this skill is never multiskill
		r = this->A[dt->k][i][j];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->A[dt->k][i][j];
		} else { // collect multiple \pi[i]
			NUMBER* q = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				q[l] = this->A[k][i][j];
			}
			r = this->pu->unite(q, n_skills, NULL, 0); //r = squishing(q, n_skills);
			free(q);
		} // definitely multiskill
	} // end sometimes multiskill
	return r;
}

// Stateful, depends on whether the flag 'cnt' is on/off and the skill is_multi and the other skill parameters
NUMBER HMMProblemComp::getB(struct data* dt, NPAR i, NPAR m) {
	NUMBER r = 0;
	// dt->cnt==0 means, get dt->k'th skill only
	// if param->tag1==1, try and see whether there are other "neighbours"
	// dt->t should be set to the global index of the row in question

	// special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
	// in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
	if(m<0) {
		r = 1;
	} else if( this->p->tag1==0 || this->is_multi[dt->k] == 0) { // we ask for single skill specifically or this skill is never multiskill
		r = this->B[dt->k][i][m];
	} else {
		NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
		if( n_skills==1 ) { // this skill is alone on row dt->t
			r = this->B[dt->k][i][m];
		} else { // collect multiple \pi[i]
			NUMBER* q = Calloc(NUMBER, (size_t)n_skills );
			NCAT k=0;
			for(NCAT l=0; l<n_skills; l++) {
				k = this->p->dat_skill_stacked[ this->p->dat_skill_rix[ dt->t ]+l ];
				q[l] = this->B[k][i][m];
			}
			r = this->pu->unite(q, n_skills, NULL, 0); //r = squishing(q, n_skills);
			free(q);
		} // definitely multiskill
	} // end sometimes multiskill
	return r;
}

//
void HMMProblemComp::setGradPI(FitBit *fb){
	if(this->p->block_fitting[0]>0) return;
	NDAT t = 0, ndat = 0;
	NPAR i, o;
//	NUMBER combined, deriv_logit;
	struct data* dt;
//	NDAT tag1 = this->p->tag1; // use global tag1 setting
//	this->p->tag1 = 1;
//	this->p->tag1 = 0; // only get parameter values not united with "neighburs", the only difference
	for(NCAT x=0; x<fb->xndat; x++) {
		dt = fb->x_data[x];
		if( dt->cnt!=0 ) continue;
		ndat += dt->n;
		o = this->p->dat_obs[ dt->ix[t] ];//->get( dt->ix[t] );
		dt->t = dt->ix[t]; // statefullness
		for(i=0; i<this->p->nS; i++) {
			NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
			if( this->is_multi[dt->k] == 0 || n_skills == 1) {
				fb->gradPI[i] -= dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
			} else {
//				combined = getPI(dt,i);
//				deriv_logit = 1 / safe0num( fb->pi[i] * (1-fb->pi[i]) );
//				fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
				// compute with parameter unite
				fb->gradPI[i] -= this->pu->derivativeUnite(getPI(dt,i), fb->pi[i], 0.0/*NULL*/, 0 /*standard params*/) * dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
			}
		}
	}
	if( this->p->Cslices>0 ) // penalty
		fb->addL2Penalty(FBV_PI, this->p, (NUMBER)ndat);
//	this->p->tag1 = tag1;
}


void HMMProblemComp::setGradA (FitBit *fb){
	if(this->p->block_fitting[1]>0) return;
	NDAT t, ndat = 0;
	NPAR o, i, j;
//	NUMBER combined, deriv_logit;
	struct data* dt;
//	NDAT tag1 = this->p->tag1;  // use global tag1 setting
//	this->p->tag1 = 1;
//	this->p->tag1 = 0; // only get parameter values not united with "neighburs", the only difference
	for(NCAT x=0; x<fb->xndat; x++) {
		dt = fb->x_data[x];
		if( dt->cnt!=0 ) continue;
		ndat += dt->n;
		for(t=1; t<dt->n; t++) {
			o = this->p->dat_obs[ dt->ix[t] ];//->get( dt->ix[t] );
			dt->t = dt->ix[t]; // statefullness
			for(i=0; i<this->p->nS; i++) {
				for(j=0; j<this->p->nS; j++) {
					NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
					if( this->is_multi[dt->k] == 0 || n_skills == 1) {
						fb->gradA[i][j] -= dt->beta[t][j] * ((o<0)?1:fb->B[j][o]) * dt->alpha[t-1][i] / safe0num(dt->p_O_param);
					} else {
//						combined = getA(dt,i,j);
//						deriv_logit = 1 / safe0num( fb->A[i][j] * (1-fb->A[i][j]) );
//						fb->gradA[i][j] -= combined * safe0num(1-combined) * deriv_logit * dt->beta[t][j] * ((o<0)?1:fb->B[j][o]) * dt->alpha[t-1][i] / safe0num(dt->p_O_param);
						// compute with parameter unite
						fb->gradA[i][j] -= this->pu->derivativeUnite(getA(dt,i,j), fb->A[i][j], 0.0/*NULL*/, 0) * dt->beta[t][j] * ((o<0)?1:fb->B[j][o]) * dt->alpha[t-1][i] / safe0num(dt->p_O_param);
					}
				}
			}
		}
	}
	if( this->p->Cslices>0 ) // penalty
		fb->addL2Penalty(FBV_A, this->p, (NUMBER)ndat);
//	this->p->tag1 = tag1;
}


void HMMProblemComp::setGradB (FitBit *fb){
	if(this->p->block_fitting[2]>0) return;
	NDAT t, ndat = 0;
	NPAR o, o0, i;//, j;
	struct data* dt;
//	NUMBER combined, deriv_logit;
//	NDAT tag1 = this->p->tag1; // use global tag1 setting
//	this->p->tag1 = 1;
//	this->p->tag1 = 0; // only get parameter values not united with "neighburs", the only difference
	for(NCAT x=0; x<fb->xndat; x++) {
		dt = fb->x_data[x];
		if( dt->cnt!=0 ) continue;
		ndat += dt->n;
		for(t=0; t<dt->n; t++) { // Levinson MMFST
			dt->t = dt->ix[t]; // statefullness
			o  = this->p->dat_obs[ dt->ix[t] ];//->get( dt->ix[t] );
			o0 = this->p->dat_obs[ dt->ix[0] ];//->get( dt->ix[t] );
			if(o<0) // if no observation -- skip
				continue;
			for(i=0; i<this->p->nS; i++) {
				NCAT n_skills = this->p->dat_skill_rcount[ dt->t ];
				if( this->is_multi[dt->k] == 0 || n_skills == 1) {
					fb->gradB[i][o] -= dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * fb->B[i][o]);
				} else {
//					combined = getB(dt,i,o);
//					deriv_logit = 1 / safe0num( fb->B[i][o] * (1-fb->B[i][o]) );
//					fb->gradB[i][o] -= combined * (1-combined) * deriv_logit * dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * ((o<0)?1:fb->B[i][o]));
					// compute with parameter unite
					fb->gradB[i][o] -= this->pu->derivativeUnite(getB(dt,i,o), fb->B[i][o], 0.0/*NULL*/, 0) * dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * ((o<0)?1:fb->B[i][o]));
				}
			}
			
//			for(j=0; j<this->p->nS; j++)
//				if(t==0) {
//					fb->gradB[j][o] -= (o0==o) * fb->pi[j] * dt->beta[t][j];
//				} else {
//					for(i=0; i<this->p->nS; i++)
//						fb->gradB[j][o] -= ( dt->alpha[t-1][i] * fb->A[i][j] * dt->beta[t][j] /*+ (o0==o) * getPI(dt,j) * dt->beta[0][j]*/ ) / safe0num(dt->p_O_param); // Levinson MMFST
//				}
		}
	}
	if( this->p->Cslices>0 ) // penalty
		fb->addL2Penalty(FBV_B, this->p, (NUMBER)ndat);
//	this->p->tag1 = tag1;
}


void HMMProblemComp::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure == STRUCTURE_COMP) {
        loglik_rmse[0] += GradientDescent();
    } else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
	
	this->p->tag1 = 1; // for prediction needs after this
}


NUMBER HMMProblemComp::GradientDescent() {
	NCAT k, x;
    NCAT nK = this->p->nK;

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
    if( this->p->single_skill!=2 ) { // if not "force single skill"
		NCAT first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
		NCAT iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
		NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK); //
		NCAT* iter_fail_skill = Calloc(NCAT, (size_t)nK); // counting consecutive number of failures to converge for skills
		int skip_k = 0;
		NDAT sum_is_multi = 0;
		//
		// Main fit Stage 1. Single-skills, skills that only show up on single-skill rows.
		//

//		int parallel_now = this->p->parallel==1; //PAR
//		#pragma omp parallel if(parallel_now) shared(iter_qual_skill, iter_fail_skill)//PAR
//		{//PAR
//		#pragma omp for schedule(dynamic) //PAR
		for(x=0; x<nK; x++) {
//			#pragma omp critical(update_sum_is_multi)//PAR
//			{//PAR
			sum_is_multi += (this->is_multi[x] == 1);
//			}//PAR
			if(this->is_multi[x] == 0) {  // READ WARNING BELOW IF YOU COMMENT
				// vvvv Extract from HMMProblem.cpp
				NCAT xndat;
				struct data** x_data;
				xndat = this->p->k_numg[x];
				x_data = this->p->k_g_data[x];
                k = x_data[0]->k; // grab the first one, since all of them are the same k
				FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
				fb->link( this->getPI(x), this->getA(x), this->getB(x), xndat, x_data);
				if(this->p->block_fitting[0]!=0) fb->pi = NULL;
				if(this->p->block_fitting[1]!=0) fb->A  = NULL;
				if(this->p->block_fitting[2]!=0) fb->B  = NULL;
				
				FitResult fr;
				fb->init(FBS_PARm1);
				fb->init(FBS_GRAD);
				if(this->p->solver==METHOD_CGD) {
					fb->init(FBS_DIR);
					fb->init(FBS_DIRm1);
					fb->init(FBS_GRADm1);
				}
				if(this->p->solver==METHOD_GBB) {
					fb->init(FBS_GRADm1);
				}
				fb->init(FBS_PARm2); // do this for all in order to capture oscillation, e.g. if new param at time t is close to param at t-2 (tolerance)
				
//				fr = GradientDescentBit(fb);
				NDAT tag1 = this->p->tag1;   // FORCE
				this->p->tag1 = 1;           // FORCE
				fr = GradientDescentBit(fb); // FORCE
				this->p->tag1 = tag1;        // FORCE
			
				// count number of concecutive failures
				if(fr.iter==this->p->maxiter) {
					iter_fail_skill[k]++;
				} else {
					iter_fail_skill[k]=0;
				}

				delete fb;
				
				if( ( (fr.conv || fr.iter==this->p->maxiter) )) {
					if(!this->p->quiet)
						printf("skill %5d, seq %5d, dat %8d, iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", x, xndat, fr.ndat, fr.iter,fr.pO0,fr.pO,fr.conv);
				}
				// ^^^^
//				#pragma omp critical(update_skip_k)//PAR
//				{//PAR
				skip_k++;
//				}//PAR
			} // if(this->is_multi[x] == 0) // READ WARNING BELOW IF YOU COMMENT
		}
		/**/
//        }//PAR
		
		//
		// Main fit Stage 2. Multi-skills (skills that at least sometimes show up in multi-skil rows
		//
//		skip_k -= sum_is_multi; // update skills showing up in multi-skill rows // WATCHOUT
		// WARNING!! if you comment out `if(this->is_multi[x] == 0)` above, uncomment the previous line, otherwise the previous line should stay commented
        int i = 0; // count runs
		NDAT tag1 = this->p->tag1;
		this->p->tag1 = 1;
//        #pragma omp parallel if(parallel_now) shared(iter_qual_skill,iter_fail_skill)//PAR
//        {//PAR
        while(skip_k<nK || i < this->p->maxiter) { // add maxiter as cap
            //
            // Multi-skills, a version of what iBKT was doing
            //
            if(skip_k<nK) {
//                #pragma omp for schedule(dynamic) //PAR
                for(k=0; k<nK; k++) { // for all A,B-by-skill
                    if(this->is_multi[k]==0 || iter_qual_skill[k]==iterations_to_qualify || iter_fail_skill[k]==iterations_to_qualify)
                        continue;
                    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                    // link
                    fb->link( this->getPI(k), this->getA(k), this->getB(k), this->p->k_numg[k], this->p->k_g_data[k]);
                    if(this->p->block_fitting[0]!=0) fb->pi = NULL;
                    if(this->p->block_fitting[1]!=0) fb->A  = NULL;
                    if(this->p->block_fitting[2]!=0) fb->B  = NULL;

                    fb->Cslice = 0;
                    fb->init(FBS_PARm1);
                    fb->init(FBS_PARm2);
                    fb->init(FBS_GRAD);
                    if(this->p->solver==METHOD_CGD) {
						fb->init(FBS_DIR);
                        fb->init(FBS_GRADm1);
                        fb->init(FBS_DIRm1);
                    }
					if(this->p->solver==METHOD_GBB) {
						fb->init(FBS_GRADm1);
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
                        if(fr.iter==1 /*e<=this->p->tol*/ || fb->xndat==0) { // converged quick, or don't care (others all converged
                            iter_qual_skill[k]++;
                            if(fb->xndat==0) {
                                iter_qual_skill[k]=iterations_to_qualify;
                                fr.conv = 1;
                            }
                            if(iter_qual_skill[k]==iterations_to_qualify ) {// criterion met, or don't care (others all converged)
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
//            #pragma omp single//PAR
//            {//PAR
            i++;
//            }//PAR
        }
//        }//PAR
		this->p->tag1 = tag1;
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_fail_skill != NULL ) free(iter_fail_skill);

    } // if not "force single skill
        
    // compute loglik
    return getSumLogPOPara(this->p->nSeq, this->p->k_data);
}

void HMMProblemComp::toFile(const char *filename) {
	switch(this->p->structure)
	{
		case STRUCTURE_COMP:
			toFileSkill(filename);
			break;
			//        case STRUCTURE_GROUP:
			//            toFileGroup(filename);
			//            break;
		default:
			fprintf(stderr,"Solver specified is not supported.\n");
			break;
	}
}

void HMMProblemComp::toFileSkill(const char *filename) {
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
	
	NCAT k;
	std::map<NCAT,std::string>::iterator it;
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		NPAR i,j,m;
		fprintf(fid,"PI\t");
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

void HMMProblemComp::producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NDAT t) {
	NPAR m, i, l;
	NCAT k, g, nS = this->p->nS, nO = this->p->nO;
	NUMBER *a_L = init1D<NUMBER>((NDAT)nS);
	NUMBER **a_B = init2D<NUMBER>((NDAT)nS, (NDAT)nO);
	
	g = this->p->dat_group[t];
	NCAT *ks;
	int nks;
	if(this->p->multiskill==0) {
		fprintf(stderr,"WARNING! Multi-skill flag should have been set!");
		k = this->p->dat_skill[t];
		ks = &k;
		nks = 1;
	} else {
		NDAT rix = this->p->dat_skill_rix[t];
		k = this->p->dat_skill_stacked[ rix ];
		ks = &this->p->dat_skill_stacked[ rix ];
		nks = this->p->dat_skill_rcount[t];
	}
	
	if(nks==1) {
		for(i=0; i<nS; i++) {
			a_L[i] = group_skill_map[g][ ks[0] ][i];
			for(m=0; m<nO; m++)
				a_B[i][m] = this->B[ ks[0] ][i][m];
		}
	} else {
		NUMBER *q = NULL;
		for(i=0; i<nS; i++) {
			q = Calloc(NUMBER, (size_t)nks);
			for(l=0;l<nks;l++) q[l] = group_skill_map[g][ ks[l] ][i];
			a_L[i] = this->pu->unite(q, nks, NULL, 0); //a_L[i] = squishing(q, nks);
			free(q);

			for(m=0; m<nO; m++) {
				q = Calloc(NUMBER, (size_t)nks);
				for(l=0;l<nks;l++) q[l] = this->B[ ks[l] ][i][m];
				free(q);
				a_B[i][m] = this->pu->unite(q, nks, NULL, 0); //a_B[i][m] = squishing(q, nks);
			}
		}
	}

	for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
	for(m=0; m<this->p->nO; m++)
		for(i=0; i<this->p->nS; i++)
			local_pred[m] += a_L[i] * a_B[i][m];
	
	free(a_L);
	free2D<NUMBER>(a_B, nS);
}


