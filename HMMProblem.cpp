/*
 *  HMMProblem.cpp
 *  HMM
 *
 *  Created by Mikhail Yudelson on 5/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include "FitBit.h"
#include <math.h>
#include "HMMProblem.h"
#include <map>

//#include "liblinear/linear.h"

HMMProblem::HMMProblem() {
}

HMMProblem::HMMProblem(struct param *param) {
    NPAR i;
    switch (param->structure) {
        case STRUCTURE_SKILL: // Expectation Maximization (Baum-Welch)
            //            this->sizes[3] = {param->nK, param->nK, param->nK};
            for(i=0; i<3; i++) this->sizes[i] = param->nK;
            this->n_params = param->nK * 4;
            break;
        case STRUCTURE_GROUP: // Gradient Descent by group
            //            this->sizes = {param->nG, param->nG, param->nG};
            for(i=0; i<3; i++) this->sizes[i] = param->nG;
            this->n_params = param->nG * 4;
            break;
        default:
            fprintf(stderr,"Structure specified is not supported and should have been caught earlier\n");
            break;
    }
    init(param);
}

void HMMProblem::init(struct param *param) {
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, (size_t)this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    if( this->p->solver == METHOD_CGD && this->p->solver_setting == -1)
        this->p->solver_setting = 1; // default Fletcher-Reeves
    
    NPAR nS = this->p->nS, nO = this->p->nO;
    NUMBER *a_PI, ** a_A, ** a_B;
    init3Params(a_PI, a_A, a_B, nS, nO);
    
    //
    // setup params
    //
	NPAR i, j, idx, offset;
	NUMBER sumPI = 0;
	NUMBER sumA[this->p->nS];
	NUMBER sumB[this->p->nS];
	for(i=0; i<this->p->nS; i++) {
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
	offset = nS-1;
	for(i=0; i<nS; i++) {
		for(j=0; j<((nS)-1); j++) {
			idx = offset + i*((nS)-1) + j;
			a_A[i][j] = this->p->init_params[idx];
			sumA[i]  += this->p->init_params[idx];
		}
		a_A[i][((nS)-1)]  = 1 - sumA[i];
	}
	// polupale B
	offset = (nS-1) + nS*(nS-1);
	for(i=0; i<nS; i++) {
		for(j=0; j<((nO)-1); j++) {
			idx = offset + i*((nO)-1) + j;
			a_B[i][j] = this->p->init_params[idx];
			sumB[i] += this->p->init_params[idx];
		}
		a_B[i][((nO)-1)]  = 1 - sumB[i];
	}
    
    // mass produce PI's, A's, B's
	if( checkPIABConstraints(a_PI, a_A, a_B) ) {
		this->PI = init2D<NUMBER>((NDAT)this->sizes[0], (NDAT)nS);
		this->A =  init3D<NUMBER>((NDAT)this->sizes[1], (NDAT)nS, (NDAT)nS);
		this->B =  init3D<NUMBER>((NDAT)this->sizes[2], (NDAT)nS, (NDAT)nO);
        NCAT x;
		for(x=0; x<this->sizes[0]; x++)
			cpy1D<NUMBER>(a_PI, this->PI[x], (NDAT)nS);
		for(x=0; x<this->sizes[1]; x++)
			cpy2D<NUMBER>(a_A, this->A[x], (NDAT)nS, (NDAT)nS);
		for(x=0; x<this->sizes[2]; x++)
			cpy2D<NUMBER>(a_B, this->B[x], (NDAT)nS, (NDAT)nO);
	} else {
		fprintf(stderr,"params do not meet constraints.\n");
		exit(1);
	}
    // destroy setup params
	free(a_PI);
	free2D<NUMBER>(a_A, (NDAT)nS);
	free2D<NUMBER>(a_B, (NDAT)nS);
    
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
			idx = offset + i*nS + j;
			lbA[i][j] = this->p->param_lo[idx];
			ubA[i][j] = this->p->param_hi[idx];
		}
	// *B
	offset = nS + nS*nS;
	for(i=0; i<nS; i++)
		for(j=0; j<nO; j++) {
			idx = offset + i*nS + j;
			lbB[i][j] = this->p->param_lo[idx];
			ubB[i][j] = this->p->param_hi[idx];
		}

}

HMMProblem::~HMMProblem() {
    destroy();
}

void HMMProblem::destroy() {
	// destroy model data
    if(this->null_obs_ratio != NULL) free(this->null_obs_ratio);
	if(this->PI != NULL) free2D<NUMBER>(this->PI, this->sizes[0]);
	if(this->A  != NULL) free3D<NUMBER>(this->A,  this->sizes[1], this->p->nS);
	if(this->B  != NULL) free3D<NUMBER>(this->B,  this->sizes[2], this->p->nS);
	if(this->lbPI!=NULL) free(this->lbPI);
	if(this->ubPI!=NULL) free(this->ubPI);
	if(this->lbA!=NULL) free2D<NUMBER>(this->lbA, this->p->nS);
	if(this->ubA!=NULL) free2D<NUMBER>(this->ubA, this->p->nS);
	if(this->lbB!=NULL) free2D<NUMBER>(this->lbB, this->p->nS);
	if(this->ubB!=NULL) free2D<NUMBER>(this->ubB, this->p->nS);
}// ~HMMProblem

bool HMMProblem::hasNon01Constraints() {
	return this->non01constraints;
}

NUMBER** HMMProblem::getPI() {
	return this->PI;
}

NUMBER*** HMMProblem::getA() {
	return this->A;
}

NUMBER*** HMMProblem::getB() {
	return this->B;
}

NUMBER* HMMProblem::getPI(NCAT x) {
	if( x > (this->sizes[0]-1) ) {
		fprintf(stderr,"While accessing PI, skill index %d exceeded last index of the data %d.\n", x, this->sizes[0]-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER** HMMProblem::getA(NCAT x) {
	if( x > (this->sizes[1]-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->sizes[1]-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblem::getB(NCAT x) {
	if( x > (this->sizes[2]-1) ) {
		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->sizes[2]-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER* HMMProblem::getLbPI() {
	if( !this->non01constraints ) return NULL;
	return this->lbPI;
}

NUMBER** HMMProblem::getLbA() {
	if( !this->non01constraints ) return NULL;
	return this->lbA;
}

NUMBER** HMMProblem::getLbB() {
	if( !this->non01constraints ) return NULL;
	return this->lbB;
}

NUMBER* HMMProblem::getUbPI() {
	if( !this->non01constraints ) return NULL;
	return this->ubPI;
}

NUMBER** HMMProblem::getUbA() {
	if( !this->non01constraints ) return NULL;
	return this->ubA;
}

NUMBER** HMMProblem::getUbB() {
	if( !this->non01constraints ) return NULL;
	return this->ubB;
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblem::getPI(struct data* dt, NPAR i) {
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->PI[dt->k][i];
            break;
        case STRUCTURE_GROUP:
            return this->PI[dt->g][i];
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    }
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblem::getA (struct data* dt, NPAR i, NPAR j) {
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->A[dt->k][i][j];
            break;
        case STRUCTURE_GROUP:
            return this->A[dt->g][i][j];
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    }
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblem::getB (struct data* dt, NPAR i, NPAR m) {
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->B[dt->k][i][m];
            break;
        case STRUCTURE_GROUP:
            return this->B[dt->g][i][m];
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    }
}

bool HMMProblem::checkPIABConstraints(NUMBER* a_PI, NUMBER** a_A, NUMBER** a_B) {
	NPAR i, j;
	// check values
	NUMBER sum_pi = 0.0;
	NUMBER sum_a_row[this->p->nS];
	NUMBER sum_b_row[this->p->nS];
	for(i=0; i<this->p->nS; i++) {
		sum_a_row[i] = 0.0;
		sum_b_row[i] = 0.0;
	}
	
	for(i=0; i<this->p->nS; i++) {
		if( a_PI[i]>1.0 || a_PI[i]<0.0)
			return false;
		sum_pi += a_PI[i];
		for(j=0; j<this->p->nS; j++) {
			if( a_A[i][j]>1.0 || a_A[i][j]<0.0)
				return false;
			sum_a_row[i] += a_A[i][j];
		}// all states 2
		for(int m=0; m<this->p->nO; m++) {
			if( a_B[i][m]>1.0 || a_B[i][m]<0.0)
				return false;
			sum_b_row[i] += a_B[i][m];
		}// all observations
	}// all states
	if(sum_pi!=1.0)
		return false;
	for(i=0; i<this->p->nS; i++)
		if( sum_a_row[i]!=1.0 || sum_b_row[i]!=1.0)
			return false;
	return true;
}

NUMBER HMMProblem::getSumLogPOPara(NCAT xndat, struct data** x_data) {
	NUMBER result = 0.0;
	for(NCAT x=0; x<xndat; x++) result += (x_data[x]->cnt==0)?x_data[x]->loglik:0;
	return result;
}

void HMMProblem::initAlpha(NCAT xndat, struct data** x_data) {
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->alpha == NULL ) {
			x_data[x]->alpha = Calloc(NUMBER*, (size_t)x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->alpha[t] = Calloc(NUMBER, (size_t)nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->alpha[t][i] = 0.0;
		}
		// p_O_param
		x_data[x]->p_O_param = 0.0;
		x_data[x]->loglik = 0.0;
        // c - scaling
		if( x_data[x]->c == NULL ) {
            x_data[x]->c = Calloc(NUMBER, (size_t)x_data[x]->n);
        } else {
			for(t=0; t<x_data[x]->n; t++)
                x_data[x]->c[t] = 0.0;
        }
	} // for all groups in skill
}

void HMMProblem::initGamma(NCAT xndat, struct data** x_data) {
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->gamma == NULL ) {
			x_data[x]->gamma = Calloc(NUMBER*, (size_t)x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->gamma[t] = Calloc(NUMBER, (size_t)nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->gamma[t][i] = 0.0;
		}
	} // for all groups in skill
}

void HMMProblem::initXi(NCAT xndat, struct data** x_data) {
	NCAT x;
	NDAT t;
	NPAR i,j, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->xi == NULL ) {
			x_data[x]->xi = Calloc(NUMBER**, (size_t)x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++) {
				x_data[x]->xi[t] = Calloc(NUMBER*, (size_t)nS);
				for(i=0; i<nS; i++)
					x_data[x]->xi[t][i] = Calloc(NUMBER, (size_t)nS);
			}
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					for(j=0; j<nS; j++)
						x_data[x]->xi[t][i][j] = 0.0;
		}
	} // for all groups in skill
}

void HMMProblem::initBeta(NCAT xndat, struct data** x_data) {
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// beta
		if( x_data[x]->beta == NULL ) {
			x_data[x]->beta = Calloc(NUMBER*, (size_t)x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->beta[t] = Calloc(NUMBER, (size_t)nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->beta[t][i] = 0.0;
		}
	} // for all groups in skill
} // initBeta

void HMMProblem::computeAlphaAndPOParam(NCAT xndat, struct data** x_data) {
	NCAT x;
	NDAT t;
	NPAR i, j, o, nS = this->p->nS;
    
    //    NUMBER mult_c, old_pOparam, neg_sum_log_c;
	initAlpha(xndat, x_data);
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
        //        mult_c = 1;
        //        neg_sum_log_c = 0.0;
        x_data[x]->p_O_param = 0; // 0 for non-scaled
		for(t=0; t<x_data[x]->n; t++) {
            //			o = x_data[x]->obs[t];
			o = this->p->dat_obs->get( x_data[x]->ix[t] );
			if(t==0) { // it's alpha(1,i)
                // compute \alpha_1(i) = \pi_i b_i(o_1)
				for(i=0; i<nS; i++) {
					x_data[x]->alpha[t][i] = getPI(x_data[x],i) * ((o<0)?1:getB(x_data[x],i,o)); // if observatiob unknown use 1
                    x_data[x]->c[t] += x_data[x]->alpha[t][i];
                }
			} else { // it's alpha(t,i)
				// compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
				for(i=0; i<nS; i++) {
					for(j=0; j<nS; j++) {
						x_data[x]->alpha[t][i] += x_data[x]->alpha[t-1][j] * getA(x_data[x],j,i);
					}
					x_data[x]->alpha[t][i] *= ((o<0)?1:getB(x_data[x],i,o)); // if observatiob unknown use 1
                    x_data[x]->c[t] += x_data[x]->alpha[t][i];
				}
			}
            // scale \alpha_{t}(i) - same for t=1 or otherwise
            x_data[x]->c[t] = 1/x_data[x]->c[t];//safe0num();
			// compute elements of p(O|param) as a sum of alpha's of last observations in sequences
			if( t==(x_data[x]->n-1) ) {
				for(i=0; i<nS; i++)
					x_data[x]->p_O_param += x_data[x]->alpha[t][i];
                x_data[x]->loglik = -safelog(x_data[x]->p_O_param);
			}
		} // for all observations within skill-group
	} // for all groups in skill
}

void HMMProblem::computeBeta(NCAT xndat, struct data** x_data) {
	NCAT x;
	int t;
	NPAR i, j, o, nS = this->p->nS;
	initBeta(xndat, x_data);
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		for(t=(NDAT)(x_data[x]->n)-1; t>=0; t--) {
			if( t==(x_data[x]->n-1) ) { // last \beta
				// \beta_T(i) = 1
				for(i=0; i<nS; i++)
					x_data[x]->beta[t][i] = 1;// x_data[x]->c[t]; // was 1
			} else {
				// \beta_t(i) = \sum_{j=1}^N{beta_{t+1}(j) a_{ij} b_j(o_{t+1})}
                //				o = x_data[x]->obs[t+1]; // next observation
                o = this->p->dat_obs->get( x_data[x]->ix[t+1] );
				for(i=0; i<nS; i++) {
					for(j=0; j<nS; j++)
						x_data[x]->beta[t][i] += x_data[x]->beta[t+1][j] * getA(x_data[x],i,j) * ((o<0)?1:getB(x_data[x],j,o)); // if observatiob unknown use 1
                    // scale
                    //                    x_data[x]->beta[t][i] *= x_data[x]->c[t];
                }
			}
		} // for all observations, starting with last one
	} // for all groups within skill
}

void HMMProblem::computeXi(NCAT xndat, struct data** x_data){
	HMMProblem::initXi(xndat, x_data);
	NCAT x;
	NDAT t;
	NPAR i, j, o_tp1, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		for(t=0; t<(x_data[x]->n-1); t++) { // -1 is important
            //			o_tp1 = x_data[x]->obs[t+1];
            o_tp1 = this->p->dat_obs->get( x_data[x]->ix[t+1] );
            //            NUMBER denom = 0;
            //            for(i=0; i<a_nS; i++)
            //				for(j=0; j<a_nS; j++) {
            //                    denom += x_data[x]->alpha[t][i] * a_A[i][j] * x_data[x]->beta[t+1][j] * a_B[j][o_tp1];
            //                }
            //            for(i=0; i<a_nS; i++)
            //				for(j=0; j<a_nS; j++) {
            //                    x_data[x]->xi[t][i][j] = x_data[x]->alpha[t][i] * a_A[i][j] * x_data[x]->beta[t+1][j] * a_B[j][o_tp1] / safe0num(denom);
            //                    if( (x_data[x]->xi[t][i][j])!=(x_data[x]->xi[t][i][j]) ) {
            //                        int z = 0;
            //                    }
            //                }
            
			for(i=0; i<nS; i++)
				for(j=0; j<nS; j++) {
					x_data[x]->xi[t][i][j] = x_data[x]->alpha[t][i] * getA(x_data[x],i,j) * x_data[x]->beta[t+1][j] * getB(x_data[x],j,o_tp1) / safe0num(x_data[x]->p_O_param);
                }
		} // for all observations within skill-group
	} // for all groups in skill
}

void HMMProblem::computeGamma(NCAT xndat, struct data** x_data) {
	HMMProblem::initGamma(xndat, x_data);
	NCAT x;
	NDAT t;
	NPAR i, j, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		for(t=0; t<(x_data[x]->n-1); t++) { // -1 is important
			for(i=0; i<nS; i++)
				for(j=0; j<nS; j++) {
					x_data[x]->gamma[t][i] += x_data[x]->xi[t][i][j];
                }
		} // for all observations within skill-group
	} // for all groups in skill
}

void HMMProblem::setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0;
    NPAR i, o;
    //    if(kg_flag == 0) { // k THIS PARAM DOESN'T MATTER HERE
    //        o = dt->obs[t];
    o = this->p->dat_obs->get( dt->ix[t] );
    for(i=0; i<this->p->nS; i++) {
        fb->gradPI[i] -= dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,getPI(dt,i)); // PENALTY
    }
    //    }
}

void HMMProblem::setGradA (struct data* dt, FitBit *fb, NPAR kg_flag){
    if(this->p->block_fitting[1]>0) return;
    NDAT t;
    NPAR o, i, j;
    //    if(kg_flag == 0) { // k THIS PARAM DOESN'T MATTER HERE
    for(t=1; t<dt->n; t++) {
        //            o = dt->obs[t];
        o = this->p->dat_obs->get( dt->ix[t] );
        for(i=0; i<this->p->nS /*&& fitparam[1]>0*/; i++)
            for(j=0; j<this->p->nS; j++)
                fb->gradA[i][j] -= dt->beta[t][j] * ((o<0)?1:getB(dt,j,o)) * dt->alpha[t-1][i] / safe0num(dt->p_O_param) + L2penalty(this->p,getA(dt,i,j)); // PENALTY
    }
    //    }
}

void HMMProblem::setGradB (struct data* dt, FitBit *fb, NPAR kg_flag){
    if(this->p->block_fitting[2]>0) return;
    NDAT t;
    NPAR o, i;
    //    if(kg_flag == 0) { // k THIS PARAM DOESN'T MATTER HERE
    for(t=0; t<dt->n; t++) {
        //            o = dt->obs[t];
        o = this->p->dat_obs->get( dt->ix[t] );
        if(o<0) // if no observation -- skip
            continue;
        for(i=0; i<this->p->nS /*&& fitparam[2]>0*/; i++)
            fb->gradB[i][o] -= dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * getB(dt,i,o)) + L2penalty(this->p,getB(dt,i,o)); // PENALTY
    }
    //    }
}

void HMMProblem::computeGradients(NCAT xndat, struct data** x_data, FitBit *fb, NPAR kg_flag){//,  NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB)
    //    NPAR nS = this->p->nS;
    fb->toZero(FBS_GRAD);
    
//    clock_t tm = clock();
//    for(int i=0;i<100; i++)
    computeAlphaAndPOParam(xndat, x_data);
//    printf("Oldtime %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));

//    tm = clock();
//    for(int i=0;i<100; i++)
    computeBeta(xndat, x_data);
//    printf("Oldtime %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
	NCAT x;
    //	NDAT t;
    //	NPAR i, j, o;
    struct data *dt;
    
	for(x=0; x<xndat; x++) {
        dt = x_data[x];
		if( dt->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
        
        if(fb->PI != NULL && this->p->block_fitting[0]==0) setGradPI(dt, fb, kg_flag);
        if(fb->A  != NULL && this->p->block_fitting[1]==0) setGradA(dt, fb, kg_flag);
        if(fb->B  != NULL && this->p->block_fitting[2]==0) setGradB(dt, fb, kg_flag);
    }// for all sequences
    
    //	for(x=0; x<xndat; x++) {
    //        dt = x_data[x];
    //		if( dt->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
    //		// Mask for gradPI is handled differently
    //		// Gradient with respect to PI
    //		t = 0;
    //		o = x_data[x]->obs[t];
    //		for(i=0; i<nS /*&& fitparam[0]>0*/; i++) {
    //			fb->gradPI[i] -= x_data[x]->beta[t][i] * getB(x_data[x],i,o) / safe0num(x_data[x]->p_O_param);
    //        }
    //
    //		for(t=0; t<x_data[x]->n; t++) {
    //			o = x_data[x]->obs[t];
    //			// Gradient with respect to A
    //			// \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
    //			if( t>0 ) {
    //				for(i=0; i<nS /*&& fitparam[1]>0*/; i++)
    //					for(j=0; j<nS; j++) {
    //						fb->gradA[i][j] -= x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
    //                    }
    //
    //			}// if not first obs in sequence
    //			// Gradient with respect to B
    //			for(i=0; i<nS /*&& fitparam[2]>0*/; i++) {
    //				fb->gradB[i][o] -= x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * getB(x_data[x],i,o));
    //            }
    //		} // for all observations within skill-group
    //	} // for all groups in skill
} // computeGradients()

void HMMProblem::toFile(const char *filename) {
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            toFileSkill(filename);
            break;
        case STRUCTURE_GROUP:
            toFileGroup(filename);
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
}

void HMMProblem::toFileSkill(const char *filename) {
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
			fprintf(fid,"%10.8f%s",this->PI[k][i],(i==(this->p->nS-1))?"\n":"\t");
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

void HMMProblem::toFileGroup(const char *filename) {
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
	NCAT g;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<this->p->nG;g++) {
		it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		NPAR i,j,m;
		fprintf(fid,"PI\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PI[g][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%10.8f%s",this->A[g][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%10.8f%s",this->B[g][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblem::producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt) {
    NPAR m, i;
    NCAT k;
    NUMBER *local_pred_inner = init1D<NUMBER>(this->p->nO);
    for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
    for(int l=0; l<nks; l++) {
        for(m=0; m<this->p->nO; m++) local_pred_inner[m] = 0.0;
        k = ks[l];
        dt->k = k;
        for(m=0; m<this->p->nO; m++)
            for(i=0; i<this->p->nS; i++)
                local_pred_inner[m] += group_skill_map[dt->g][k][i] * getB(dt,i,m);//B[i][m];
        for(m=0; m<this->p->nO; m++)
            local_pred[m] += local_pred_inner[m]; // local_pred[m] = 0.0;
    }
    if(nks>1) {
        for(m=0; m<this->p->nO; m++)
            local_pred[m] /= nks;
        //            projectsimplex(local_pred, this->p->nO);
    }
    free(local_pred_inner);
}

void HMMProblem::predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill, bool only_unlabeled) {
	NDAT t;
	NCAT g, k;
	NPAR i, j, m, o, isTarget;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
	NUMBER *local_pred = init1D<NUMBER>(nO); // local prediction
//	char local_know[1024];
	NUMBER pLe[nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3D<NUMBER>(nG, nK, nS);
    NUMBER ll = 0.0, rmse = 0.0, rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NUMBER p;
    FILE *fid; // file for storing prediction should that be necessary
    bool output_this; // flag for turning on/off the writing out
    if(this->p->predictions>0) {
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr,"Can't write output model file %s\n",filename);
            exit(1);
        }
    }
	// initialize
    struct data* dt = new data;
	
    // temporary experimental: IRT-like for fitting pLo in liblinear
    /*
    FILE *fid0 = fopen("uopx_irt.txt","w");
    NPAR **group_skill_mask = init2D<NPAR>(nG, nK);
    NCAT g_k;
    data *dat;
    NPAR obs;
	for(g=0; g<nG; g++) {
        g_k = this->p->g_numk[g];
		for(k=0; k<g_k; k++) {
            dat = this->p->g_k_data[g][ k ];
            t = dat->ix[0];
            NCAT *ar;
            int n = 0;
            if(this->p->multiskill==0) {
                k = dat_skill->get(t);
                ar = &k;
                n = 1;
            } else {
                ar = &dat_multiskill->get(t)[1];
                n = dat_multiskill->get(t)[0];
                qsortNcat(ar, n);
            }
            obs = this->p->dat_obs->get( dat->ix[0] );
            NPAR count = 0; // 557687 -> 499117
            for(int l=0; l<n; l++)
                count += group_skill_mask[g][ ar[l] ] == 1;
            if(count<n) {
                fprintf(fid0,"%s %u:1", ((1-obs)==0)?"-1":"+1",dat->g+1);
                
                for(int l=0; l<n; l++) {
                    fprintf(fid0, " %u:1",ar[l]+nG+1);
                    group_skill_mask[g][ ar[l] ] = 1;
                }
                fprintf(fid0,"\n");
            }
        }
    }
    fclose(fid0);
    free2D(group_skill_mask, nG);
    // ^^
    */
    
	for(g=0; g<nG; g++)
		for(k=0; k<nK; k++) {
            dt->k = k;
            dt->g = g;
			for(i=0; i<nO; i++)
                group_skill_map[g][k][i] =  getPI(dt,i);//PI[i];
		}
	
	for(t=0; t<this->p->N; t++) {
        output_this = true;
		o = dat_obs->get(t);//[t];
        if( only_unlabeled && o>-1 ) // if we only output predictions for unlabelled, it's labelled - turn off
            output_this = false;
		g = dat_group->get(t);//[t];
        dt->g = g;
        if(!only_unlabeled) // this means we were not predicting in the first place
            isTarget = this->p->metrics_target_obs == o;
        NCAT *ar;
        int n;
        if(this->p->multiskill==0) {
            k = dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &dat_multiskill->get(t)[1];
            n = dat_multiskill->get(t)[0];
        }
        // deal with null skill
        if(ar[0]<0) { // if no skill label
            
            if(!only_unlabeled) { // this means we were not predicting in the first place
                isTarget = this->null_skill_obs==o;
                rmse += pow(isTarget - this->null_skill_obs_prob,2);
                accuracy += isTarget == (this->null_skill_obs_prob>=0.5);
                ll -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1 - this->null_skill_obs_prob);
            }
            if(this->p->predictions>0 && output_this) // write predictions file if it was opened
                for(m=0; m<nO; m++)
                    fprintf(fid,"%10.8f%s",this->null_obs_ratio[m],(m<(nO-1))?"\t":"\n");
            continue;
        }
        // produce prediction and copy to result
        producePCorrect(group_skill_map, local_pred, ar, n, dt);
        // update pL
        for(int l=0; l<n; l++) {
            //for(m=0; m<nO; m++) local_pred_inner[m] = 0.0;
            k = ar[l];
            dt->k = k;
            
            if(o>-1) { // known observations
                // update p(L)
                pLe_denom = 0.0;
                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
                for(i=0; i<nS; i++) pLe_denom += group_skill_map[g][k][i] * getB(dt,i,o);
                for(i=0; i<nS; i++) pLe[i] = group_skill_map[g][k][i] * getB(dt,i,o) / safe0num(pLe_denom);
                // 2. L = (pLe'*A)';
                for(i=0; i<nS; i++) group_skill_map[g][k][i] = 0.0;
                for(j=0; j<nS; j++)
                    for(i=0; i<nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * getA(dt,i,j);//A[i][j];
            } else { // unknown observation
                // 2. L = (pL'*A)';
                for(i=0; i<nS; i++) pLe[i] = group_skill_map[g][k][i]; // copy first;
                for(i=0; i<nS; i++) group_skill_map[g][k][i] = 0.0; // erase old value
                for(j=0; j<nS; j++)
                    for(i=0; i<nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * getA(dt,i,j);
            }// ibservations
        }
//        local_know[0] = 0;
//        for(int l=0; l<n; l++)
//            sprintf(local_know,"%s%s%10.8f",local_know,(strlen(local_know)>0)?"\t":"",group_skill_map[g][ ar[l] ][0]);
        if(this->p->predictions>0 && output_this) { // write predictions file if it was opened
            for(m=0; m<nO; m++)
                fprintf(fid,"%10.8f%s",local_pred[m],(m<(nO-1))?"\t": ((this->p->predictions==1)?"\n":"\t") );// if we print states of KCs, continut
            if(this->p->predictions==2) { // if we print out states of KC's as welll
                for(int l=0; l<n; l++) { // all KC here
//                    for(NPAR i=0; i<nS; i++) // all states
//                        fprintf(fid,"%10.8f%s",group_skill_map[g][ ar[l] ][i], (l==(n-1) && i==(nS-1))?"\n":"\t"); // if end of all states: end line
                    fprintf(fid,"%10.8f%s",group_skill_map[g][ ar[l] ][0], (l==(n-1) && l==(n-1))?"\n":"\t"); // if end of all states: end line
                }
            }
        }
        if(!only_unlabeled) { // this means we were not predicting in the first place
            rmse += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
            rmse_no_null += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
            accuracy += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5);
            accuracy_no_null += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5);
            p = safe01num(local_pred[this->p->metrics_target_obs]);
            ll -= safelog(  p)*   isTarget  +  safelog(1-p)*(1-isTarget);
            
            // temporary experimental
            for(int l=0; this->p->per_kc_rmse_acc && this->p->kc_counts!=NULL && l<n; l++) {
                //for(m=0; m<nO; m++) local_pred_inner[m] = 0.0;
                k = ar[l];
                this->p->kc_counts[k]++;
                this->p->kc_rmse[k] += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
                this->p->kc_acc[k]  += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5);
            }
        }
	} // for all data
    
    // temporary experimental
    for(int k=0; this->p->per_kc_rmse_acc && this->p->kc_counts!=NULL && k<this->p->nK; k++) {
        this->p->kc_rmse[k] = sqrt(this->p->kc_rmse[k] / this->p->kc_counts[k]);
        this->p->kc_acc[k]  =      this->p->kc_acc[k]  / this->p->kc_counts[k];
    }
    
    
    delete(dt);
	free(local_pred);
    //	free(local_pred_inner);
    free3D<NUMBER>(group_skill_map, nG, nK);
    if(!only_unlabeled) { // this means we were not predicting in the first place
        rmse = sqrt(rmse / this->p->N);
        rmse_no_null = sqrt(rmse_no_null / (this->p->N - this->p->N_null));
        if(metrics != NULL) {
            metrics[0] = ll;
            metrics[1] = 2*this->n_params + 2*ll;
            metrics[2] = this->n_params*safelog(this->p->N) + 2*ll;
            metrics[3] = rmse;
            metrics[4] = rmse_no_null;
            metrics[5] = accuracy/this->p->N;
            metrics[6] = accuracy_no_null/(this->p->N-this->p->N_null);
        }
    }
    if(this->p->predictions>0) // close predictions file if it was opened
        fclose(fid);
}

//void HMMProblem::computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data) {
//    NCAT x;
//    NPAR o,m,i,j, isTarget;
//    NPAR nS = this->p->nS, nO = this->p->nO;
//    NDAT t;
//    NUMBER *local_pred = init1D<NUMBER>(nO);
//    NUMBER *pL = init1D<NUMBER>(nS);
//    NUMBER *pLe = init1D<NUMBER>(nS);
//    NUMBER pLe_denom = 0;
//    NUMBER prob;
//    NDAT N = 0;
//    for(x=0; x<xndat; x++) {
//		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
//        N += x_data[x]->n;
//        for(i=0;i<nO; i++) pL[i] = safe01num(getPI(x_data[x],i)); // /*safe01num*/(a_PI[i]); // init pL
//        for(t=0; t<x_data[x]->n; t++) { // for all
////            o = x_data[x]->obs[t];
//            o = this->p->dat_obs->get( x_data[x]->ix[t] );
//            isTarget = (this->p->metrics_target_obs == o);
//            // produce prediction
//            for(m=0; m<nO; m++) local_pred[m] = 0.0; // init pCorr
//            for(m=0; m<nO; m++)
//                for(i=0; i<nS; i++)
//                    local_pred[m] += pL[i] * getB(x_data[x],i,m);//a_B[i][m];
//            prob = safe01num(local_pred[this->p->metrics_target_obs]);
//            loglik_rmse[0] -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
//            loglik_rmse[1] += pow(isTarget - prob, 2);
//            loglik_rmse[2] += pow(isTarget - prob, 2); // for RMSE without null skill
//            loglik_rmse[5] += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5); // acccuracy all
//            loglik_rmse[6] += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5); // acccuracy no null
//            // update p(L)
//            pLe_denom = 0;
//            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//            for(i=0; i<nS; i++)
//                pLe_denom += (pL[i]) * ((o<0)?1:getB(x_data[x],i,o)); // if observatiob unknown use 1
//            for(i=0; i<nS; i++)
//                pLe[i] = pL[i] * ((o<0)?1:getB(x_data[x],i,o)) / pLe_denom; // if observatiob unknown use 1
////            projectsimplex(pLe, param->nS);
//            // 2. L = (pLe'*A)';
//            for(i=0; i<nS; i++) pL[i] = 0.0;
//            for(j=0; j<nS; j++)
//                for(i=0; i<nS; i++)
//                    pL[j] += safe01num(pLe[i] * getA(x_data[x],i,j));//a_A[i][j]);
//        } // for all data
//    }
//    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
//    if(!keep_SE) loglik_rmse[2] = sqrt(loglik_rmse[2]/N);
//    free(local_pred);
//    free(pLe);
//    free(pL);
//}
//
NUMBER HMMProblem::getLogLik() { // get log likelihood of the fitted model
    return neg_log_lik;
}

NCAT HMMProblem::getNparams() {
    return this->n_params;
}

NUMBER HMMProblem::getNullSkillObs(NPAR m) {
    return this->null_obs_ratio[m];
}

void HMMProblem::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do RMSE*/);
    switch(this->p->solver)
    {
        case METHOD_BW: // Conjugate Gradient Descent
            loglik_rmse[0] += BaumWelchSkill();
            break;
        case METHOD_GD: // Gradient Descent
        case METHOD_CGD: // Gradient Descent
            loglik_rmse[0] += GradientDescent();
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

//void HMMProblem::computeMetrics(NUMBER* metrics) {
//    computeLogLikRMSENullSkill(metrics, true /* only SE*/);
//    // despite cycling on k-skill, would work for all
//    for(NCAT k=0; k<this->p->nK; k++)
//        computeLogLikRMSE(metrics, true /* only SE*/, this->p->k_numg[k], this->p->k_g_data[k]);
//    metrics[3] = metrics[1]; // move Squared Errors from position 2
//    metrics[3] = sqrt(metrics[3]/this->p->N);  // convert SE to RMSE
//    metrics[4] = metrics[2]; // move Squared Errors from position 2
//    metrics[4] = sqrt(metrics[4]/(this->p->N - this->p->N_null));  // convert SE to RMSE
//    metrics[5] = metrics[5] / this->p->N; // Accuracy
//    metrics[6] = metrics[6] / (this->p->N - this->p->N_null ); // Accuracy no null
//    
//    metrics[1] = 2*this->n_params + 2*metrics[0]/*loglik*/;  // AIC
//    metrics[2] = this->n_params*safelog(this->p->N) + 2*metrics[0]/*loglik*/;  // BIC
//}

void HMMProblem::FitNullSkill(NUMBER* loglik_rmse, bool keep_SE) {
    if(this->p->n_null_skill_group==0) {
        this->null_obs_ratio[0] = 1; // set first obs to 1, simplex preserved
        return; // 0 loglik
    }
    NDAT count_all_null_skill = 0;
    struct data *dat; // used as pointer
    // count occurrences
    NCAT g;
    NDAT t;
    NPAR isTarget, o, m;
    for(g=0; g<this->p->n_null_skill_group; g++) {
        dat = &this->p->null_skills[g];
        if(dat->cnt != 0)
            continue; // observe block
        count_all_null_skill += dat->n;
        for(t=0; t<dat->n; t++) {
            //            o = dat->obs[t];
            o = this->p->dat_obs->get( dat->ix[t] );
            if(((int)o)>=0) // -1 we skip \xff in char notation
                this->null_obs_ratio[ o ]++;
        }
    }
    // produce means
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    for(m=0; m<this->p->nO; m++) {
        this->null_obs_ratio[m] /= count_all_null_skill;
        if( this->null_obs_ratio[m] > this->null_skill_obs_prob ) {
            this->null_skill_obs_prob = this->null_obs_ratio[m];
            this->null_skill_obs = m;
        }
    }
    this->null_skill_obs_prob = safe01num(this->null_skill_obs_prob); // safety for logging
    // compute loglik
    NDAT N = 0;
    for(g=0; g<this->p->n_null_skill_group; g++) {
        dat = &this->p->null_skills[g];
        if(dat->cnt != 0)
            continue; // observe block
        for(t=0; t<dat->n; t++) {
            N += dat->n;
            //            isTarget = dat->obs[t] == this->null_skill_obs;
            isTarget = this->p->dat_obs->get( dat->ix[t] ) == this->null_skill_obs;
            loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
            loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
        }
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
}

//void HMMProblem::computeLogLikRMSENullSkill(NUMBER* loglik_rmse, bool keep_SE) {
//    // compute loglik
//    NDAT N = 0;
//    NPAR isTarget;
//    for(NCAT g=0; g<this->p->n_null_skill_group; g++) {
//        struct data *dat = &this->p->null_skills[g];
//        N += dat->n;
//        for(NDAT t=0; t<dat->n; t++) {
//            //            isTarget = dat->obs[t] == this->null_skill_obs;
//            isTarget = this->p->dat_obs->get( dat->ix[t] ) == this->null_skill_obs;
//            loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
//            loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
//            loglik_rmse[5] += isTarget == (this->null_skill_obs_prob>=0.5); // acccuracy all (nulls here)
//        }
//    }
//    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
//}

void HMMProblem::init3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO) {
    PI = init1D<NUMBER>((NDAT)nS);
    A  = init2D<NUMBER>((NDAT)nS, (NDAT)nS);
    B  = init2D<NUMBER>((NDAT)nS, (NDAT)nO);
}

void HMMProblem::toZero3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO) {
    toZero1D<NUMBER>(PI, (NDAT)nS);
    toZero2D<NUMBER>(A,  (NDAT)nS, (NDAT)nS);
    toZero2D<NUMBER>(B,  (NDAT)nS, (NDAT)nO);
}

void HMMProblem::cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO) {
    cpy1D<NUMBER>(soursePI, targetPI, (NDAT)nS);
    cpy2D<NUMBER>(sourseA,  targetA,  (NDAT)nS, (NDAT)nS);
    cpy2D<NUMBER>(sourseB,  targetB,  (NDAT)nS, (NDAT)nO);
}

void HMMProblem::free3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS) {
    free(PI);
    free2D<NUMBER>(A, (NDAT)nS);
    free2D<NUMBER>(B, (NDAT)nS);
    PI = NULL;
    A  = NULL;
    B  = NULL;
}

FitResult HMMProblem::GradientDescentBit(NCAT x, NCAT xndat, struct data** x_data, NPAR kg_flag, FitBit *fb, bool is1SkillForAll) {
    FitResult fr;
    fr.iter = 1;
    fr.pO0  = 0.0;
    fr.pO   = 0.0;
    fr.conv = 0; // converged
    while( !fr.conv && fr.iter<=this->p->maxiter ) {
//        clock_t tm = clock();
//        for(int i=0;i<100; i++)
        computeGradients(xndat, x_data, fb, kg_flag);//a_gradPI, a_gradA, a_gradB);
//        printf("Oldtime %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));

        if(fr.iter==1)
            fr.pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
        // copy parameter values
        fb->copy(FBS_PAR, FBS_PARm1);
        // make step
        if( fr.iter==1 || this->p->solver!=METHOD_CGD)
            doLinearStep(xndat, x_data, fb, (is1SkillForAll)?0:-1/*copy KC0*/); // step for linked skill 0
        else
            doConjugateLinearStep(xndat, x_data, fb, (is1SkillForAll)?0:-1/*copy KC0*/);
        // if we fit one skill for all
        if(is1SkillForAll) {
            NCAT nY = (kg_flag==0/*by skill*/)?this->p->nK:this->p->nG;
            for(NCAT y=0; y<nY; y++) { // copy the rest
                if(nY==x) continue;
                NUMBER *aPI = this->getPI(y);
                NUMBER **aA = this->getA(y);
                NUMBER **aB = this->getB(y);
                cpy3Params(fb->PI, fb->A, fb->B, aPI, aA, aB, this->p->nS, this->p->nO);
            }
        }
        // converge?
        fr.conv = fb->checkConvergence();
        // report if converged
        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
            computeAlphaAndPOParam(xndat, x_data);
            fr.pO = HMMProblem::getSumLogPOPara(xndat, x_data);
        } else if (this->p->solver==METHOD_CGD) {
            fb->copy(FBS_GRAD, FBS_GRADm1);
        }
        fr.iter ++;
    }// single skill loop
    RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
    fr.iter--;
    return fr;
}

NUMBER HMMProblem::GradientDescent() {
	NCAT x;
    NUMBER loglik = 0;
    NCAT nX;
    if(this->p->structure==STRUCTURE_SKILL)
        nX = this->p->nK;
    else if (this->p->structure==STRUCTURE_GROUP)
        nX = this->p->nG;
    else
        exit(1);
	FitResult fr;
    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    if(this->p->solver==METHOD_CGD) {
        fb->init(FBS_GRADm1);
        fb->init(FBS_DIRm1);
    }
	
//    // time testing
//    clock_t tm = clock();
//    for(int i=0;i<100; i++)
//        computeAlphaAndPOParam(this->p->na, this->p->k_data);
//    printf("Oldtime %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
    
    
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill>0) {
        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(0/*use skill 0*/, this->p->nSeq, this->p->k_data, 0/* by skill*/, fb, true /*is1SkillForAll*/);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
	}
	//
	// Main fit
	//
	for(x=0; x<nX && this->p->single_skill!=2; x++) { // if not "force single skill" too
        NCAT xndat;
        struct data** x_data;
        if(this->p->structure==STRUCTURE_SKILL) {
            xndat = this->p->k_numg[x];
            x_data = this->p->k_g_data[x];
        } else if(this->p->structure==STRUCTURE_GROUP) {
            xndat = this->p->g_numk[x];
            x_data = this->p->g_k_data[x];
        } else {
            xndat = 0;
            x_data = NULL;
        }
        fb->linkPar( this->getPI(x), this->getA(x), this->getB(x));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(x/*use skill x*/, xndat, x_data, (this->p->structure!=STRUCTURE_SKILL), fb, false /*is1SkillForAll*/);
        if( !this->p->quiet )
            printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", x, fr.iter,fr.pO0,fr.pO,fr.conv);
	} // for all skills
    delete fb; // that takes care of *m1, and *GRAD
    return loglik;
}

NUMBER HMMProblem::BaumWelchSkill() {
	NCAT k;
    NUMBER loglik = 0;
	
    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
    fb->init(FBS_PARm1);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    //    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
    //	if(this->p->fit_single_skill) {
    //		for(k=0; k<this->p->nK; k++) this->p->mask_skill[k] = true; // mask with k'th==true to be computed
    //		iter = 1;
    //		pO0 = 0.0;
    //		conv = 0; // converged
    //		while( !conv && iter<=this->p->maxiter ) {
    //			if(iter>1) {
    //				toZero1DNumber(gradPI, this->p->nS);
    //				toZero2DNumber(gradA,  this->p->nS, this->p->nS);
    //				toZero2DNumber(gradB,  this->p->nS, this->p->nO);
    //			}
    //			hmm->computeGradients();
    //			if(iter==1)
    //				pO0 = hmm->getSumLogPOParaK(flog);
    //			// add gradients
    //			for(k=0; k<this->p->nK; k++) {
    //				add1DNumbersWeighted(hmm->getGradPI(k), gradPI, this->p->nS, 1.0);
    //				add2DNumbersWeighted(hmm->getGradA(k),  gradA,  this->p->nS, this->p->nS, 1.0);
    //				add2DNumbersWeighted(hmm->getGradB(k),  gradB,  this->p->nS, this->p->nO, 1.0);
    //			}
    //			// copy old SAVED! values for params, just for skill #0 is enough
    //			cpy1DNumber(hmm->getPI(0), PI_m1, this->p->nS);
    //			cpy2DNumber(hmm->getA(0),  A_m1,  this->p->nS, this->p->nS);
    //			cpy2DNumber(hmm->getB(0),  B_m1,  this->p->nS, this->p->nO);
    //
    //			// make step
    //			for(k=0; k<this->p->nK; k++)
    //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
    //			// check convergence, on any skill, e.g. #0
    //			conv = checkConvergence(hmm->getPI(k), hmm->getA(k), hmm->getB(k), PI_m1, A_m1, B_m1, param);
    //
    //			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||/**/ (conv || iter==this->p->maxiter) )) {
    //				hmm->computeAlphaAndPOParam();
    //				NUMBER pO = hmm->getSumLogPOParaK(flog);
    //				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
    //			}
    //			iter ++;
    //		}// single skill loop
    //
    //		free(gradPI);				   // from now on, just a pointer
    //		free2DNumber(gradA, this->p->nS); // from now on, just a pointer
    //		free2DNumber(gradB, this->p->nS); // from now on, just a pointer
    //	}
	
	//
	// Main fit
	//
    //	for(k=0; k<this->p->nK; k++) this->p->mask_skill[k] = false; // mask with k'th==true to be computed
	for(k=0; k<this->p->nK; k++) {  // for(k=218; k<219; k++) { //
        NCAT xndat = this->p->k_numg[k];
        struct data** x_data =  this->p->k_g_data[k];
        
		conv = 0; // converged
		iter = 1; // iteration count
		pO0 = 0.0;
        pO = 0.0;
		
        fb->linkPar(this->getPI(k), this->getA(k), this->getB(k));
		
		while( !conv && iter<=this->p->maxiter ) {
			if(iter==1) {
                computeAlphaAndPOParam(xndat, x_data);
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
            fb->copy(FBS_PAR, FBS_PARm1);
			
            //			hmm->zeroLabelsK(k); // reset blocking labels // THIS IS NOT DONE HERE
			doBaumWelchStep(xndat, this->p->k_g_data[k], fb);// PI, A, B);
            
			// check convergence
            conv = fb->checkConvergence();
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                computeAlphaAndPOParam(xndat, x_data);
                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                loglik += pO*(pO>0);
                if(!this->p->quiet)
                    printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",k,iter,pO0,pO,conv);
			}
			
			iter ++;
		} // main solver loop
		// recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
        RecycleFitData(xndat, x_data, this->p);
		// recycle
	} // for all skills
    delete fb;
    return loglik;
}

NUMBER HMMProblem::doLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy) {//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO;
	// first scale down gradients
    fb->doLog10ScaleGentle(FBS_GRAD);
    //	doLog10Scale1DGentle(fb->gradPI, fb->PI, nS);
    //	doLog10Scale2DGentle(fb->gradA,  fb->A,  nS, nS);
    //	doLog10Scale2DGentle(fb->gradB,  fb->B,  nS, nO);
	
    //    NUMBER *PI_cpy, ** A_cpy, ** B_cpy; // replace with PARcopy in fb
    //    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
    fb->init(FBS_PARcopy);
    
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
    //	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
    fb->copy(FBS_PAR, FBS_PARcopy);
    //    cpy3Params(fb->PI, fb->A, fb->B, PI_cpy, A_cpy, B_cpy, nS, nO);
    //	cpy1DNumber(fb->PI, PI_cpy, nS); // save copy
    //	cpy2DNumber(fb->A,  A_cpy,  nS, nS); // save copy
    //	cpy2DNumber(fb->B,  B_cpy,  nS, nO); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++)
	{
		if(fb->PI != NULL) p_k_by_neg_p_k -= fb->gradPI[i]*fb->gradPI[i];
		if(fb->A  != NULL) for(j=0; j<nS; j++) p_k_by_neg_p_k -= fb->gradA[i][j]*fb->gradA[i][j];
		if(fb->B  != NULL) for(m=0; m<nO; m++) p_k_by_neg_p_k -= fb->gradB[i][m]*fb->gradB[i][m];
	}
	int iter = 0; // limit iter steps to 20, via ArmijoMinStep (now 10)
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			if(fb->PI != NULL) fb->PI[i] = fb->PIcopy[i] - e * fb->gradPI[i];
            if(fb->A  != NULL)
                for(j=0; j<nS; j++)
                    fb->A[i][j] = fb->Acopy[i][j] - e * fb->gradA[i][j];
            if(fb->B  != NULL)
                for(m=0; m<nO; m++)
                    fb->B[i][m] = fb->Bcopy[i][m] - e * fb->gradB[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			if(fb->PI != NULL) projectsimplex(fb->PI, nS);
			for(i=0; i<nS; i++) {
				if(fb->A  != NULL) projectsimplex(fb->A[i], nS);
				if(fb->B  != NULL) projectsimplex(fb->B[i], nS);
			}
		} else {
			if(fb->PI != NULL) projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
			for(i=0; i<nS; i++) {
				if(fb->A  != NULL) projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
				if(fb->B  != NULL) projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
			}
		}
        if(copy >= 0) { // copy parameters from position 'copy' to others
            for(NCAT x=0; (fb->PI != NULL) && x<sizes[0]; x++)
                if(x!=copy)
                    cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
            for(NCAT x=0; (fb->A  != NULL) && x<sizes[1]; x++)
                if(x!=copy)
                    cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
            for(NCAT x=0; (fb->B  != NULL) && x<sizes[2]; x++)
                if(x!=copy)
                    cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
        }
		// recompute alpha and p(O|param)
		computeAlphaAndPOParam(xndat, x_data);
		// compute f(x_{k+1})
		f_xkplus1 = HMMProblem::getSumLogPOPara(xndat, x_data);
		// compute Armijo compliance
		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
		iter++;
	} // armijo loop
    if(!compliesArmijo) { // we couldn't step away from current, copy the inital point back
        e = 0;
        fb->copy(FBS_PARcopy, FBS_PAR);
    }
    //    RecycleFitData(xndat, x_data, this->p);
    fb->destroy(FBS_PARcopy);
    return e;
} // doLinearStep

NUMBER HMMProblem::doConjugateLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO;
	// first scale down gradients
    fb->doLog10ScaleGentle(FBS_GRAD);
    
    // compute beta_gradient_direction
    NUMBER beta_grad = 0, beta_grad_num = 0, beta_grad_den = 0;
    
    switch (this->p->solver_setting) {
        case 1: // Fletcher-Reeves
            for(i=0; i<nS; i++)
            {
                if(fb->PI != NULL) {
                    beta_grad_num = fb->gradPI  [i]*fb->gradPI   [i];
                    beta_grad_den = fb->gradPIm1[i]*fb->gradPIm1[i];
                }
                if(fb->A  != NULL)
                    for(j=0; j<nS; j++) {
                        beta_grad_num = fb->gradA  [i][j]*fb->gradA   [i][j];
                        beta_grad_den = fb->gradAm1[i][j]*fb->gradAm1[i][j];
                    }
                if(fb->B  != NULL)
                    for(m=0; m<nO; m++) {
                        beta_grad_num = fb->gradB  [i][m]*fb->gradB  [i][m];
                        beta_grad_den = fb->gradBm1[i][m]*fb->gradBm1[i][m];
                    }
            }
            break;
        case 2: // Polak–Ribiere
            for(i=0; i<nS; i++)
            {
                if(fb->PI != NULL) {
                    beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                    beta_grad_den =  fb->gradPIm1[i]*fb->gradPIm1[i];
                }
                if(fb->A != NULL)
                    for(j=0; j<nS; j++) {
                        beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                        beta_grad_den =  fb->gradAm1[i][j]*fb->gradAm1[i][j];
                    }
                if(fb->B  != NULL)
                    for(m=0; m<nO; m++) {
                        beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                        beta_grad_den =  fb->gradBm1[i][m]*fb->gradBm1[i][m];
                    }
            }
            break;
        case 3: // Hestenes-Stiefel
            for(i=0; i<nS; i++)
            {
                if(fb->PI != NULL) {
                    beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                    beta_grad_den =  fb->dirPIm1[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                }
                if(fb->A  != NULL)
                    for(j=0; j<nS; j++) {
                        beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                        beta_grad_den =  fb->dirAm1[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                    }
                if(fb->B  != NULL)
                    for(m=0; m<nO; m++) {
                        beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                        beta_grad_den =  fb->dirBm1[i][m]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                    }
            }
            break;
        default:
            fprintf(stderr,"Wrong stepping algorithm (%d) specified for Conjugate Gradient Descent\n",this->p->solver_setting);
            break;
    }
    beta_grad = beta_grad_num / safe0num(beta_grad_den);
    beta_grad = (beta_grad>=0)?beta_grad:0;
	
    // compute new direction (in place of old)
    fb->toZero(FBS_DIRm1);
	for(i=0; i<nS; i++)
	{
		if(fb->PI != NULL) fb->dirPIm1[i] = -fb->gradPI[i] + beta_grad * fb->dirPIm1[i];
		if(fb->A  != NULL) for(j=0; j<nS; j++) fb->dirAm1[i][j] = -fb->gradA[i][j] + beta_grad * fb->dirAm1[i][j];
		if(fb->B  != NULL) for(m=0; m<nO; m++) fb->dirBm1[i][m] = -fb->gradB[i][m] + beta_grad * fb->dirBm1[i][m];
	}
	// scale down direction
    fb->doLog10ScaleGentle(FBS_DIRm1);
    
    fb->init(FBS_PARcopy);
    
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
    fb->copy(FBS_PAR, FBS_PARcopy);
	// compute p_k * -p_k >>>> now current gradient by current direction
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++)
	{
		if(fb->PI != NULL) p_k_by_neg_p_k = fb->gradPI[i]*fb->dirPIm1[i];
		if(fb->A  != NULL) for(j=0; j<nS; j++) p_k_by_neg_p_k = fb->gradA[i][j]*fb->dirAm1[i][j];
		if(fb->B  != NULL) for(m=0; m<nO; m++) p_k_by_neg_p_k = fb->gradB[i][m]*fb->dirBm1[i][m];
	}
	int iter = 0; // limit iter steps to 20, actually now 10, via ArmijoMinStep
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			if(fb->PI != NULL) fb->PI[i] = fb->PIcopy[i] + e * fb->dirPIm1[i];
            if(fb->A  != NULL)
                for(j=0; j<nS; j++)
                    fb->A[i][j] = fb->Acopy[i][j] + e * fb->dirAm1[i][j];
            if(fb->B  != NULL)
                for(m=0; m<nO; m++)
                    fb->B[i][m] = fb->Bcopy[i][m] + e * fb->dirBm1[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			if(fb->PI != NULL) projectsimplex(fb->PI, nS);
			for(i=0; i<nS; i++) {
				if(fb->A != NULL) projectsimplex(fb->A[i], nS);
				if(fb->B != NULL) projectsimplex(fb->B[i], nS);
			}
		} else {
			if(fb->PI != NULL) projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
			for(i=0; i<nS; i++) {
				if(fb->A != NULL) projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
				if(fb->B != NULL) projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
			}
		}
        // copy parameters from position 'copy' to others
        if(copy >= 0) {
            if(fb->PI != NULL)
                for(NCAT x=0; x<sizes[0]; x++)
                    if(x!=copy)
                        cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
            if(fb->A  != NULL)
                for(NCAT x=0; x<sizes[1]; x++)
                    if(x!=copy)
                        cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
            if(fb->B  != NULL)
                for(NCAT x=0; x<sizes[2]; x++)
                    if(x!=copy)
                        cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
        }
		// recompute alpha and p(O|param)
		computeAlphaAndPOParam(xndat, x_data);
		// compute f(x_{k+1})
		f_xkplus1 = HMMProblem::getSumLogPOPara(xndat, x_data);
		// compute Armijo compliance
		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
		iter++;
	} // armijo loop
    if(!compliesArmijo) { // failed to step away from initial, reinstate the inital parameters
        e = 0;
        fb->copy(FBS_PARcopy, FBS_PAR);
    }
    //    RecycleFitData(xndat, x_data, this->p);
    fb->destroy(FBS_PARcopy);
    return e;
} // doLinearStep

NUMBER HMMProblem::doBarzalaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
	doLog10Scale2DGentle(a_gradB,  a_B,  nS, nO);
    
    // compute s_k_m1
  	NUMBER *s_k_m1_PI = init1D<NUMBER>((NDAT)nS);
	NUMBER **s_k_m1_A = init2D<NUMBER>((NDAT)nS, (NDAT)nS);
	NUMBER **s_k_m1_B = init2D<NUMBER>((NDAT)nS, (NDAT)nS);
	for(i=0; i<nS; i++)
	{
		s_k_m1_PI[i] = a_PI[i] - a_PI_m1[i];
		for(j=0; j<nS; j++) s_k_m1_A[i][j] = a_A[i][j] - a_A_m1[i][j];
		for(m=0; m<this->   p->nO; m++) s_k_m1_B[i][m] = a_B[i][m] - a_B_m1[i][m];
	}
    // compute alpha_step
    NUMBER alpha_step = 0, alpha_step_num = 0, alpha_step_den = 0;
    // Barzalai Borweig: s' * s / ( s' * (g-g_m1) )
	for(i=0; i<nS; i++)
	{
		alpha_step_num = s_k_m1_PI[i]*s_k_m1_PI[i];
		alpha_step_den = s_k_m1_PI[i]*(a_gradPI[i] - a_gradPI_m1[i]);
		for(j=0; j<nS; j++) {
            alpha_step_num = s_k_m1_A[i][j]*s_k_m1_A[i][j];
            alpha_step_den = s_k_m1_A[i][j]*(a_gradA[i][j] - a_gradA_m1[i][j]);
        }
		for(m=0; m<nO; m++) {
            alpha_step_num = s_k_m1_B[i][m]*s_k_m1_B[i][m];
            alpha_step_den = s_k_m1_B[i][m]*(a_gradB[i][m] - a_gradB_m1[i][m]);
        }
	}
    alpha_step = alpha_step_num / safe0num(alpha_step_den);
    
    // step
    for(i=0; i<nS; i++) {
        a_PI[i] = a_PI[i] - alpha_step * a_gradPI[i];
        for(j=0; j<nS; j++)
            a_A[i][j] = a_A[i][j] - alpha_step * a_gradA[i][j];
        for(m=0; m<nO; m++)
            a_B[i][m] = a_B[i][m] - alpha_step * a_gradB[i][m];
    }
    // scale
    if( !this->hasNon01Constraints() ) {
        projectsimplex(a_PI, nS);
        for(i=0; i<nS; i++) {
            projectsimplex(a_A[i], nS);
            projectsimplex(a_B[i], nS);
        }
    } else {
        projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), nS);
        for(i=0; i<nS; i++) {
            projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], nS);
            projectsimplexbounded(a_B[i], this->getLbB()[i], this->getUbB()[i], nS);
        }
    }
	free(s_k_m1_PI);
	free2D<NUMBER>(s_k_m1_B, nS);
	free2D<NUMBER>(s_k_m1_A, nS);
    return alpha_step;
}
void HMMProblem::doBaumWelchStep(NCAT xndat, struct data** x_data, FitBit *fb) {//, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B) {
	NCAT x;
	NPAR i,j,m, o;
	NDAT t;
    //	if(!this->p->mask_skill[k]) return; // if mask is set ...
    
    //    HMMProblem::computeAlphaAndPOParam(xndat, x_data, a_PI, a_A, a_B, this->p->nS);
    //	HMMProblem::computeBeta(xndat, x_data, a_A, a_B, this->p->nS);
    computeAlphaAndPOParam(xndat, x_data);
	computeBeta(xndat, x_data);
	HMMProblem::computeXi(xndat, x_data);
	HMMProblem::computeGamma(xndat, x_data);
	
	NUMBER * b_PI = init1D<NUMBER>((NDAT)this->p->nS);
	NUMBER ** b_A_num = init2D<NUMBER>((NDAT)this->p->nS, (NDAT)this->p->nS);
	NUMBER ** b_A_den = init2D<NUMBER>((NDAT)this->p->nS, (NDAT)this->p->nS);
	NUMBER ** b_B_num = init2D<NUMBER>((NDAT)this->p->nS, (NDAT)this->p->nO);
	NUMBER ** b_B_den = init2D<NUMBER>((NDAT)this->p->nS, (NDAT)this->p->nO);
	// compute sums PI
	NUMBER sum_p_O_param = 0;
	for(x=0; x<xndat; x++) {
        if( x_data[x]->cnt!=0 ) continue;
		sum_p_O_param += x_data[x]->p_O_param;
    }
	for(x=0; x<xndat; x++) {
        if( x_data[x]->cnt!=0 ) continue;
		for(i=0; i<this->p->nS; i++)
			fb->PI[i] += x_data[x]->gamma[0][i] / xndat;
		
		for(t=0;t<(x_data[x]->n-1);t++) {
            //			o = x_data[x]->obs[t];
            o = this->p->dat_obs->get( x_data[x]->ix[t] );
			for(i=0; i<this->p->nS; i++) {
				for(j=0; j<this->p->nS; j++){
					b_A_num[i][j] += x_data[x]->xi[t][i][j];
					b_A_den[i][j] += x_data[x]->gamma[t][i];
				}
				for(m=0; m<this->p->nO; m++) {
					b_B_num[i][m] += (m==o) * x_data[x]->gamma[t][i];
					b_B_den[i][m] += x_data[x]->gamma[t][i];
				}
			}
		}
	} // for all groups within a skill
	// set params
	for(i=0; i<this->p->nS; i++) {
		fb->PI[i] = b_PI[i];
		for(j=0; j<this->p->nS; j++)
			fb->A[i][j] = b_A_num[i][j] / safe0num(b_A_den[i][j]);
		for(m=0; m<this->p->nO; m++)
			fb->B[i][m] = b_B_num[i][m] / safe0num(b_B_den[i][m]);
	}
    // scale
    if( !this->hasNon01Constraints() ) {
        projectsimplex(fb->PI, this->p->nS);
        for(i=0; i<this->p->nS; i++) {
            projectsimplex(fb->A[i], this->p->nS);
            projectsimplex(fb->B[i], this->p->nS);
        }
    } else {
        projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), this->p->nS);
        for(i=0; i<this->p->nS; i++) {
            projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], this->p->nS);
            projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], this->p->nS);
        }
    }
    // free mem
    //    RecycleFitData(xndat, x_data, this->p);
	free(b_PI);
	free2D<NUMBER>(b_A_num, this->p->nS);
	free2D<NUMBER>(b_A_den, this->p->nS);
	free2D<NUMBER>(b_B_num, this->p->nS);
	free2D<NUMBER>(b_B_den, this->p->nS);
}

//void HMMProblem::createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space) {
//    int target_state = 0; // known
//
//	long long elements, q;
//
//	prob.l = this->p->N - this->p->N_null;
//	elements = (long long)prob.l * ((long long)this->p->nK*this->p->nK + 1); // +1 more is for ?, but it's there
//
//	prob.bias= 1; // no bias
////	prob.y = Malloc(double,prob.l);
////	prob.x = Malloc(struct feature_node *,prob.l);
////	x_space = Malloc(struct feature_node,elements+prob.l); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
//    //	max_index = 0; // max number of columns
//    //
//    // vvvvv init HMM
//    NDAT /*max_index, inst_max_index,*/ t, tidx;
//	NCAT g, k;
//	NPAR i, j, m, o, isTarget;
//	NUMBER *local_pred = init1DNumber(this->p->nO); // local prediction
//	NUMBER pLe[this->p->nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3DNumber(this->p->nG, this->p->nK, this->p->nS);
//    struct data dt;
//	for(g=0; g<this->p->nG; g++) {
//        dt.g = g;
//		for(k=0; k<this->p->nK; k++) {
//            dt.k = k;
//			for(i=0; i<this->p->nO; i++)
//                group_skill_map[g][k][i] =  getPI(&dt,i);//PI[i];
//		}
//    }
//    // ^^^^^ init HMM
//    //
//    t = 0;
//	q = 0;
//    for(tidx=0; tidx<this->p->N; tidx++) {
//        // vvvvvvv
//        k = this->p->dat_skill->get(tidx);
//        if(k<0) // null skill
//            continue;
//        g = this->p->dat_group->get(tidx);
//        dt.k = k;
//        dt.g = g;
//        o = this->p->dat_obs->get(tidx); //dat_obs->get(t);//[t];
//        // ^^^^^^^
////        fprintf(stdout,"..%d:ll,",tidx);
//        //		inst_max_index = 0; // strtol gives 0 if wrong format
//        //		readline(fp); // not needed
//        prob.x[t] = &x_space[q];
//        //		label = strtok(line," \t\n"); // not needed
//        //		if(label == NULL) // empty line  // not needed
//        //			exit_input_error(i+1);       // not needed
//
//        prob.y[t] = (double)o;//dat_obs->get(t); // strtod(label,&  );
//        //		if(endptr == label || *endptr != '\0')
//        //			exit_input_error(i+1);
//        for(NCAT r=0; r<(this->p->nK); r++) //while(1)
//        {
//            // k*this->p->nK - shift, r position, 1 - idexes start with 1
//            NUMBER value = ( group_skill_map[g][r][target_state] );
//            if( value!= 0 ) {
//                x_space[q].index = k*(this->p->nK+1) + r+1; // idexes are 1-starting //(int) strtol(idx,&endptr,10);
//                x_space[q].value = logit(value); // this student, all r skills, just the target state // strtod(val,&endptr);
//            }
//            ++q;
//        }
//        ++t;
//        //		if(inst_max_index > max_index) // we know the count
//        //			max_index = inst_max_index;
//        if(prob.bias >= 0) // leave it, let bias variable control it
//            x_space[q++].value = prob.bias; // copy bias and step over
//        x_space[q++].index = -1; // set index to -1???
//        //
//        //
//        // vvvvvv make an update to the HMM
////        fprintf(stdout,"bkt");
//        // deal with null skill
//        if(k<0) continue;// if no skill label
//        isTarget = this->p->metrics_target_obs == o;
//        for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
//        // produce prediction and copy to result
//        for(m=0; m<this->p->nO; m++)
//            for(i=0; i<this->p->nS; i++)
//                local_pred[m] += group_skill_map[g][k][i] * getB(&dt,i,m);//B[i][m];
//        // update p(L)
//        pLe_denom = 0.0;
//        // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//        for(i=0; i<this->p->nS; i++) pLe_denom += group_skill_map[g][k][i] * getB(&dt,i,o);//B[i][o];
//        for(i=0; i<this->p->nS; i++) pLe[i] = group_skill_map[g][k][i] * getB(&dt,i,o)/*B[i][o]*/ / safe0num(pLe_denom);
//        // 2. L = (pLe'*A)';
//        for(i=0; i<this->p->nS; i++) group_skill_map[g][k][i] = 0.0;
//        for(j=0; j<this->p->nS; j++)
//            for(i=0; i<this->p->nS; i++)
//                group_skill_map[g][k][j] += pLe[i] * getA(&dt,i,j);//A[i][j];
//        // ^^^^^^ make an update to the HMM
//        //
//        //
////        fprintf(stdout,";\n");
//    } // for all t in {G,K}
//    //
//    // vvvvv Recycle prediction stuff
//	free(local_pred);
//    free3DNumber(group_skill_map, this->p->nG, this->p->nK);
//    // ^^^^^ Recycle prediction stuff
//    //
//	if(prob.bias >= 0) // taken care of preemptively
//	{
//        prob.n = (long long)this->p->nK*this->p->nK  + 1;
//		for(t=1;t<prob.l;t++)
//			(prob.x[t]-2)->index = prob.n;
//		x_space[q-2].index = prob.n;
//	}
//	else
//        prob.n = (long long)this->p->nK*this->p->nK ;
//    //	fclose(fp);
//    //
//    // now set up parameter
//    //
//    this->p->solver_type = L1R_LR;
//    this->p->C = 1;
//    this->p->eps = 0.01;
//    this->p->weight_label = NULL;
//    this->p->weight = NULL;
//
//}
//
//void HMMProblem::createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space, NCAT k) {
//    if(k<0 || k>=this->p->nK) {
//        fprintf(stderr, "KC specified is out of range!\n");
//        exit(1);
//    }
//
//    int target_obs = 0; // known
//
//    NDAT t, tidx;
//	long long q;
//    fprintf(stdout,"prob.l=%d\n",prob.l);
//	prob.bias= 1; // has bias
//    // vvvvv init HMM
//	NCAT g, kk;
//	NPAR i, j, m, o, isTarget;
//	NUMBER *local_pred = init1DNumber(this->p->nO); // local prediction
//	NUMBER pLe[this->p->nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3DNumber(this->p->nG, this->p->nK, this->p->nS);
//    struct data dt;
//	for(g=0; g<this->p->nG; g++) {
//        dt.g = g;
//		for(kk=0; kk<this->p->nK; kk++) {
//            dt.k = kk;
//			for(i=0; i<this->p->nO; i++)
//                group_skill_map[g][kk][i] =  getPI(&dt,i);//PI[i];
//		}
//    }
//    // ^^^^^ init HMM
//    //
//	q = 0; // position in linearized array of features
//    t = 0; // row of feature matrix
//    for(tidx=0; tidx<this->p->N; tidx++) {
//        // vvvvvvv
//        kk = this->p->dat_skill->get(tidx);
//        if(kk<0) continue;
//        g = this->p->dat_group->get(tidx);
//        o = this->p->dat_obs->get(tidx);
//        dt.k = kk;
//        dt.g = g;
//        // ^^^^^^^
//        if(kk==k) { // add this to the data
//            prob.x[t] = &x_space[q];
//            prob.y[t] = (double)o;//dat_obs->get(t); // strtod(label,&  );
//            for(NCAT r=0; r<(this->p->nK); r++) //while(1)
//            {
//                NUMBER pCorr = 0;
//                dt.k = r; // set to regression KC(s)
//                for(i=0; i<this->p->nS; i++)
//                    pCorr += group_skill_map[g][r][i] * getB(&dt,i,target_obs);//B[i][m];
//                if( pCorr!= 0 ) {
//                    x_space[q].index = r+1; // idexes are 1-starting //(int) strtol(idx,&endptr,10);
//                    x_space[q].value = logit(pCorr); // this student, all r skills, just the target state // strtod(val,&endptr);
//                }
//                ++q;
//            }
//            ++t;
//            dt.k = kk; //back to main KC
//            if(prob.bias >= 0) // leave it, let bias variable control it
//                x_space[q++].value = prob.bias; // copy bias and step over
//            x_space[q++].index = -1; // set index to -1???
//        }
//        //
//        //
//        // vvvvvv make an update to the HMM
//        // deal with null skill
//        isTarget = this->p->metrics_target_obs == o;
//        for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
//        // produce prediction and copy to result
//        for(m=0; m<this->p->nO; m++)
//            for(i=0; i<this->p->nS; i++)
//                local_pred[m] += group_skill_map[g][kk][i] * getB(&dt,i,m);//B[i][m];
//        // update p(L)
//        pLe_denom = 0.0;
//        // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//        for(i=0; i<this->p->nS; i++) pLe_denom += group_skill_map[g][kk][i] * getB(&dt,i,o);//B[i][o];
//        for(i=0; i<this->p->nS; i++) pLe[i] = group_skill_map[g][kk][i] * getB(&dt,i,o) / safe0num(pLe_denom);
//        // 2. L = (pLe'*A)';
//        for(i=0; i<this->p->nS; i++) group_skill_map[g][k][i] = 0.0;
//        for(j=0; j<this->p->nS; j++)
//            for(i=0; i<this->p->nS; i++)
//                group_skill_map[g][kk][j] += pLe[i] * getA(&dt,i,j);//A[i][j];
//        // ^^^^^^ make an update to the HMM
//        //
//        //
//    } // for all t in {G,K}
//
//    if(prob.bias >= 0) // taken care of preemptively
//	{
//        prob.n = this->p->nK + 1;
//		for(t=1;t<prob.l;t++)
//			(prob.x[t]-2)->index = prob.n;
//		x_space[q-2].index = prob.n;
//	}
//	else
//        prob.n = this->p->nK;
//    //
//    // vvvvv Recycle prediction stuff
//	free(local_pred);
//    free3DNumber(group_skill_map, this->p->nG, this->p->nK);
//    // ^^^^^ Recycle prediction stuff
//    //
//    //	fclose(fp);
//    //
//    // now set up parameter
//    //
//    this->p->solver_type = L2R_LR;
//    this->p->C = 1;
//    this->p->eps = 0.01;
//    this->p->weight_label = NULL;
//    this->p->weight = NULL;
//    this->p->nr_weight = 0;
//
//}
//
//void HMMProblem::recycleLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space) {
//	destroy_param(&param);
//	free(prob.y);
//	free(prob.x);
//	free(x_space);
//}

void HMMProblem::readNullObsRatio(FILE *fid, struct param* param, NDAT *line_no) {
	NPAR i;
	//
	// read null skill ratios
	//
    fscanf(fid, "Null skill ratios\t");
    this->null_obs_ratio =Calloc(NUMBER, (size_t)this->p->nO);
    this->null_skill_obs      = 0;
    this->null_skill_obs_prob = 0;
	for(i=0; i<param->nO; i++) {
        if( i==(param->nO-1) ) // end
            fscanf(fid,"%lf\n",&this->null_obs_ratio[i] );
        else
            fscanf(fid,"%lf\t",&this->null_obs_ratio[i] );
		
        if( this->null_obs_ratio[i] > this->null_skill_obs_prob ) {
            this->null_skill_obs_prob = this->null_obs_ratio[i];
            this->null_skill_obs = i;
        }
	}
    (*line_no)++;
}

void HMMProblem::readModel(const char *filename, bool overwrite) {
	FILE *fid = fopen(filename,"r");
	if(fid == NULL)
	{
		fprintf(stderr,"Can't read model file %s\n",filename);
		exit(1);
	}
	int max_line_length = 1024;
	char *line = Malloc(char,(size_t)max_line_length);
	NDAT line_no = 0;
    struct param initparam;
    set_param_defaults(&initparam);
    
    //
    // read solver info
    //
    if(overwrite)
        readSolverInfo(fid, this->p, &line_no);
    else
        readSolverInfo(fid, &initparam, &line_no);
    //
    // read model
    //
    readModelBody(fid, &initparam, &line_no, overwrite);
		
	fclose(fid);
	free(line);
}

void HMMProblem::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
	NPAR i,j,m;
	NCAT k = 0, idxk = 0;
    std::map<std::string,NCAT>::iterator it;
	string s;
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
	// read skills
	//
	for(k=0; k<param->nK; k++) {
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
        fscanf(fid,"PI\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->PI[idxk][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->PI[idxk][i] = atof(col);
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
