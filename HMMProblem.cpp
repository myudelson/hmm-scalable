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

#include "liblinear/linear.h"

HMMProblem::HMMProblem() {
}

HMMProblem::HMMProblem(struct param *param) {
    NPAR i;
    switch (param->solver) {
        case BKT_CGD: // Conjugate Gradient Descent
        case BKT_GD: // Gradient Descent
        case BKT_BW: // Expectation Maximization (Baum-Welch)
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            for(i=0; i<3; i++) this->sizes[i] = param->nK;
            this->n_params = param->nK * 4;
            break;
        case BKT_GD_G: // Gradient Descent by group
            for(i=0; i<3; i++) this->sizes[i] = param->nG;
            this->n_params = param->nG * 4;
            break;
        default:
            fprintf(stderr,"Method specified is not supported and should have been caught earlier\n");
            break;
    }
    init(param);
}

void HMMProblem::init(struct param *param) {
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    if( this->p->solver == BKT_CGD && this->p->solver_settting == -1)
        this->p->solver_settting = 1; // default Fletcher-Reeves
        
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
		this->PI = init2DNumber(this->sizes[0], nS);
		this->A =  init3DNumber(this->sizes[1], nS, nS);
		this->B =  init3DNumber(this->sizes[2], nS, nO);
        NCAT x;
		for(x=0; x<this->sizes[0]; x++)
			cpy1DNumber(a_PI, this->PI[x], nS);
		for(x=0; x<this->sizes[1]; x++)
			cpy2DNumber(a_A, this->A[x], nS, nS);
		for(x=0; x<this->sizes[2]; x++)
			cpy2DNumber(a_B, this->B[x], nS, nO);
	} else {
		fprintf(stderr,"params do not meet constraints.\n");
		exit(1);
	}
    // destroy setup params
	free(a_PI);
	free2DNumber(a_A, nS);
	free2DNumber(a_B, nS);
	
    
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
    //	this->gradPI = NULL;
    //	this->gradA = NULL;
    //	this->gradB = NULL;
}

HMMProblem::~HMMProblem() {
    destroy();
}

void HMMProblem::destroy() {
	// destroy model data
    if(this->null_obs_ratio != NULL) free(this->null_obs_ratio);
	if(this->PI != NULL) free2DNumber(this->PI, this->sizes[0]);
	if(this->A  != NULL) free3DNumber(this->A,  this->sizes[1], this->p->nS);
	if(this->B  != NULL) free3DNumber(this->B,  this->sizes[2], this->p->nS);
	if(this->lbPI!=NULL) free(this->lbPI);
	if(this->ubPI!=NULL) free(this->ubPI);
	if(this->lbA!=NULL) free2DNumber(this->lbA, this->p->nS);
	if(this->ubA!=NULL) free2DNumber(this->ubA, this->p->nS);
	if(this->lbB!=NULL) free2DNumber(this->lbB, this->p->nS);
	if(this->ubB!=NULL) free2DNumber(this->ubB, this->p->nS);
	// destroy fitting data
//    NPAR i;
//	// gradPI
//	if ( this->gradPI != NULL) { // allocate
//		free(this->gradPI);
//		this->gradPI = NULL;
//	}
//	// gradA
//	if ( this->gradA != NULL) { // allocate
//		for(i=0; i<this->p->nS; i++) {
//			free(this->gradA[i]);
//		}
//		free(this->gradA);
//		this->gradA = NULL;
//	}
//	// gradB
//	if ( this->gradB != NULL) { // allocate
//		for(i=0; i<this->p->nS; i++) {
//			free(this->gradB[i]);
//		}
//		free(this->gradB);
//		this->gradB = NULL;
//	}
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
    switch(this->p->solver)
    {
        case BKT_CGD: // Conjugate Gradient Descent
        case BKT_GD: // Gradient Descent
        case BKT_BW: // Expectation Maximization (Baum-Welch)
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            return this->PI[dt->k][i];
            break;
        case BKT_GD_G: // Gradient Descent by group
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
    switch(this->p->solver)
    {
        case BKT_CGD: // Conjugate Gradient Descent
        case BKT_GD: // Gradient Descent
        case BKT_BW: // Expectation Maximization (Baum-Welch)
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            return this->A[dt->k][i][j];
            break;
        case BKT_GD_G: // Gradient Descent by group
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
    switch(this->p->solver)
    {
        case BKT_CGD: // Conjugate Gradient Descent
        case BKT_GD: // Gradient Descent
        case BKT_BW: // Expectation Maximization (Baum-Welch)
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            return this->B[dt->k][i][m];
            break;
        case BKT_GD_G: // Gradient Descent by group
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
			x_data[x]->alpha = Calloc(NUMBER*, x_data[x]->ndat);
			for(t=0; t<x_data[x]->ndat; t++)
				x_data[x]->alpha[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->ndat; t++)
				for(i=0; i<nS; i++)
					x_data[x]->alpha[t][i] = 0.0;
		}
		// p_O_param
		x_data[x]->p_O_param = 0.0;
		x_data[x]->loglik = 0.0;
        // c - scaling
		if( x_data[x]->c == NULL ) {
            x_data[x]->c = Calloc(NUMBER, x_data[x]->ndat);
        } else {
			for(t=0; t<x_data[x]->ndat; t++)
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
			x_data[x]->gamma = Calloc(NUMBER*, x_data[x]->ndat);
			for(t=0; t<x_data[x]->ndat; t++)
				x_data[x]->gamma[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->ndat; t++)
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
			x_data[x]->xi = Calloc(NUMBER**, x_data[x]->ndat);
			for(t=0; t<x_data[x]->ndat; t++) {
				x_data[x]->xi[t] = Calloc(NUMBER*, nS);
				for(i=0; i<nS; i++)
					x_data[x]->xi[t][i] = Calloc(NUMBER, nS);
			}
			
		} else {
			for(t=0; t<x_data[x]->ndat; t++)
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
			x_data[x]->beta = Calloc(NUMBER*, x_data[x]->ndat);
			for(t=0; t<x_data[x]->ndat; t++)
				x_data[x]->beta[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->ndat; t++)
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
		for(t=0; t<x_data[x]->ndat; t++) {
			o = x_data[x]->obs[t];
			if(t==0) { // it's alpha(1,i)
                // compute \alpha_1(i) = \pi_i b_i(o_1)
				for(i=0; i<nS; i++) {
					x_data[x]->alpha[t][i] = getPI(x_data[x],i) * getB(x_data[x],i,o);
                    x_data[x]->c[t] += x_data[x]->alpha[t][i];
                }
			} else { // it's alpha(t,i)
				// compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
				for(i=0; i<nS; i++) {
					for(j=0; j<nS; j++) {
						x_data[x]->alpha[t][i] += x_data[x]->alpha[t-1][j] * getA(x_data[x],j,i);
					}
					x_data[x]->alpha[t][i] *= getB(x_data[x],i,o);
                    x_data[x]->c[t] += x_data[x]->alpha[t][i];
				}
			}
            // scale \alpha_{t}(i) - same for t=1 or otherwise
            x_data[x]->c[t] = 1/x_data[x]->c[t];//safe0num();
			// compute elements of p(O|param) as a sum of alpha's of last observations in sequences
			if( t==(x_data[x]->ndat-1) ) {
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
		for(t=(NDAT)(x_data[x]->ndat)-1; t>=0; t--) {
			if( t==(x_data[x]->ndat-1) ) { // last \beta
				// \beta_T(i) = 1
				for(i=0; i<nS; i++)
					x_data[x]->beta[t][i] = 1;// x_data[x]->c[t]; // was 1
			} else {
				// \beta_t(i) = \sum_{j=1}^N{beta_{t+1}(j) a_{ij} b_j(o_{t+1})}
				o = x_data[x]->obs[t+1]; // next observation
				for(i=0; i<nS; i++) {
					for(j=0; j<nS; j++)
						x_data[x]->beta[t][i] += x_data[x]->beta[t+1][j] * getA(x_data[x],i,j) * getB(x_data[x],j,o);
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
		for(t=0; t<(x_data[x]->ndat-1); t++) { // -1 is important
			o_tp1 = x_data[x]->obs[t+1];
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
		for(t=0; t<(x_data[x]->ndat-1); t++) { // -1 is important
			for(i=0; i<nS; i++)
				for(j=0; j<nS; j++) {
					x_data[x]->gamma[t][i] += x_data[x]->xi[t][i][j];
                }
		} // for all observations within skill-group
	} // for all groups in skill
}

void HMMProblem::computeGradients(NCAT xndat, struct data** x_data, FitBit *fb){//,  NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB)
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
//	toZero1DNumber(a_gradPI, nS);
//	toZero2DNumber(a_gradA , nS, nS);
//	toZero2DNumber(a_gradB , nS, nO);
    fb->toZero(FBS_GRAD);
    
	computeAlphaAndPOParam(xndat, x_data);
	computeBeta(xndat, x_data);
    //    computeAlphaAndPOParam(xndat, x_data, this->PI[x_data[0]->k], this->A[x_data[0]->k], this->B[x_data[0]->k], nS);
    //	computeBeta(xndat, x_data, this->A[x_data[0]->k], this->B[x_data[0]->k], nS);
	NCAT x;
	NDAT t;
	NPAR i, j, o;
    
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// Mask for gradPI is handled differently
		// Gradient with respect to PI
		t = 0;
		o = x_data[x]->obs[t];
		for(i=0; i<nS /*&& fitparam[0]>0*/; i++) {
			fb->gradPI[i] -= x_data[x]->beta[t][i] * getB(x_data[x],i,o) / safe0num(x_data[x]->p_O_param);
            //			a_gradPI[i] -= x_data[x]->beta[t][i] * this->B[x_data[x]->k][i][o] / safe0num(x_data[x]->p_O_param);
        }
        
		for(t=0; t<x_data[x]->ndat; t++) {
			o = x_data[x]->obs[t];
			// Gradient with respect to A
			// \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
			if( t>0 ) {
				for(i=0; i<nS /*&& fitparam[1]>0*/; i++)
					for(j=0; j<nS; j++) {
						fb->gradA[i][j] -= x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                        //						a_gradA[i][j] -= x_data[x]->beta[t][j] * this->B[x_data[x]->k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                    }
				
			}// if not first obs in sequence
			// Gradient with respect to B
			for(i=0; i<nS /*&& fitparam[2]>0*/; i++) {
				fb->gradB[i][o] -= x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * getB(x_data[x],i,o));
                //				a_gradB[i][o] -= x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * this->B[x_data[x]->k][i][o]);
            }
		} // for all observations within skill-group
	} // for all groups in skill
} // computeGradients()

void HMMProblem::toFile(const char *filename) {
    switch(this->p->solver)
    {
        case BKT_CGD: // Conjugate Gradient Descent
        case BKT_GD: // Gradient Descent
        case BKT_BW: // Expectation Maximization (Baum-Welch)
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            toFileSkill(filename);
            break;
        case BKT_GD_G: // pLo - per stident, pT,pS,pG - per skill - alternating
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
    NUMBER *local_pred_inner = init1DNumber(this->p->nO);
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

void HMMProblem::predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill) {
	NDAT t;
	NCAT g, k;
	NPAR i, j, m, o, isTarget;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
	NUMBER *local_pred = init1DNumber(nO); // local prediction
	char local_know[1024];
	NUMBER pLe[nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(nG, nK, nS);
    NUMBER ll = 0.0, rmse = 0.0, rmse_no_null = 0.0;
    NUMBER p;
    FILE *fid; // file for storing prediction should that be necessary
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
	
	for(g=0; g<nG; g++)
		for(k=0; k<nK; k++) {
            dt->k = k;
            dt->g = g;
			for(i=0; i<nO; i++)
                group_skill_map[g][k][i] =  getPI(dt,i);//PI[i];
		}
	
	for(t=0; t<this->p->N; t++) {
		o = dat_obs->get(t);//[t];
		g = dat_group->get(t);//[t];
        dt->g = g;
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
            isTarget = this->null_skill_obs==o;
            rmse += pow(isTarget - this->null_skill_obs_prob,2);
            ll -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1 - this->null_skill_obs_prob);
            if(this->p->predictions>0) // write predictions file if it was opened
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
            // update p(L)
            pLe_denom = 0.0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++) pLe_denom += group_skill_map[g][k][i] * getB(dt,i,o);//B[i][o];
            for(i=0; i<nS; i++) pLe[i] = group_skill_map[g][k][i] * getB(dt,i,o)/*B[i][o]*/ / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) group_skill_map[g][k][i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    group_skill_map[g][k][j] += pLe[i] * getA(dt,i,j);//A[i][j];
        }
        local_know[0] = 0;
        for(int l=0; l<n; l++)
            sprintf(local_know,"%s%s%10.8f",local_know,(strlen(local_know)>0)?",":"",group_skill_map[g][ ar[l] ][0]);
        if(this->p->predictions>0) { // write predictions file if it was opened
            for(m=0; m<nO; m++)
                fprintf(fid,"%10.8f%s",local_pred[m],(m<(nO-1))?"\t":"\t");
            fprintf(fid,"%s\n",local_know);
            //            for(i=0; i<nS; i++)
            //                fprintf(fid,"%10.8f%s",local_know[i],(i<(nS-1))?"\t":"\n");
        }
        rmse += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
        p = safe01num(local_pred[this->p->metrics_target_obs]);
        ll -= safelog(  p)*   isTarget  +  safelog(1-p)*(1-isTarget);
	} // for all data
    delete(dt);
	free(local_pred);
    //	free(local_pred_inner);
    free3DNumber(group_skill_map, nG, nK);
    rmse = sqrt(rmse / this->p->N);
    rmse_no_null = sqrt(rmse_no_null / (this->p->N - this->p->N_null));
    metrics[0] = ll;
    metrics[1] = 2*nK*4 + 2*ll;
    metrics[2] = nK*4*safelog(this->p->N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    
    if(this->p->predictions>0) // close predictions file if it was opened
        fclose(fid);
}

void HMMProblem::computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data) {
    NCAT x;
    NPAR o,m,i,j, isTarget;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    NDAT t;
    NUMBER *local_pred = init1DNumber(nO);
    NUMBER *pL = init1DNumber(nS);
    NUMBER *pLe = init1DNumber(nS);
    NUMBER pLe_denom = 0;
    NUMBER prob;
    NDAT N = 0;
    for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
        N += x_data[x]->ndat;
        for(i=0;i<nO; i++) pL[i] = safe01num(getPI(x_data[x],i)); // /*safe01num*/(a_PI[i]); // init pL
        for(t=0; t<x_data[x]->ndat; t++) { // for all
            o = x_data[x]->obs[t];//[t];
            isTarget = (this->p->metrics_target_obs == o);
            // produce prediction
            for(m=0; m<nO; m++) local_pred[m] = 0.0; // init pCorr
            for(m=0; m<nO; m++)
                for(i=0; i<nS; i++)
                    local_pred[m] += pL[i] * getB(x_data[x],i,m);//a_B[i][m];
            prob = safe01num(local_pred[this->p->metrics_target_obs]);
            loglik_rmse[0] -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
            //            if( loglik_rmse[0]!=loglik_rmse[0] ) {
            //                int z = 0;
            //            }
            loglik_rmse[1] += pow(isTarget - prob, 2);
            loglik_rmse[2] += pow(isTarget - prob, 2); // for RMSE without null skill
            // update p(L)
            pLe_denom = 0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++)
                pLe_denom += (pL[i]) * getB(x_data[x],i,o);//a_B[i][o];
            for(i=0; i<nS; i++)
                pLe[i] = pL[i] * getB(x_data[x],i,o)/*a_B[i][o]*/ / pLe_denom;
            //            projectsimplex(pLe, param->nS);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) pL[i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    pL[j] += safe01num(pLe[i] * getA(x_data[x],i,j));//a_A[i][j]);
        } // for all data
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
    if(!keep_SE) loglik_rmse[2] = sqrt(loglik_rmse[2]/N);
    free(local_pred);
    free(pLe);
    free(pL);
}

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
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*do RMSE*/);
    switch(this->p->solver)
    {
        case BKT_CGD: // Conjugate Gradient Descent
            loglik_rmse[0] += ConjugateGradientDescent(true /*by skill*/);
            break;
        case BKT_GD: // Gradient Descent
            loglik_rmse[0] += GradientDescent(true /*by skill*/);
            break;
        case BKT_BW: // Expectation Maximization (Baum-Welch)
            loglik_rmse[0] += BaumWelchSkill();
            break;
        case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
            GradientDescent(true /*by skill*/);
            loglik_rmse[0] += BaumWelchSkill();
            break;
        case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            BaumWelchSkill();
            loglik_rmse[0] += GradientDescent(true /*by skill*/);
            break;
        case 6: // pLo - per stident, pT,pS,pG - per skill - alternating
            loglik_rmse[0] += GradientDescent(false /*by group*/);
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

void HMMProblem::computeMetrics(NUMBER* metrics) {
    computeLogLikRMSENullSkill(metrics, true /* only SE*/);
    // despite cycling on k-skill, would work for all
    for(NCAT k=0; k<this->p->nK; k++)
        computeLogLikRMSE(metrics, true /* only SE*/, this->p->k_numg[k], this->p->k_g_data[k]);
    metrics[3] = metrics[1]; // move Squared Errors from position 2
    metrics[3] = sqrt(metrics[3]/this->p->N);  // convert SE to RMSE
    metrics[4] = metrics[2]; // move Squared Errors from position 2
    metrics[4] = sqrt(metrics[4]/(this->p->N - this->p->N_null));  // convert SE to RMSE
    metrics[1] = 2*this->n_params + 2*metrics[0]/*loglik*/;  // AIC
    metrics[2] = this->n_params*safelog(this->p->N) + 2*metrics[0]/*loglik*/;  // BIC
}

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
    NPAR isTarget;
    for(g=0; g<this->p->n_null_skill_group; g++) {
        dat = &this->p->null_skills[g];
        if(dat->cnt != 0)
            continue; // observe block
        count_all_null_skill += dat->ndat;
        for(t=0; t<dat->ndat; t++)
            this->null_obs_ratio[ dat->obs[t] ]++;
    }
    // produce means
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    for(NPAR o=0; o<this->p->nO; o++) {
        this->null_obs_ratio[o] /= count_all_null_skill;
        if( this->null_obs_ratio[o] > this->null_skill_obs_prob ) {
            this->null_skill_obs_prob = this->null_obs_ratio[o];
            this->null_skill_obs = o;
        }
    }
    this->null_skill_obs_prob = safe01num(this->null_skill_obs_prob); // safety for logging
    // compute loglik
    NDAT N = 0;
    for(g=0; g<this->p->n_null_skill_group; g++) {
        dat = &this->p->null_skills[g];
        if(dat->cnt != 0)
            continue; // observe block
        for(t=0; t<dat->ndat; t++) {
            N += dat->ndat;
            isTarget = dat->obs[t] == this->null_skill_obs;
            loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
            loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
        }
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
}

void HMMProblem::computeLogLikRMSENullSkill(NUMBER* loglik_rmse, bool keep_SE) {
    // compute loglik
    NDAT N = 0;
    NPAR isTarget;
    for(NCAT g=0; g<this->p->n_null_skill_group; g++) {
        struct data *dat = &this->p->null_skills[g];
        N += dat->ndat;
        for(NDAT t=0; t<dat->ndat; t++) {
            isTarget = dat->obs[t] == this->null_skill_obs;
            loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
            loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
        }
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/N);
}

void HMMProblem::init3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO) {
    PI = init1DNumber(nS);
    A  = init2DNumber(nS, nS);
    B  = init2DNumber(nS, nO);
}

void HMMProblem::toZero3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO) {
    toZero1DNumber(PI, nS);
    toZero2DNumber(A,  nS, nS);
    toZero2DNumber(B,  nS, nO);
}

void HMMProblem::cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO) {
    cpy1DNumber(soursePI, targetPI, nS);
    cpy2DNumber(sourseA,  targetA,  nS, nS);
    cpy2DNumber(sourseB,  targetB,  nS, nO);
}

void HMMProblem::free3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS) {
    free(PI);
    free2DNumber(A, nS);
    free2DNumber(B, nS);
    PI = NULL;
    A  = NULL;
    B  = NULL;
}

NUMBER HMMProblem::GradientDescent(bool bySkill) {
	NCAT x;
    NUMBER loglik = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG, nX;
    
    if(bySkill)
        nX = nK;
    else
        nX = nG;
    
    FitBit *fb = new FitBit(this->p);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
//        fb->init(FBS_GRADsum);
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));
		while( !conv && iter<=this->p->maxiter ) {
//			if(iter>1)
//                fb->toZero(FBS_GRADsum);
			// add gradients
//			for(k=0; k<nK; k++) {
                computeGradients(this->p->ndata, this->p->k_data, fb);//a_gradPI, a_gradA, a_gradB);
                if(iter==1)
                    pO0 = HMMProblem::getSumLogPOPara(this->p->ndata, this->p->k_data);
//                fb->add(FBS_GRAD, FBS_GRADsum);
//			}
            fb->copy(FBS_PAR, FBS_PARm1);
			// make step
//            fb->copy(FBS_GRADsum, FBS_GRAD);
            doLinearStep(this->p->ndata, this->p->k_data, fb, 0/*copy KC0*/); // step for linked skill 0
			for(x=1; x<nX; x++) { // copy the rest
                NUMBER *aPI = this->getPI(x);
                NUMBER **aA = this->getA(x);
                NUMBER **aB = this->getB(x);
                cpy3Params(fb->PI, fb->A, fb->B, aPI, aA, aB, nS, nO);
            }
            conv = fb->checkConvergence(this->p, conv_flags);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
//                for(k=0; k<nK; k++) {
                    //                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, nS);
                    //                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                    computeAlphaAndPOParam(this->p->ndata, this->p->k_data);
                    pO = HMMProblem::getSumLogPOPara(this->p->ndata, this->p->k_data);
//                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
        RecycleFitData(this->p->ndata, this->p->k_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
//        free3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS);
//        fb->destroy(FBS_GRADsum);
	}
	
	//
	// Main fit
	//
    //	for(k=0; k<nK; k++) this->p->mask_skill[k] = false; // mask with k'th==true to be computed
    
	for(x=0; x<nX; x++) {  // for(k=218; k<219; k++) { //
        NCAT xndat;
        struct data** x_data;
        if(bySkill) {
            xndat = this->p->k_numg[x];
            x_data = this->p->k_g_data[x];
        } else {
            xndat = this->p->g_numk[x];
            x_data = this->p->g_k_data[x];
        }
		
		conv = 0; // converged
		iter = 1; // iteration count
		pO0 = 0.0;
        pO= 0.0;
		
//		PI = this->getPI(k);         // pointers stay same through fitting
//		A  = this->getA(k);
//		B  = this->getB(k);
        fb->linkPar( this->getPI(x), this->getA(x), this->getB(x));
		
		while( !conv && iter<=this->p->maxiter ) {
			computeGradients(xndat, x_data, fb);// a_gradPI, a_gradA, a_gradB);
			if(iter==1) {
                //                computeAlphaAndPOParam(xndat, x_data);
                //				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                //                computeAlphaAndPOParam(xndat, x_data, PI, A, B, nS);
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
//            cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
            fb->copy(FBS_PAR, FBS_PARm1);
			
            doLinearStep(xndat, x_data, fb, -1/*co copy*/);//PI, A, B, a_gradPI, a_gradA, a_gradB);
            
			// check convergence
//			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
            conv = fb->checkConvergence(this->p, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, nS);
                computeAlphaAndPOParam(xndat, x_data);
                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                //                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                loglik += pO*(pO>0);
                if(!this->p->quiet)
                    printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",x,iter,pO0,pO,conv);
			}
			iter ++;
		} // main solver loop
        RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
		// recycle
	} // for all skills
    
    //    // report pO per group
    //    for(NCAT g=0; g<nG; g++) {
    //        //        computeAlphaAndPOParamG(g);
    //        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
    //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        printf("group %3d p(O|param)= %15.7f\n",g,pO);
    //    }
//    free3Params(PI_m1, A_m1, B_m1, nS);
//    free3Params(a_gradPI, a_gradA, a_gradB, nS);
    delete fb; // that takes care of *m1, and grad*
    return loglik;
}

NUMBER HMMProblem::ConjugateGradientDescent(bool bySkill) {
	NCAT x;
    NUMBER loglik = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG, nX;
    if(bySkill)
        nX = nK;
    else
        nX = nG;

//	NUMBER *PI, **A, **B; // just pointer
//	NUMBER *PI_m1, **A_m1, **B_m1; // just pointer
//	NUMBER *a_gradPI, **a_gradA, **a_gradB;
//	NUMBER *a_gradPI_m1, **a_gradA_m1, **a_gradB_m1;
//	NUMBER *a_dirPI_m1, **a_dirA_m1, **a_dirB_m1;
//    init3Params(PI_m1, A_m1, B_m1, nS, nO);
//    init3Params(a_gradPI, a_gradA, a_gradB, nS, nO);
//    init3Params(a_gradPI_m1, a_gradA_m1, a_gradB_m1, nS, nO);
//    init3Params(a_dirPI_m1, a_dirA_m1, a_dirB_m1, nS, nO);
    FitBit *fb = new FitBit(this->p);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    fb->init(FBS_GRADm1);
    fb->init(FBS_DIRm1);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
//        NUMBER *a_gradPI_sum, **a_gradA_sum, **a_gradB_sum;
//        NUMBER *a_gradPI_sum_m1, **a_gradA_sum_m1, **a_gradB_sum_m1;
//        init3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS, nO);
//        init3Params(a_gradPI_sum_m1, a_gradA_sum_m1, a_gradB_sum_m1, nS, nO);
//        fb->init(FBS_GRADsum);
//        fb->init(FBS_GRADsumm1);
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
//        PI = this->getPI(0);
//        A = this->getA(0);
//        B = this->getB(0);
        fb->linkPar(this->getPI(0), this->getA(0), this->getB(0));
		while( !conv && iter<=this->p->maxiter ) {
//			if(iter>1) {
////                toZero3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS, nO);
//                fb->toZero(FBS_GRADsum);
//			}
			// add gradients
//			for(k=0; k<nK; k++) {
//            computeGradients(this->p->k_numg[k], this->p->k_g_data[k], fb);// a_gradPI, a_gradA, a_gradB);
                computeGradients(this->p->ndata, this->p->k_data, fb);// a_gradPI, a_gradA, a_gradB);
                if(iter==1)
                    pO0 = HMMProblem::getSumLogPOPara(this->p->ndata, this->p->k_data);
//				add1DNumbersWeighted(a_gradPI, a_gradPI_sum, nS, 1.0);
//				add2DNumbersWeighted(a_gradA,  a_gradA_sum,  nS, nS, 1.0);
//				add2DNumbersWeighted(a_gradB,  a_gradB_sum,  nS, nO, 1.0);
//                fb->add(FBS_GRAD, FBS_GRADsum);
//			}
			// copy old SAVED! values for params, just for skill #0 is enough
//            cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
            fb->copy(FBS_PAR, FBS_PARm1);
			
			// make step
//			for(k=0; k<nK; k++) {
                //                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, a_gradPI, a_gradA, a_gradB);
//                NCAT xndat = this->p->k_numg[k];
//                struct data** x_data = this->p->k_g_data[k];
                if(iter==1)
                    doLinearStep(this->p->ndata, this->p->k_data, fb, 0/*copy KC0*/);//PI, A, B, a_gradPI_sum, a_gradA_sum, a_gradB_sum);
                else {
                    doConjugateLinearStep(this->p->ndata, this->p->k_data, fb, 0/*copy KC0*/);//PI, A, B, a_gradPI_sum_m1, a_gradA_sum_m1, a_gradB_sum_m1, a_gradPI_sum, a_gradA_sum, a_gradB_sum, a_dirPI_m1, a_dirA_m1, a_dirB_m1);
                }
//            }
			// check convergence, on any skill, e.g. #0
//			conv = checkConvergence(this->getPI(k), this->getA(k), this->getB(k), PI_m1, A_m1, B_m1, conv_flags);
			conv = fb->checkConvergence(this->p, conv_flags);
            
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
//                for(k=0; k<nK; k++) {
//                    NCAT xndat = this->p->k_numg[k];
//                    struct data** x_data = this->p->k_g_data[k];
                    //                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, nS);
                    //                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                    computeAlphaAndPOParam(this->p->ndata, this->p->k_data);
                    pO = HMMProblem::getSumLogPOPara(this->p->ndata, this->p->k_data);
//                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			} else {
                iter ++;
//                cpy3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, a_gradPI_sum_m1, a_gradA_sum_m1, a_gradB_sum_m1, nS, nO);
                fb->copy(FBS_GRAD, FBS_GRADm1);
            }
			iter ++;
		}// single skill loop
//        free3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS);
//        free3Params(a_gradPI_sum_m1, a_gradA_sum_m1, a_gradB_sum_m1, nS);
//        fb->destroy(FBS_GRADsum);
//        fb->destroy(FBS_GRADsumm1);
	}
	
	//
	// Main fit
	//
    //	for(k=0; k<nK; k++) this->p->mask_skill[k] = false; // mask with k'th==true to be computed
    
	for(x=0; x<nK; x++) {  // for(k=218; k<219; k++) { //
        NCAT xndat;
        struct data** x_data;
        if(bySkill) {
            xndat = this->p->k_numg[x];
            x_data = this->p->k_g_data[x];
        } else {
            xndat = this->p->g_numk[x];
            x_data = this->p->g_k_data[x];
        }
		
		conv = 0; // converged
		iter = 1; // iteration count
		pO0 = 0.0;
        pO= 0.0;
		
//		PI = this->getPI(k);         // pointers stay same through fitting
//		A  = this->getA(k);
//		B  = this->getB(k);
        fb->linkPar(this->getPI(x), this->getA(x), this->getB(x));
		
		while( !conv && iter<=this->p->maxiter ) {
			computeGradients(xndat, x_data, fb);// a_gradPI, a_gradA, a_gradB);
			if(iter==1) {
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
//            cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
            fb->copy(FBS_PAR, FBS_PARm1);
            
            if(iter==1)
                doLinearStep(xndat, x_data, fb, -1/*no copy*/);//PI, A, B, a_gradPI, a_gradA, a_gradB);
            else {
                doConjugateLinearStep(xndat, x_data, fb, -1/*no copy*/); //PI, A, B, a_gradPI_m1, a_gradA_m1, a_gradB_m1, a_gradPI, a_gradA, a_gradB, a_dirPI_m1, a_dirA_m1, a_dirB_m1);
            }
            
			// check convergence
//			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
            conv = fb->checkConvergence(this->p, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, nS);
                computeAlphaAndPOParam(xndat, x_data);
                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                //                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                loglik += pO*(pO>0);
                if(!this->p->quiet)
                    printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",x,iter,pO0,pO,conv);
			} else {
                iter ++;
//                cpy3Params(a_gradPI, a_gradA, a_gradB, a_gradPI_m1, a_gradA_m1, a_gradB_m1, nS, nO);
                fb->copy(FBS_GRAD, FBS_GRADm1);
            }
		} // main solver loop
        RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
		// recycle
	} // for all skills
    
    //    // report pO per group
    //    for(NCAT g=0; g<this->p->nG; g++) {
    //        //        computeAlphaAndPOParamG(g);
    //        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
    //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        printf("group %3d p(O|param)= %15.7f\n",g,pO);
    //    }
    
//    free3Params(PI_m1, A_m1, B_m1, nS);
//    free3Params(a_gradPI, a_gradA, a_gradB, nS);
//    free3Params(a_gradPI_m1, a_gradA_m1, a_gradB_m1, nS);
//    free3Params(a_dirPI_m1, a_dirA_m1, a_dirB_m1, nS);
    delete fb; // takes care of everything
    return loglik;
}

NUMBER HMMProblem::BaumWelchSkill() {
	NCAT k;
    NUMBER loglik = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK;
	
    FitBit *fb = new FitBit(this->p);
    fb->init(FBS_PARm1);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
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
                //				hmm->zeroLabelsK(k); // reset blocking labels // THIS IS NOT DONE HERE
//                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, this->p->nS);
//				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                computeAlphaAndPOParam(xndat, x_data);
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
            fb->copy(FBS_PAR, FBS_PARm1);
			
            //			hmm->zeroLabelsK(k); // reset blocking labels // THIS IS NOT DONE HERE
			doBaumWelchStep(xndat, this->p->k_g_data[k], fb);// PI, A, B);
            
            
			// check convergence
//			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
            conv = fb->checkConvergence(this->p, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
//                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, this->p->nS);
//                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
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
//	free(PI_m1);
//	free2DNumber(A_m1, this->p->nS);
//	free2DNumber(B_m1, this->p->nS);
    delete fb;
    return loglik;
}

NUMBER HMMProblem::doLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy) {//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK;
	// first scale down gradients
	doLog10Scale1DGentle(fb->gradPI, fb->PI, nS);
	doLog10Scale2DGentle(fb->gradA,  fb->A,  nS, nS);
	doLog10Scale2DGentle(fb->gradB,  fb->B,  nS, nO);
	
    NUMBER *PI_cpy, ** A_cpy, ** B_cpy;
    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
    //	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
    cpy3Params(fb->PI, fb->A, fb->B, PI_cpy, A_cpy, B_cpy, nS, nO);
//	cpy1DNumber(fb->PI, PI_cpy, nS); // save copy
//	cpy2DNumber(fb->A,  A_cpy,  nS, nS); // save copy
//	cpy2DNumber(fb->B,  B_cpy,  nS, nO); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++)
	{
		p_k_by_neg_p_k -= fb->gradPI[i]*fb->gradPI[i];
		for(j=0; j<nS; j++) p_k_by_neg_p_k -= fb->gradA[i][j]*fb->gradA[i][j];
		for(m=0; m<nO; m++) p_k_by_neg_p_k -= fb->gradB[i][m]*fb->gradB[i][m];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			fb->PI[i] = PI_cpy[i] - e * fb->gradPI[i];
			for(j=0; j<nS; j++)
				fb->A[i][j] = A_cpy[i][j] - e * fb->gradA[i][j];
			for(m=0; m<nO; m++)
				fb->B[i][m] = B_cpy[i][m] - e * fb->gradB[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			projectsimplex(fb->PI, nS);
			for(i=0; i<nS; i++) {
				projectsimplex(fb->A[i], nS);
				projectsimplex(fb->B[i], nS);
			}
		} else {
			projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
			for(i=0; i<nS; i++) {
				projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
				projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
			}
		}
        if(copy >= 0) { // copy parameters from position 'copy' to others
            for(NCAT x=0; x<sizes[0]; x++)
                if(x!=copy)
                    cpy1DNumber(fb->PI, this->PI[x], nS);
            for(NCAT x=0; x<sizes[1]; x++)
                if(x!=copy)
                    cpy2DNumber(fb->A, this->A[x], nS, nS);
            for(NCAT x=0; x<sizes[2]; x++)
                if(x!=copy)
                    cpy2DNumber(fb->B, this->B[x], nS, nO);
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
    if(!compliesArmijo) {
        e = 0;
        cpy1DNumber(PI_cpy, fb->PI, nS); // save copy
        cpy2DNumber(A_cpy,  fb->A,  nS, nS); // save copy
        cpy2DNumber(B_cpy,  fb->B,  nS, nO); // save copy
    }
//    RecycleFitData(xndat, x_data, this->p);
	free(PI_cpy);
	free2DNumber(A_cpy, nS);
	free2DNumber(B_cpy, nS);
    return e;
} // doLinearStep

NUMBER HMMProblem::doConjugateLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy) {//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK;
	// first scale down gradients
	doLog10Scale1DGentle(fb->gradPI, fb->PI, nS);
	doLog10Scale2DGentle(fb->gradA,  fb->A,  nS, nS);
	doLog10Scale2DGentle(fb->gradB,  fb->B,  nS, nO);
    
    // compute beta_gradient_direction
    NUMBER beta_grad = 0, beta_grad_num = 0, beta_grad_den = 0;
    
    switch (this->p->solver_settting) {
        case 1: // Fletcher-Reeves
            for(i=0; i<nS; i++)
            {
                beta_grad_num = fb->gradPI  [i]*fb->gradPI   [i];
                beta_grad_den = fb->gradPIm1[i]*fb->gradPIm1[i];
                for(j=0; j<nS; j++) {
                    beta_grad_num = fb->gradA  [i][j]*fb->gradA   [i][j];
                    beta_grad_den = fb->gradAm1[i][j]*fb->gradAm1[i][j];
                }
                for(m=0; m<nO; m++) {
                    beta_grad_num = fb->gradB  [i][m]*fb->gradB  [i][m];
                    beta_grad_den = fb->gradBm1[i][m]*fb->gradBm1[i][m];
                }
            }
            break;
        case 2: // Polak–Ribiere
            for(i=0; i<nS; i++)
            {
                beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                beta_grad_den =  fb->gradPIm1[i]*fb->gradPIm1[i];
                for(j=0; j<nS; j++) {
                    beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                    beta_grad_den =  fb->gradAm1[i][j]*fb->gradAm1[i][j];
                }
                for(m=0; m<nO; m++) {
                    beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                    beta_grad_den =  fb->gradBm1[i][m]*fb->gradBm1[i][m];
                }
            }
            break;
        case 3: // Hestenes-Stiefel
            for(i=0; i<nS; i++)
            {
                beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                beta_grad_den =  fb->dirPIm1[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
                for(j=0; j<nS; j++) {
                    beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                    beta_grad_den =  fb->dirAm1[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
                }
                for(m=0; m<nO; m++) {
                    beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                    beta_grad_den =  fb->dirBm1[i][m]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
                }
            }
            break;
        default:
            fprintf(stderr,"Wrong stepping algorithm (%d) specified for Conjugate Gradient Descent\n",this->p->solver_settting);
            break;
    }
    beta_grad = beta_grad_num / safe0num(beta_grad_den);
    beta_grad = (beta_grad>=0)?beta_grad:0;
	
    // compute new direction (in place of old)
    fb->toZero(FBS_DIRm1);
//    NUMBER *a_dirPI, ** a_dirA, ** a_dirB;
//    init3Params(a_dirPI, a_dirA, a_dirB, nS, nO);
	for(i=0; i<nS; i++)
	{
		fb->dirPIm1[i] = -fb->gradPI[i] + beta_grad * fb->dirPIm1[i];
		for(j=0; j<nS; j++) fb->dirAm1[i][j] = -fb->gradA[i][j] + beta_grad * fb->dirAm1[i][j];
		for(m=0; m<nO; m++) fb->dirBm1[i][m] = -fb->gradB[i][m] + beta_grad * fb->dirBm1[i][m];
	}
	// scale down direction
	doLog10Scale1DGentle(fb->dirPIm1, fb->PI, nS);
	doLog10Scale2DGentle(fb->dirAm1,  fb->A,  nS, nS);
	doLog10Scale2DGentle(fb->dirBm1,  fb->B,  nS, nO);
    
    
    NUMBER *PI_cpy, ** A_cpy, ** B_cpy;
    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
    //	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
//	cpy1DNumber(fb->PI, PI_cpy, nS); // save copy
//	cpy2DNumber(fb->A,  A_cpy,  nS, nS); // save copy
//	cpy2DNumber(fb->B,  B_cpy,  nS, nO); // save copy
    cpy3Params(fb->PI, fb->A, fb->B, PI_cpy, A_cpy, B_cpy, nS, nO);
	// compute p_k * -p_k >>>> now current gradient by current direction
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++)
	{
		p_k_by_neg_p_k = fb->gradPI[i]*fb->dirPIm1[i];
		for(j=0; j<nS; j++) p_k_by_neg_p_k = fb->gradA[i][j]*fb->dirAm1[i][j];
		for(m=0; m<nO; m++) p_k_by_neg_p_k = fb->gradB[i][m]*fb->dirBm1[i][m];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			fb->PI[i] = PI_cpy[i] + e * fb->dirPIm1[i];
			for(j=0; j<nS; j++)
				fb->A[i][j] = A_cpy[i][j] + e * fb->dirAm1[i][j];
			for(m=0; m<nO; m++)
				fb->B[i][m] = B_cpy[i][m] + e * fb->dirBm1[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			projectsimplex(fb->PI, nS);
			for(i=0; i<nS; i++) {
				projectsimplex(fb->A[i], nS);
				projectsimplex(fb->B[i], nS);
			}
		} else {
			projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
			for(i=0; i<nS; i++) {
				projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
				projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
			}
		}
        // copy parameters from position 'copy' to others
        if(copy >= 0) {
            for(NCAT x=0; x<sizes[0]; x++)
                if(x!=copy)
                    cpy1DNumber(fb->PI, this->PI[x], nS);
            for(NCAT x=0; x<sizes[1]; x++)
                if(x!=copy)
                    cpy2DNumber(fb->A, this->A[x], nS, nS);
            for(NCAT x=0; x<sizes[2]; x++)
                if(x!=copy)
                    cpy2DNumber(fb->B, this->B[x], nS, nO);
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
    if(!compliesArmijo) {
        e = 0;
//        cpy1DNumber(PI_cpy, fb->PI, nS); // save copy
//        cpy2DNumber(A_cpy,  fb->A,  nS, nS); // save copy
//        cpy2DNumber(B_cpy,  fb->B,  nS, nO); // save copy
        cpy3Params(PI_cpy, A_cpy, B_cpy, fb->PI, fb->A, fb->B, nS, nO);
    }
//    RecycleFitData(xndat, x_data, this->p);
//	free(PI_cpy);
//	free2DNumber(A_cpy, nS);
//	free2DNumber(B_cpy, nS);
    free3Params(PI_cpy, A_cpy, B_cpy, nS);
    
//    // copy directions - DONE
//    cpy1DNumber(a_dirPI, a_dirPI_m1, nS); // save copy
//    cpy2DNumber(a_dirA,  a_dirA_m1,  nS, nS); // save copy
//    cpy2DNumber(a_dirB,  a_dirB_m1,  nS, nO); // save copy
//    // free directions
//	free(a_dirPI);
//	free2DNumber(a_dirA, nS);
//	free2DNumber(a_dirB, nS);
    return e;
} // doLinearStep

NUMBER HMMProblem::doBarzalaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
	doLog10Scale2DGentle(a_gradB,  a_B,  nS, nO);
    
    // compute s_k_m1
  	NUMBER *s_k_m1_PI = init1DNumber(nS);
	NUMBER **s_k_m1_A = init2DNumber(nS,nS);
	NUMBER **s_k_m1_B = init2DNumber(nS,nS);
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
	free2DNumber(s_k_m1_B, nS);
	free2DNumber(s_k_m1_A, nS);
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
	
	NUMBER * b_PI = init1DNumber(this->p->nS);
	NUMBER ** b_A_num = init2DNumber(this->p->nS, this->p->nS);
	NUMBER ** b_A_den = init2DNumber(this->p->nS, this->p->nS);
	NUMBER ** b_B_num = init2DNumber(this->p->nS, this->p->nO);
	NUMBER ** b_B_den = init2DNumber(this->p->nS, this->p->nO);
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
		
		for(t=0;t<(x_data[x]->ndat-1);t++) {
			o = x_data[x]->obs[t];
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
	free2DNumber(b_A_num, this->p->nS);
	free2DNumber(b_A_den, this->p->nS);
	free2DNumber(b_B_num, this->p->nS);
	free2DNumber(b_B_den, this->p->nS);
}

bool HMMProblem::checkConvergence(NUMBER* PI, NUMBER** A, NUMBER** B, NUMBER* PI_m1, NUMBER** A_m1, NUMBER** B_m1, bool flags[3]) {
	NUMBER critetion = 0;
	for(NPAR i=0; i<this->p->nS; i++)
	{
		if(flags[0]) critetion += /*(param->fitparam[0]>0)?*/pow( PI[i]-PI_m1[i], 2 )/*:0*/;
		for(NPAR j=0; flags[1] && j<this->p->nS; j++) {
			critetion += pow(A[i][j] - A_m1[i][j],2);
		}
		for(NPAR k=0; flags[2] && k<this->p->nO; k++) {
			critetion += pow(B[i][k] - B_m1[i][k],2);
		}
	}
	return sqrt(critetion) < this->p->tol;
}

void HMMProblem::createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space) {
    int target_state = 0; // known
    
	long long elements, q;

	prob.l = this->p->N - this->p->N_null;
	elements = (long long)prob.l * ((long long)this->p->nK*this->p->nK + 1); // +1 more is for ?, but it's there

	prob.bias= 1; // no bias
//	prob.y = Malloc(double,prob.l);
//	prob.x = Malloc(struct feature_node *,prob.l);
//	x_space = Malloc(struct feature_node,elements+prob.l); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
    //	max_index = 0; // max number of columns
    //
    // vvvvv init HMM
    NDAT /*max_index, inst_max_index,*/ t, tidx;
	NCAT g, k;
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1DNumber(this->p->nO); // local prediction
	NUMBER pLe[this->p->nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(this->p->nG, this->p->nK, this->p->nS);
    struct data dt;
	for(g=0; g<this->p->nG; g++) {
        dt.g = g;
		for(k=0; k<this->p->nK; k++) {
            dt.k = k;
			for(i=0; i<this->p->nO; i++)
                group_skill_map[g][k][i] =  getPI(&dt,i);//PI[i];
		}
    }
    // ^^^^^ init HMM
    //
    t = 0;
	q = 0;
    for(tidx=0; tidx<this->p->N; tidx++) {
        // vvvvvvv
        k = this->p->dat_skill->get(tidx);
        if(k<0) // null skill
            continue;
        g = this->p->dat_group->get(tidx);
        dt.k = k;
        dt.g = g;
        o = this->p->dat_obs->get(tidx); //dat_obs->get(t);//[t];
        // ^^^^^^^
//        fprintf(stdout,"..%d:ll,",tidx);
        //		inst_max_index = 0; // strtol gives 0 if wrong format
        //		readline(fp); // not needed
        prob.x[t] = &x_space[q];
        //		label = strtok(line," \t\n"); // not needed
        //		if(label == NULL) // empty line  // not needed
        //			exit_input_error(i+1);       // not needed
        
        prob.y[t] = (double)o;//dat_obs->get(t); // strtod(label,&  );
        //		if(endptr == label || *endptr != '\0')
        //			exit_input_error(i+1);
        for(NCAT r=0; r<(this->p->nK); r++) //while(1)
        {
            // k*this->p->nK - shift, r position, 1 - idexes start with 1
            NUMBER value = ( group_skill_map[g][r][target_state] );
            if( value!= 0 ) {
                x_space[q].index = k*(this->p->nK+1) + r+1; // idexes are 1-starting //(int) strtol(idx,&endptr,10);
                x_space[q].value = logit(value); // this student, all r skills, just the target state // strtod(val,&endptr);
            }
            ++q;
        }
        ++t;
        //		if(inst_max_index > max_index) // we know the count
        //			max_index = inst_max_index;
        if(prob.bias >= 0) // leave it, let bias variable control it
            x_space[q++].value = prob.bias; // copy bias and step over
        x_space[q++].index = -1; // set index to -1???
        //
        //
        // vvvvvv make an update to the HMM
//        fprintf(stdout,"bkt");
        // deal with null skill
        if(k<0) continue;// if no skill label
        isTarget = this->p->metrics_target_obs == o;
        for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
        // produce prediction and copy to result
        for(m=0; m<this->p->nO; m++)
            for(i=0; i<this->p->nS; i++)
                local_pred[m] += group_skill_map[g][k][i] * getB(&dt,i,m);//B[i][m];
        // update p(L)
        pLe_denom = 0.0;
        // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
        for(i=0; i<this->p->nS; i++) pLe_denom += group_skill_map[g][k][i] * getB(&dt,i,o);//B[i][o];
        for(i=0; i<this->p->nS; i++) pLe[i] = group_skill_map[g][k][i] * getB(&dt,i,o)/*B[i][o]*/ / safe0num(pLe_denom);
        // 2. L = (pLe'*A)';
        for(i=0; i<this->p->nS; i++) group_skill_map[g][k][i] = 0.0;
        for(j=0; j<this->p->nS; j++)
            for(i=0; i<this->p->nS; i++)
                group_skill_map[g][k][j] += pLe[i] * getA(&dt,i,j);//A[i][j];
        // ^^^^^^ make an update to the HMM
        //
        //
//        fprintf(stdout,";\n");
    } // for all t in {G,K}
    //
    // vvvvv Recycle prediction stuff
	free(local_pred);
    free3DNumber(group_skill_map, this->p->nG, this->p->nK);
    // ^^^^^ Recycle prediction stuff
    //
	if(prob.bias >= 0) // taken care of preemptively
	{
        prob.n = (long long)this->p->nK*this->p->nK  + 1;
		for(t=1;t<prob.l;t++)
			(prob.x[t]-2)->index = prob.n;
		x_space[q-2].index = prob.n;
	}
	else
        prob.n = (long long)this->p->nK*this->p->nK ;
    //	fclose(fp);
    //
    // now set up parameter
    //
    param.solver_type = L1R_LR;
    param.C = 1;
    param.eps = 0.01;
    param.weight_label = NULL;
    param.weight = NULL;
    
}

void HMMProblem::createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space, NCAT k) {
    if(k<0 || k>=this->p->nK) {
        fprintf(stderr, "KC specified is out of range!\n");
        exit(1);
    }
    
    int target_obs = 0; // known
    
    NDAT t, tidx;
	long long q;
    fprintf(stdout,"prob.l=%d\n",prob.l);
	prob.bias= 1; // has bias
    // vvvvv init HMM
	NCAT g, kk;
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1DNumber(this->p->nO); // local prediction
	NUMBER pLe[this->p->nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(this->p->nG, this->p->nK, this->p->nS);
    struct data dt;
	for(g=0; g<this->p->nG; g++) {
        dt.g = g;
		for(kk=0; kk<this->p->nK; kk++) {
            dt.k = kk;
			for(i=0; i<this->p->nO; i++)
                group_skill_map[g][kk][i] =  getPI(&dt,i);//PI[i];
		}
    }
    // ^^^^^ init HMM
    //
	q = 0; // position in linearized array of features
    t = 0; // row of feature matrix
    for(tidx=0; tidx<this->p->N; tidx++) {
        // vvvvvvv
        kk = this->p->dat_skill->get(tidx);
        if(kk<0) continue;
        g = this->p->dat_group->get(tidx);
        o = this->p->dat_obs->get(tidx);
        dt.k = kk;
        dt.g = g;
        // ^^^^^^^
        if(kk==k) { // add this to the data
            prob.x[t] = &x_space[q];
            prob.y[t] = (double)o;//dat_obs->get(t); // strtod(label,&  );
            for(NCAT r=0; r<(this->p->nK); r++) //while(1)
            {
                NUMBER pCorr = 0;
                dt.k = r; // set to regression KC(s)
                for(i=0; i<this->p->nS; i++)
                    pCorr += group_skill_map[g][r][i] * getB(&dt,i,target_obs);//B[i][m];
                if( pCorr!= 0 ) {
                    x_space[q].index = r+1; // idexes are 1-starting //(int) strtol(idx,&endptr,10);
                    x_space[q].value = logit(pCorr); // this student, all r skills, just the target state // strtod(val,&endptr);
                }
                ++q;
            }
            ++t;
            dt.k = kk; //back to main KC
            if(prob.bias >= 0) // leave it, let bias variable control it
                x_space[q++].value = prob.bias; // copy bias and step over
            x_space[q++].index = -1; // set index to -1???
        }
        //
        //
        // vvvvvv make an update to the HMM
        // deal with null skill
        isTarget = this->p->metrics_target_obs == o;
        for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
        // produce prediction and copy to result
        for(m=0; m<this->p->nO; m++)
            for(i=0; i<this->p->nS; i++)
                local_pred[m] += group_skill_map[g][kk][i] * getB(&dt,i,m);//B[i][m];
        // update p(L)
        pLe_denom = 0.0;
        // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
        for(i=0; i<this->p->nS; i++) pLe_denom += group_skill_map[g][kk][i] * getB(&dt,i,o);//B[i][o];
        for(i=0; i<this->p->nS; i++) pLe[i] = group_skill_map[g][kk][i] * getB(&dt,i,o) / safe0num(pLe_denom);
        // 2. L = (pLe'*A)';
        for(i=0; i<this->p->nS; i++) group_skill_map[g][k][i] = 0.0;
        for(j=0; j<this->p->nS; j++)
            for(i=0; i<this->p->nS; i++)
                group_skill_map[g][kk][j] += pLe[i] * getA(&dt,i,j);//A[i][j];
        // ^^^^^^ make an update to the HMM
        //
        //
    } // for all t in {G,K}

    if(prob.bias >= 0) // taken care of preemptively
	{
        prob.n = this->p->nK + 1;
		for(t=1;t<prob.l;t++)
			(prob.x[t]-2)->index = prob.n;
		x_space[q-2].index = prob.n;
	}
	else
        prob.n = this->p->nK;
    //
    // vvvvv Recycle prediction stuff
	free(local_pred);
    free3DNumber(group_skill_map, this->p->nG, this->p->nK);
    // ^^^^^ Recycle prediction stuff
    //
    //	fclose(fp);
    //
    // now set up parameter
    //
    param.solver_type = L2R_LR;
    param.C = 1;
    param.eps = 0.01;
    param.weight_label = NULL;
    param.weight = NULL;
    param.nr_weight = 0;
    
}

void HMMProblem::recycleLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space) {
	destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
}

