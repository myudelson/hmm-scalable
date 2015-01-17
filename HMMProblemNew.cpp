/*
 
 Copyright (c) 2012-2014, Michael (Mikhail) Yudelson
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
    a_PI = init1D<NUMBER>((NDAT)nS);
    a_A  = init2D<NUMBER>((NDAT)nS, (NDAT)nS);
    a_B  = init2D<NUMBER>((NDAT)nS, (NDAT)nO);
    
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
		this->PI = init2D<NUMBER>(this->sizes[0], (NDAT)nS);
		this->A =  init3D<NUMBER>(this->sizes[1], (NDAT)nS, (NDAT)nS);
		this->B =  init3D<NUMBER>(this->sizes[2], (NDAT)nS, (NDAT)nO);
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
	
    
    // populate boundaries
	// populate lb*/ub*
	// *PI
    this->lbPI = init1D<NUMBER>(nS);
    this->lbA  = init2D<NUMBER>(nS, nS);
    this->lbB  = init2D<NUMBER>(nS, nO);
    this->lbPI = init1D<NUMBER>(nS);
    this->lbA  = init2D<NUMBER>(nS, nS);
    this->lbB  = init2D<NUMBER>(nS, nO);
    this->ubPI = init1D<NUMBER>(nS);
    this->ubA  = init2D<NUMBER>(nS, nS);
    this->ubB  = init2D<NUMBER>(nS, nO);
    this->ubPI = init1D<NUMBER>(nS);
    this->ubA  = init2D<NUMBER>(nS, nS);
    this->ubB  = init2D<NUMBER>(nS, nO);
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
    // hmm variables
    this->alpha = NULL;
    this->beta  = NULL;
    this->pL    = NULL;
    this->c     = NULL;
}

HMMProblem::~HMMProblem() {
    destroy();
}

void HMMProblem::destroy() {
    NPAR nS = this->p->nS, NN = this->p->NN;
	// destroy model data
    if(this->null_obs_ratio != NULL) free(this->null_obs_ratio);
	if(this->PI != NULL) free2D<NUMBER>(this->PI, this->sizes[0]);
	if(this->A  != NULL) free3D<NUMBER>(this->A,  this->sizes[1], this->p->nS);
	if(this->B  != NULL) free3D<NUMBER>(this->B,  this->sizes[2], this->p->nS);
	if(this->lbPI!=NULL) free(this->lbPI);
	if(this->ubPI!=NULL) free(this->ubPI);
	if(this->lbA!=NULL) free2D<NUMBER>(this->lbA, nS);
	if(this->ubA!=NULL) free2D<NUMBER>(this->ubA, nS);
	if(this->lbB!=NULL) free2D<NUMBER>(this->lbB, nS);
	if(this->ubB!=NULL) free2D<NUMBER>(this->ubB, nS);
    
    if(this->alpha != NULL) free2D<NUMBER>(this->alpha, NN);
    if(this->beta  != NULL) free2D<NUMBER>(this->beta,  NN);
    if(this->c     != NULL) free(this->c);
    if(this->pL    != NULL) free2D<NUMBER>(this->pL, NN);
    
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

// getters for computing alpha, beta, gamma - also works as get pL :)
NUMBER HMMProblem::getPI(NDAT t, NPAR i) {
    // combination formula for a,b,c
    // \frac{a b c}{ a b c + (1-a)(1-b)(1-c) }
    // the computation is done by keeping numerator and denominator of the innter fraction
    NUMBER prod = 1;   // product of probabilities
    NUMBER prod1m = 1; // product of q minus probabilities
    NUMBER val;
    NDAT ix; 
    // skill components
    if(this->p->multiskill == 0) {
        NCAT k = this->p->dat_skill->get(t);
        val = this->PI[ k ][ i ];
        prod1m *= (1-val);
        prod *= val;
    } else {
        NCAT *ms  = this->p->dat_multiskill->get(t); // multiskills
        data **msq = this->p->dat_multiskill_seq[t]; // multiskill data sequences
        val = 0;
        for(NPAR l=0; l<ms[0]; l++) {
            if( msq[l]->cnt == 0 ) // if we start the sequence. cnt marks the current sequence position, and it's "start"
                val = this->PI[ ms[l+1] ][ i ];
            else { // grab pL of prior step instead
                ix = msq[l]->ix1st + msq[l]->cnt;
                val = this->pL[ ix-1 ][ i ];
            }
            prod1m *= (1-val);
            prod *= val;
        }
    }
    /*
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->PI[dt->id][i];
            break;
//        case STRUCTURE_GROUP:
//            return this->PI[dt->g][i];
//            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    } */
    return prod / safe0num( prod + prod1m );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblem::getA (NDAT t, NPAR i, NPAR j) {
    // combination formula for a,b,c
    // \frac{1}{ 1 + \frac{ (1-a)(1-b)(1-c) }{ a b c } }
    // the computation is done by keeping numerator and denominator of the innter fraction
    NUMBER prod   = 1; // product of probabilities
    NUMBER prod1m = 1; // product of q minus probabilities
    // skill components
    if(this->p->multiskill == 0) {
        NDAT ix = this->p->dat_skill->get(t);
        NUMBER val = this->A[ ix ][i][j];
        prod1m *= (1-val);
        prod   *= val;
    } else {
        NCAT *ms = this->p->dat_multiskill->get(t);
        data **dats = this->p->dat_multiskill_seq[t];
        NUMBER val = 0;
        for(NPAR k=0; k<ms[0]; k++) {
            if( dats[k]->cnt ==0 ) // if this is a first element of sequence, no transition is relevant
                continue;
            val = this->A[ ms[k+1] ][i][j];
            prod1m *= (1-val);
            prod   *= val;
        }
    }
    /*
    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->A[dt->id][i][j];
            break;
//        case STRUCTURE_GROUP:
//            return this->A[dt->g][i][j];
//            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    }*/
    return prod / safe0num( prod + prod1m );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblem::getB (NDAT t, NPAR i, NPAR m) {
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;

    // combination formula for a,b,c
    // \frac{1}{ 1 + \frac{ (1-a)(1-b)(1-c) }{ a b c } }
    // the computation is done by keeping numerator and denominator of the innter fraction
    NUMBER prod   = 1; // product of probabilities
    NUMBER prod1m = 1; // product of q minus probabilities
    // skill components
    if(this->p->multiskill == 0) {
        NDAT ix = this->p->dat_skill->get(t);
        NUMBER val = this->B[ ix ][i][m];
        prod1m *= (1-val);
        prod   *= val;
    } else {
        NCAT *ms = this->p->dat_multiskill->get(t);
        data **dats = this->p->dat_multiskill_seq[t];
        NUMBER val = 0;
        for(NPAR k=0; k<ms[0]; k++) {
            if( dats[k]->cnt ==0 ) // if this is a first element of sequence, no emission is relevant
                continue;
            val = this->B[ ms[k+1] ][i][m];
            prod1m *= (1-val);
            prod   *= val;
        }
    }
    /*

    switch(this->p->structure)
    {
        case STRUCTURE_SKILL:
            return this->B[dt->id][i][m];
            break;
//        case STRUCTURE_GROUP:
//            return this->B[dt->g][i][m];
//            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            exit(1);
            break;
    }*/
    return prod / safe0num( prod + prod1m );
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

// still works in skill-ordered array
NUMBER HMMProblem::getSumLogPOPara(NCAT ix1, NCAT n, struct data* data) {
	NUMBER result = 0.0;
    //	for(NCAT x=0; x<xndat; x++) result += (x_data[x]->tag==0)?x_data[x]->p_O_param:0; << WAS JUST PLANE WRONG (due to massive find replace most surely)
	for(NCAT x=0; x<n; x++) result -= safelog(data[ix1+x].p_O_param); // tag != 0, but p_O_param still summed (it's 0 by default, or pre-computed already)
    // take neg log here
	return result;
}

NUMBER HMMProblem::getSumLogPOPara(NCAT k) {
	NUMBER result = 0.0;
    // tag != 0, but p_O_param still summed (it's 0 by default, or pre-computed already)
    // take neg log here
	for(NCAT q=0; q<this->p->k_nG[k]; q++)
        result -= safelog( this->p->all_seq[ q+this->p->k_ix1stSeq[k] ].p_O_param );
    
	return result;
}


// just init internal variable, including c
void HMMProblem::initAlpha() { //NCAT xndat, struct data** x_data) {
    if(this->alpha != NULL || this->c != NULL || this->pL != NULL) {
        fprintf(stderr,"Error when initing `alpha`, `pL`, or scaling `c` which should be NULL but it is not.\n");
        return;
    }
    this->alpha = init2D<NUMBER>(this->p->NN, this->p->nS);
    this->c     = init1D<NUMBER>(this->p->NN);
    this->pL    = init2D<NUMBER>(this->p->NN, this->p->nS);
    /*
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->alpha == NULL ) {
			x_data[x]->alpha = Calloc(NUMBER*, x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->alpha[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->alpha[t][i] = 0.0;
		}
		// p_O_param
		x_data[x]->p_O_param = 0.0;
//		x_data[x]->n = 0.0;
        // c - scaling
		if( x_data[x]->c == NULL ) {
            x_data[x]->c = Calloc(NUMBER, x_data[x]->n);
        } else {
			for(t=0; t<x_data[x]->n; t++)
                x_data[x]->c[t] = 0.0;
        }
	} // for all groups in skill
    */
}

void HMMProblem::initBeta() {//NCAT xndat, struct data** x_data) {
    if(this->beta != NULL) {
        fprintf(stderr,"Error when initing beta which should be NULL but it is not.\n");
        return;
    }
    this->beta = init2D<NUMBER>(this->p->NN, this->p->nS);
    /*
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// beta
		if( x_data[x]->beta == NULL ) {
			x_data[x]->beta = Calloc(NUMBER*, x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->beta[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->beta[t][i] = 0.0;
		}
	} // for all groups in skill
    */
} // initBeta

void HMMProblem::initGamma() { //(NCAT xndat, struct data** x_data) {
    fprintf(stderr,"initGamma is not implemented yet.\n");
    /*
	NCAT x;
	NDAT t;
	NPAR i, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->gamma == NULL ) {
			x_data[x]->gamma = Calloc(NUMBER*, x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++)
				x_data[x]->gamma[t] = Calloc(NUMBER, nS);
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					x_data[x]->gamma[t][i] = 0.0;
		}
	} // for all groups in skill
    */
}

void HMMProblem::initXi() { //(NCAT xndat, struct data** x_data) {
    fprintf(stderr,"initXi is not implemented yet.\n");
    /*
	NCAT x;
	NDAT t;
	NPAR i,j, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// alpha
		if( x_data[x]->xi == NULL ) {
			x_data[x]->xi = Calloc(NUMBER**, x_data[x]->n);
			for(t=0; t<x_data[x]->n; t++) {
				x_data[x]->xi[t] = Calloc(NUMBER*, nS);
				for(i=0; i<nS; i++)
					x_data[x]->xi[t][i] = Calloc(NUMBER, nS);
			}
			
		} else {
			for(t=0; t<x_data[x]->n; t++)
				for(i=0; i<nS; i++)
					for(j=0; j<nS; j++)
						x_data[x]->xi[t][i][j] = 0.0;
		}
	} // for all groups in skill
     */
}

void HMMProblem::computeAlphaEtAl0() { //NCAT xndat, struct data** x_data) { // fastest, no subfunction
    NDAT ix; // t_m1 - t minus 1 (previous opportunity)
    NPAR i, j, o, n, nS = this->p->nS;
    NCAT *ar, k;
    data **dats;
    NUMBER pLe[nS];// p(L|evidence);
    NUMBER pLe_denom;// denomunator for p(L|evidence)
    bool is_multiskill = p->multiskill != 0;
    for(NDAT t=0; t<p->N; t++) {
        o = this->p->dat_obs->get(t);
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]<0 ) // empty skill
            continue;
        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            // grab index in the variable array
            ix = dats[l]->ix1st + dats[l]->cnt;
            
            // first update the pL (no reason for doing it first, just a preference)
            for(i=0; i<nS; i++)
                this->pL[ix][i] = 0;
            pLe_denom = 0.0;
            // compute prior pL - can be done via getPI(t,i) for actually any t :)
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++) pLe_denom += this->getPI(t,i) * this->getB(t,i,o);
            for(i=0; i<nS; i++) pLe[i] = this->getPI(t,i) * this->getB(t,i,o) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) this->pL[ix][i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    this->pL[ix][j] += pLe[i] * this->getA(t,i,j);
            
            //  compule alpha
            this->c[ix] = 0;
            if(dats[l]->cnt==0) { // first opportunity with the skill
                // compute \alpha_1(i) = \pi_i b_i(o_1)
                for(i=0; i<nS; i++) {
                    this->alpha[ix][i] = getPI(t,i) * ((o<0)?1:getB(t,i,o)); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            } else { // further opportunities
                // compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
                for(i=0; i<nS; i++)
                    this->alpha[ix][i] = 0;
                for(i=0; i<nS; i++) {
                    for(j=0; j<nS; j++) {
                        this->alpha[ix][i] += this->alpha[ix-1][j] * getA(t,j,i);
                    }
                    this->alpha[ix][i] *= ((o<0)?1:getB(t,i,o)); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            }
            // adjust scale
            this->c[ix] = 1/this->c[ix];//safe0num();
            //
            // advance counters
            //
            dats[l]->cnt++;
            // if last row - compute p(O|param)
            if( dats[l]->n == dats[l]->cnt ) {
                dats[l]->p_O_param = 0;
                for(i=0; i<nS; i++)
                    dats[l]->p_O_param += this->alpha[ix][i];
                dats[l]->cnt--; // roll back
            }
        }// for all row skills
    }// for all data
}

void HMMProblem::computeAlphaEtAl0(NCAT target_k, bool doZeroCount) { // for a particular skill
    NDAT ix, t; // t_m1 - t minus 1 (previous opportunity)
    NPAR i, j, o, n, nS = this->p->nS;
    NCAT *ar, k;
    data **dats;
    NUMBER pLe[nS];// p(L|evidence);
    NUMBER pLe_denom;// denomunator for p(L|evidence)
    bool is_multiskill = p->multiskill != 0;
    for(NDAT T=0; T<this->p->k_N[target_k]; T++) {
        t = this->p->k_t[target_k][T];
        o = this->p->dat_obs->get( t );
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]<0 ) // empty skill
            return;
        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            if( k!= target_k) // only compute for target skill
                continue;
            // grab index in the variable array
            ix = dats[l]->ix1st + dats[l]->cnt;
            
            // first update the pL (no reason for doing it first, just a preference)
            for(i=0; i<nS; i++)
                this->pL[ix][i] = 0;
            pLe_denom = 0.0;
            // compute prior pL - can be done via getPI(t,i) for actually any t :)
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++) pLe_denom += this->getPI(t,i) * this->getB(t,i,o);
            for(i=0; i<nS; i++) pLe[i] = this->getPI(t,i) * this->getB(t,i,o) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) this->pL[ix][i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    this->pL[ix][j] += pLe[i] * this->getA(t,i,j);
            
            //  compule alpha
            this->c[ix] = 0;
            if(dats[l]->cnt==0) { // first opportunity with the skill
                // compute \alpha_1(i) = \pi_i b_i(o_1)
                for(i=0; i<nS; i++) {
                    this->alpha[ix][i] = getPI(t,i) * ((o<0)?1:getB(t,i,o)); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            } else { // further opportunities
                // compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
                for(i=0; i<nS; i++)
                    this->alpha[ix][i] = 0;
                for(i=0; i<nS; i++) {
                    for(j=0; j<nS; j++) {
                        this->alpha[ix][i] += this->alpha[ix-1][j] * getA(t,j,i);
                    }
                    this->alpha[ix][i] *= ((o<0)?1:getB(t,i,o)); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            }
            // adjust scale
            this->c[ix] = 1/this->c[ix];//safe0num();
            //
            // advance counters
            //
            dats[l]->cnt++;
            // if last row - compute p(O|param)
            if( dats[l]->n == dats[l]->cnt ) {
                dats[l]->p_O_param = 0;
                for(i=0; i<nS; i++)
                    dats[l]->p_O_param += this->alpha[ix][i];
                dats[l]->cnt--; // roll back
            }
        }// for all row skills
    }// for all data
    // zero counts for this skill
    if(doZeroCount)
        for(NCAT q=0; q<this->p->k_nG[target_k]; q++)
            this->p->all_seq[ q + this->p->k_ix1stSeq[target_k] ].cnt = 0;
}

void HMMProblem::computeAlphaEtAl() { //NCAT xndat, struct data** x_data) { // fastest, no subfunction
//    clock_t t0;
//    NUMBER t1 = 0/*multi single*/, t2 = 0/*update pL*/, t3 = 0/*compute aplpha*/, t4 = 0/*increase counters*/;
    
    NDAT ix; // t_m1 - t minus 1 (previous opportunity)
    NPAR i, j, o, n, m, nS = this->p->nS, nO = this->p->nO;
    NCAT *ar, k;
    data **dats;
    NUMBER pLe[nS];// p(L|evidence);
    NUMBER pLe_denom;// denomunator for p(L|evidence)
    bool is_multiskill = p->multiskill != 0;
    NUMBER *locPI;
    NUMBER **locA;
    NUMBER **locB;
    if(is_multiskill) {
        locPI = init1D<NUMBER>(nS);
        locA  = init2D<NUMBER>(nS, nS);
        locB  = init2D<NUMBER>(nS, nO);
    }
    for(NDAT t=0; t<p->N; t++) {
        o = this->p->dat_obs->get(t);
        //
        // multi or single skill
        //
//        t0 = clock();
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
//        t1 += (NUMBER)(clock()-t0)/(1*CLOCKS_PER_SEC);
        if( ar[0]<0 ) // empty skill
            continue;
        
        if(is_multiskill) {
            for(i=0; i<nS; i++) {
                locPI[i] = this->getPI(t,i);
                for(j=0; j<nS; j++)
                    locA[i][j] = this->getA(t,i,j);
                for(m=0; m<nO; m++)
                    locB[i][m] = this->getB(t,i,m);
            }
        } else {
            locPI = this->PI[ ar[0] ];
            locA  = this->A [ ar[0] ];
            locB  = this->B [ ar[0] ];
        }
        
        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            // grab index in the variable array
            ix = dats[l]->ix1st + dats[l]->cnt;
            
            // first update the pL (no reason for doing it first, just a preference)
//            t0 = clock();
            for(i=0; i<nS; i++)
                this->pL[ix][i] = 0;
            pLe_denom = 0.0;
            // compute prior pL - can be done via getPI(t,i) for actually any t :)
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++) pLe_denom += locPI[i] * ((o<0)?1:locB[i][o]);
            for(i=0; i<nS; i++) pLe[i] = locPI[i] * ((o<0)?1:locB[i][o]) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) this->pL[ix][i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    this->pL[ix][j] += pLe[i] * locA[i][j];
//            t2 += (NUMBER)(clock()-t0)/(1*CLOCKS_PER_SEC);
            //  compule alpha
//            t0 = clock();
            this->c[ix] = 0;
            if(dats[l]->cnt==0) { // first opportunity with the skill
                // compute \alpha_1(i) = \pi_i b_i(o_1)
                for(i=0; i<nS; i++) {
                    this->alpha[ix][i] = locPI[i] * ((o<0)?1:locB[i][o]); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            } else { // further opportunities
                // compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
                for(i=0; i<nS; i++)
                    this->alpha[ix][i] = 0;
                for(i=0; i<nS; i++) {
                    for(j=0; j<nS; j++) {
                        this->alpha[ix][i] += this->alpha[ix-1][j] * locA[j][i];
                    }
                    this->alpha[ix][i] *= ((o<0)?1:locB[i][o]); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            }
            // adjust scale
            this->c[ix] = 1/this->c[ix];//safe0num();
//            t3 += (NUMBER)(clock()-t0)/(1*CLOCKS_PER_SEC);
            //
            // advance counters
            //
//            t0 = clock();
            dats[l]->cnt++;
            // if last row - compute p(O|param)
            if( dats[l]->n == dats[l]->cnt ) {
                dats[l]->p_O_param = 0;
                for(i=0; i<nS; i++)
                    dats[l]->p_O_param += this->alpha[ix][i];
                dats[l]->cnt--; // roll back
            }
//            t4 += (NUMBER)(clock()-t0)/(1*CLOCKS_PER_SEC);
        }// for all row skills
    }// for all data
    if(is_multiskill) {
        free(locPI);
        free2D(locA, nS);
        free2D(locB, nS);
    }
}

void HMMProblem::computeAlphaEtAl(NCAT target_k, bool doZeroCount) { // for a particular skill
    NDAT ix, t; // t_m1 - t minus 1 (previous opportunity)
    NPAR i, j, o, n, m, nS = this->p->nS, nO = this->p->nO;
    NCAT *ar, k;
    data **dats;
    NUMBER pLe[nS];// p(L|evidence);
    NUMBER pLe_denom;// denomunator for p(L|evidence)
    bool is_multiskill = p->multiskill != 0;
    NUMBER *locPI;
    NUMBER **locA;
    NUMBER **locB;
    if(is_multiskill) {
        locPI = init1D<NUMBER>(nS);
        locA  = init2D<NUMBER>(nS, nS);
        locB  = init2D<NUMBER>(nS, nO);
    }
    for(NDAT T=0; T<this->p->k_N[target_k]; T++) {
        t = this->p->k_t[target_k][T];
        o = this->p->dat_obs->get( t );
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]<0 ) // empty skill
            return;
        
        if(is_multiskill) {
            for(i=0; i<nS; i++) {
                locPI[i] = this->getPI(t,i);
                for(j=0; j<nS; j++)
                    locA[i][j] = this->getA(t,i,j);
                for(m=0; m<nO; m++)
                    locB[i][m] = this->getB(t,i,m);
            }
        } else {
            locPI = this->PI[ ar[0] ];
            locA  = this->A [ ar[0] ];
            locB  = this->B [ ar[0] ];
        }

        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            if( k!= target_k) // only compute for target skill
                continue;
            // grab index in the variable array
            ix = dats[l]->ix1st + dats[l]->cnt;
            
            // first update the pL (no reason for doing it first, just a preference)
            for(i=0; i<nS; i++)
                this->pL[ix][i] = 0;
            pLe_denom = 0.0;
            // compute prior pL - can be done via getPI(t,i) for actually any t :)
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++) pLe_denom += locPI[i] * locB[i][o];
            for(i=0; i<nS; i++) pLe[i] = locPI[i] * locB[i][o] / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++) this->pL[ix][i] = 0.0;
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    this->pL[ix][j] += pLe[i] * locA[i][j];
            
            //  compule alpha
            this->c[ix] = 0;
            if(dats[l]->cnt==0) { // first opportunity with the skill
                // compute \alpha_1(i) = \pi_i b_i(o_1)
                for(i=0; i<nS; i++) {
                    this->alpha[ix][i] = locPI[i] * ((o<0)?1:locB[i][o]); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            } else { // further opportunities
                // compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
                for(i=0; i<nS; i++)
                    this->alpha[ix][i] = 0;
                for(i=0; i<nS; i++) {
                    for(j=0; j<nS; j++) {
                        this->alpha[ix][i] += this->alpha[ix-1][j] * locA[j][i];
                    }
                    this->alpha[ix][i] *= ((o<0)?1:locB[i][o]); // if observatiob unknown use 1
                    this->c[ix] += this->alpha[ix][i];
                }
            }
            // adjust scale
            this->c[ix] = 1/this->c[ix];//safe0num();
            //
            // advance counters
            //
            dats[l]->cnt++;
            // if last row - compute p(O|param)
            if( dats[l]->n == dats[l]->cnt ) {
                dats[l]->p_O_param = 0;
                for(i=0; i<nS; i++)
                    dats[l]->p_O_param += this->alpha[ix][i];
                dats[l]->cnt--; // roll back
            }
        }// for all row skills
    }// for all data
    // zero counts for this skill
    if(doZeroCount)
        for(NCAT q=0; q<this->p->k_nG[target_k]; q++)
            this->p->all_seq[ q + this->p->k_ix1stSeq[target_k] ].cnt = 0;
    if(is_multiskill) {
        free(locPI);
        free2D(locA, nS);
        free2D(locB, nS);
    }
}


void HMMProblem::computeBeta() { //NCAT xndat, struct data** x_data) {
    NCAT *ar, k;
    NDAT n, ix, t_p1;
    data **dats;
    NPAR i, j, o, nS = this->p->nS;
    param *p = this->p;
    bool is_multiskill = p->multiskill != 0;
    bool finished = false; // NDAT t is alwasy >=0, use flag to control this
    for(NDAT t=(this->p->N-1); !finished; t--) { // for all data
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]>=0 ) { // empty skill
            //
            // for all skills in question
            //
            for(int l=0; l<n; l++) {
                // if we're still fitting this one
                if(dats[l]->tag == 0) { // not blocked seq
                    k = ar[l];
                    // grab index in the variable array
                    ix = dats[l]->ix1st + dats[l]->cnt;
                    
                    if( dats[l]->cnt==(dats[l]->n-1) ) { // last \beta
                        // \beta_T(i) = 1
                        for(i=0; i<nS; i++)
                            this->beta[ix][i] = 1;// x_data[x]->c[t]; // was 1
                    } else { // not las beta
                        // \beta_t(i) = \sum_{j=1}^N{beta_{t+1}(j) a_{ij} b_j(o_{t+1})}
                        for(i=0; i<nS; i++) // zero it first
                            this->beta[ix][i] = 0;
                        t_p1 = this->p->dat_tix[ ix + 1 ];
                        o = this->p->dat_obs->get(t_p1); // next observation actually
                        for(i=0; i<nS; i++) {
                            for(j=0; j<nS; j++)
                                this->beta[ix][i] += this->beta[ix+1][j] * getA(t,i,j) * ((o<0)?1:getB(t,j,o)); // if observatiob unknown use 1
                            // scale
                            //                    this->beta[ix][i] *= this->c[ix];
                        }
                    } // done with beta
                    //
                    // reduce counters
                    //
                    if(dats[l]->cnt>0) // do not decrease below 0
                        dats[l]->cnt--;
                } // not blocked seq
            }// all skills in question
        }// not nul; skill
        if(t==0) // last row
            finished = true;
    }// for all data
}

void HMMProblem::computeBeta(NCAT target_k) { //NCAT xndat, struct data** x_data) {
    NCAT *ar, k;
    NDAT t, n, ix, t_p1;
    data **dats;
    NPAR i, j, o, nS = this->p->nS;
    param *p = this->p;
    bool is_multiskill = p->multiskill != 0;
    bool finished = false; // NDAT t is alwasy >=0, use flag to control this
    for(NDAT T=(this->p->k_N[target_k]-1); !finished; T--) {  // for all data
        t = this->p->k_t[target_k][T];
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]>=0 ) { // non mpty skill
            //
            // for all skills in question
            //
            for(int l=0; l<n; l++) { // all skills in row
                // if we're still fitting this one
                if(dats[l]->tag != 0)
                    continue;
                k = ar[l];
                if( k == target_k) { // only compute for target skill
                    // grab index in the variable array
                    ix = dats[l]->ix1st + dats[l]->cnt;
                    
                    if( dats[l]->cnt==(dats[l]->n-1) ) { // last \beta
                        // \beta_T(i) = 1
                        for(i=0; i<nS; i++)
                            this->beta[ix][i] = 1;// x_data[x]->c[t]; // was 1
                    } else { // not las beta
                        // \beta_t(i) = \sum_{j=1}^N{beta_{t+1}(j) a_{ij} b_j(o_{t+1})}
                        for(i=0; i<nS; i++) // zero it first
                            this->beta[ix][i] = 0;
                        t_p1 = this->p->dat_tix[ ix + 1 ];
                        o = this->p->dat_obs->get(t_p1); // next observation actually
                        for(i=0; i<nS; i++) {
                            for(j=0; j<nS; j++)
                                this->beta[ix][i] += this->beta[ix+1][j] * getA(t,i,j) * ((o<0)?1:getB(t,j,o)); // if observatiob unknown use 1
                            // scale
                            //                    this->beta[ix][i] *= this->c[ix];
                        }
                    } // done with beta
                    //
                    // reduce counters
                    //
                    if(dats[l]->cnt>0) // do not decrease below 0
                        dats[l]->cnt--;
                } // only compute for target skill
            } // all skills in row
        } // non mpty skill
        if(T==0) // last row
            finished = true;
    }// for all data
}

void HMMProblem::computeXi() { // (NCAT xndat, struct data** x_data){
    fprintf(stderr,"computeXi is not implemented yet.\n");
    /*
	HMMProblem::initXi(xndat, x_data);
	NCAT x;
	NDAT t;
	NPAR i, j, o_tp1, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
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
    */
}

void HMMProblem::computeGamma() { //(NCAT xndat, struct data** x_data) {
    fprintf(stderr,"computeGamma is not implemented yet.\n");
    /*
	HMMProblem::initGamma(xndat, x_data);
	NCAT x;
	NDAT t;
	NPAR i, j, nS = this->p->nS;
	for(x=0; x<xndat; x++) {
		if( x_data[x]->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		for(t=0; t<(x_data[x]->n-1); t++) { // -1 is important
			for(i=0; i<nS; i++)
				for(j=0; j<nS; j++) {
					x_data[x]->gamma[t][i] += x_data[x]->xi[t][i][j];
                }
		} // for all observations within skill-group
	} // for all groups in skill
    */
}

void HMMProblem::setGradPI(struct data* dt, FitBit *fb){
    NPAR i, o;
    NDAT ix = dt->ix1st, t;
    NUMBER combined, separate, mult;
    t = p->dat_tix[ix+0]; // first index only
    o = p->dat_obs->get( t );
    for(i=0; i<p->nS; i++) {
        combined = getPI(t,i); // combined PI_i
        separate = fb->PI[i]; // singled out
        mult = (this->p->multiskill==0)?1:combined * (1-combined)  / safe0num( separate * (1-separate) );
        fb->gradPI[i] -=  mult * this->beta[ix][i] * ((o<0)?1:getB(t,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,separate); // PENALTY
    }
}

void HMMProblem::setGradPI(NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb){
    NPAR i;
    NUMBER combined, separate, mult;
    for(i=0; i<p->nS; i++) {
        if(is_multiskill) {
            combined = getPI(t,i); // combined PI_i
            separate = fb->PI[i]; // singled out
            mult = combined * (1-combined)  / safe0num( separate * (1-separate) );
            fb->gradPI[i] -=  mult * this->beta[ix][i] * ((o<0)?1:getB(t,i,o)) / safe0num(p_O_param) + L2penalty(this->p,separate); // PENALTY
        } else {
            fb->gradPI[i] -=  this->beta[ix][i] * ((o<0)?1:fb->B[i][o]) / safe0num(p_O_param) + L2penalty(this->p, fb->PI[i]); // PENALTY
        }
    }
}

void HMMProblem::setGradA (struct data* dt, FitBit *fb){
    NDAT t, ix;
    NPAR o, i, j;
    NPAR nS = this->p->nS;
    NUMBER combined, separate, mult;
    for(ix=(dt->ix1st+1); ix<(dt->ix1st + dt->n); ix++) { // skip first element, start with 1
        t = this->p->dat_tix[ix];
        o = this->p->dat_obs->get( t );
        for(i=0; i<nS /*&& fitparam[1]>0*/; i++)
            for(j=0; j<nS; j++) {
                combined = getA(t,i,j);
                separate = fb->A[i][j];
                mult = (this->p->multiskill==0)?1:combined * (1-combined)  / safe0num( separate * (1-separate) );
                fb->gradA[i][j] -= mult * this->beta[ix][j] * ((o<0)?1:getB(t,j,o)) * this->alpha[ix-1][i] / safe0num(dt->p_O_param) + L2penalty(this->p,separate); // PENALTY
            }
    }
}

void HMMProblem::setGradA (NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb){
    NPAR i, j;
    NPAR nS = this->p->nS;
    NUMBER combined, separate, mult;
    t = this->p->dat_tix[ix];
    o = this->p->dat_obs->get( t );
    for(i=0; i<nS /*&& fitparam[1]>0*/; i++)
        for(j=0; j<nS; j++) {
            if(is_multiskill) {
                combined = getA(t,i,j);
                separate = fb->A[i][j];
                mult = combined * (1-combined)  / safe0num( separate * (1-separate) );
                fb->gradA[i][j] -= mult * this->beta[ix][j] * ((o<0)?1:getB(t,j,o)) * this->alpha[ix-1][i] / safe0num(p_O_param) + L2penalty(this->p,separate); // PENALTY
            } else {
                fb->gradA[i][j] -= this->beta[ix][j] * ((o<0)?1:fb->B[j][o]) * this->alpha[ix-1][i] / safe0num(p_O_param) + L2penalty(this->p,fb->A[i][j]); // PENALTY
            }
        }
}

void HMMProblem::setGradB (struct data* dt, FitBit *fb){
    NDAT t, ix;
    NPAR o, i;
    NPAR nS = this->p->nS;
    NUMBER combined, separate, mult;
    for(ix=dt->ix1st; ix<(dt->ix1st + dt->n); ix++) {
        t = this->p->dat_tix[ix];
        o = this->p->dat_obs->get(t);
        if(o<0) // if no observation -- skip
            continue;
        for(i=0; i<nS /*&& fitparam[2]>0*/; i++) {
            combined = getB(t,i,o);
            separate = fb->B[i][o];
            mult = (this->p->multiskill==0)?1:combined * (1-combined)  / safe0num( separate * (1-separate) );
            fb->gradB[i][o] -= mult * this->alpha[ix][i] * this->beta[ix][i] / safe0num(dt->p_O_param * getB(t,i,o)) + L2penalty(this->p,separate); // PENALTY
        }
    }
}

void HMMProblem::setGradB (NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb){
    NPAR i;
    NPAR nS = this->p->nS;
    NUMBER combined, separate, mult;
    //o < 0 taken care of outside of this
    for(i=0; i<nS /*&& fitparam[2]>0*/; i++) {
        if(is_multiskill) {
            combined = getB(t,i,o);
            separate = fb->B[i][o];
            mult = combined * (1-combined)  / safe0num( separate * (1-separate) );
            fb->gradB[i][o] -= mult * this->alpha[ix][i] * this->beta[ix][i] / safe0num(p_O_param * getB(t,i,o)) + L2penalty(this->p,fb->B[i][o]); // PENALTY
        } else {
            fb->gradB[i][o] -= this->alpha[ix][i] * this->beta[ix][i] / safe0num(p_O_param * fb->B[i][o]) + L2penalty(this->p,separate); // PENALTY
        }
    }
}

//void HMMProblem::computeGradients(NCAT ix1, NCAT n, struct data* data, FitBit *fb, NCAT target_k){//,  NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB)
//    fb->toZero(FBS_GRAD);
////    clock_t tm = clock();
////    for(int i=0;i<100; i++)
//    computeAlphaEtAl(ix1, n, data, target_k, false /*do not zero counts*/); // TODO: maybe we don't need to recompute it sometimes, e.g. in the beginning when globals were fit
////    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
////    tm = clock();
////    for(int i=0;i<100; i++)
//    computeBeta(ix1, n, data, target_k);
////    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
//    struct data *seq;
//	for(NCAT q=0; q<n; q++) {
//        seq = &data[ix1+q];
//		if( seq->tag!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
//        if(fb->PI != NULL) setGradPI(seq, fb);
//        if(fb->A  != NULL) setGradA (seq, fb);
//        if(fb->B  != NULL) setGradB (seq, fb);
//    }// for all sequences
//} // computeGradients()

void HMMProblem::computeGradients(FitBit **fbs) {
    NPAR o, n;
    NCAT *ar, k;
    data **dats;
    NDAT ix; // index in big array of variables
    bool is_multiskill = p->multiskill != 0;
    for(k=0; k<this->p->nK; k++)
        fbs[k]->toZero(FBS_GRAD); // set gradient to zero
    
//    clock_t tm;
//    tm = clock();
//    for(int i=0;i<1; i++)
    computeAlphaEtAl();// global swipe (with necessary skips)
//    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(1*CLOCKS_PER_SEC));

//    tm = clock();
//    for(int i=0;i<1; i++)
    computeBeta();// global swipe (with necessary skips)
//    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(1*CLOCKS_PER_SEC));
//    tm = clock();
    for(NDAT t=0; t<this->p->N; t++) { // for all rows
        o = this->p->dat_obs->get(t);
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]<0 ) // empty skill
            continue;
        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {  // for all skills
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            ix = dats[l]->ix1st + dats[l]->cnt;
            // do the gradients
            if(dats[l]->cnt==0) // PI only for first row
                setGradPI(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            if(dats[l]->cnt>0) // A only starting with 2nd row
                setGradA(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            if(o>=0) // B only for non-skipped observation
                setGradB(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            //
            // advance counters
            //
            dats[l]->cnt++;
        } // for all skills
    } // for all rows
    // zero counters
//    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(1*CLOCKS_PER_SEC));
    zeroCounts(this->p);
}

void HMMProblem::computeGradients(FitBit **fbs, NCAT target_k) {
    NDAT t;
    NPAR o, n;
    NCAT *ar, k;
    data **dats;
    NDAT ix; // index in big array of variables
    bool is_multiskill = p->multiskill != 0;
    fbs[target_k]->toZero(FBS_GRAD); // set gradient to zero
    computeAlphaEtAl(target_k, false /*do not zero counts*/);// global swipe (with necessary skips)
    computeBeta(target_k); // global swipe (with necessary skips)
    
    for(NDAT T=0; T<this->p->k_N[target_k]; T++) { // for all skill rows
        t = this->p->k_t[target_k][T];
        o = this->p->dat_obs->get(t);
        //
        // multi or single skill
        //
        if( !is_multiskill ) {
            k = p->dat_skill->get(t);
            ar = &k;
            n = 1;
            dats = &p->dat_skill_seq[t];
        } else {
            ar = &p->dat_multiskill->get(t)[1];
            n = p->dat_multiskill->get(t)[0];
            dats = p->dat_multiskill_seq[t];
        }
        if( ar[0]<0 ) // empty skill
            continue;
        //
        // for all skills in question
        //
        for(NPAR l=0; l<n; l++) {  // for all skills
            // if we're still fitting this one
            if(dats[l]->tag != 0)
                continue;
            k = ar[l];
            if(k!=target_k) // if it is not the target skill
                continue;
            ix = dats[l]->ix1st + dats[l]->cnt;
            // do the gradients
            if(dats[l]->cnt==0) // PI only for first row
                setGradPI(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            if(dats[l]->cnt>0) // A only starting with 2nd row
                setGradA(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            if(o>=0) // B only for non-skipped observation
                setGradB(o, t, ix, dats[l]->p_O_param, is_multiskill, fbs[k]);
            //
            // advance counters
            //
            dats[l]->cnt++;
        } // for all skills
    } // for all skill rows
    // zero counts for this skill
    for(NCAT q=0; q<this->p->k_nG[target_k]; q++)
        this->p->all_seq[ q + this->p->k_ix1stSeq[target_k] ].cnt = 0;
}


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
		fprintf(fid," %12.10f%s",this->null_obs_ratio[m],(m==(this->p->nO-1))?"\n":"\t");
    
	NCAT k;
	std::map<NCAT,std::string>::iterator it;
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		NPAR i,j,m;
		fprintf(fid,"PI\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->PI[k][i],(i==(this->p->nS-1))?"\n":"\t");
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
		fprintf(fid," %12.10f%s",this->null_obs_ratio[m],(m==(this->p->nO-1))?"\n":"\t");
	NCAT g;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<this->p->nG;g++) {
		it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		NPAR i,j,m;
		fprintf(fid,"PI\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->PI[g][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%12.10f%s",this->A[g][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%12.10f%s",this->B[g][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

//void HMMProblem::producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt) {
//    NPAR m, i;
//    NCAT k;
//    NUMBER *local_pred_inner = init1D<NUMBER>(this->p->nO);
//    for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
//    for(int l=0; l<nks; l++) {
//        for(m=0; m<this->p->nO; m++) local_pred_inner[m] = 0.0;
//        k = ks[l];
//        dt->id = k;
//        for(m=0; m<this->p->nO; m++)
//            for(i=0; i<this->p->nS; i++)
//                local_pred_inner[m] += group_skill_map[dt->g][k][i] * getB(dt,i,m);//B[i][m];
//        for(m=0; m<this->p->nO; m++)
//            local_pred[m] += local_pred_inner[m]; // local_pred[m] = 0.0;
//    }
//    if(nks>1) {
//        for(m=0; m<this->p->nO; m++)
//            local_pred[m] /= nks;
////            projectsimplex(local_pred, this->p->nO);
//    }
//    free(local_pred_inner);
//}
//
//void HMMProblem::predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill, bool only_unlabeled) {
//	NDAT t;
//	NCAT g, k;
//	NPAR i, j, m, o, isTarget;
//    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
//	NUMBER *local_pred = init1D<NUMBER>(nO); // local prediction
//	char local_know[1024];
//	NUMBER pLe[nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3D<NUMBER>(nG, nK, nS);
//    NUMBER ll = 0.0, rmse = 0.0, rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
//    NUMBER p;
//    FILE *fid; // file for storing prediction should that be necessary
//    bool output_this; // flag for turning on/off the writing out
//    if(this->p->predictions>0) {
//        fid = fopen(filename,"w");
//        if(fid == NULL)
//        {
//            fprintf(stderr,"Can't write output model file %s\n",filename);
//            exit(1);
//        }
//    }
//	// initialize
//    struct data* dt = new data;
//	
//	for(g=0; g<nG; g++)
//		for(k=0; k<nK; k++) {
//            dt->k = k;
//            dt->g = g;
//			for(i=0; i<nO; i++)
//                group_skill_map[g][k][i] =  getPI(dt,i);//PI[i];
//		}
//	
//	for(t=0; t<this->p->N; t++) {
//        output_this = true;
//		o = dat_obs->get(t);//[t];
//        if( only_unlabeled && o>-1 ) // if we only output predictions for unlabelled, it's labelled - turn off
//            output_this = false;
//		g = dat_group->get(t);//[t];
//        dt->g = g;
//        isTarget = this->p->metrics_target_obs == o;
//        NCAT *ar;
//        int n;
//        if(this->p->multiskill==0) {
//            k = dat_skill->get(t);
//            ar = &k;
//            n = 1;
//        } else {
//            ar = &dat_multiskill->get(t)[1];
//            n = dat_multiskill->get(t)[0];
//        }
//        // deal with null skill
//        if(ar[0]<0) { // if no skill label
//            isTarget = this->null_skill_obs==o;
//            rmse += pow(isTarget - this->null_skill_obs_prob,2);
//            accuracy += isTarget == (this->null_skill_obs_prob>=0.5);
//            ll -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1 - this->null_skill_obs_prob);
//            if(this->p->predictions>0 && output_this) // write predictions file if it was opened
//                for(m=0; m<nO; m++)
//                    fprintf(fid,"%12.10f%s",this->null_obs_ratio[m],(m<(nO-1))?"\t":"\n");
//            continue;
//        }
//        // produce prediction and copy to result
//        producePCorrect(group_skill_map, local_pred, ar, n, dt);
//        // update pL
//        for(int l=0; l<n; l++) {
//            //for(m=0; m<nO; m++) local_pred_inner[m] = 0.0;
//            k = ar[l];
//            dt->k = k;
//            
//            if(o>-1) { // known observations
//            // update p(L)
//                pLe_denom = 0.0;
//                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//                for(i=0; i<nS; i++) pLe_denom += group_skill_map[g][k][i] * getB(dt,i,o);
//                for(i=0; i<nS; i++) pLe[i] = group_skill_map[g][k][i] * getB(dt,i,o) / safe0num(pLe_denom);
//                // 2. L = (pLe'*A)';
//                for(i=0; i<nS; i++) group_skill_map[g][k][i] = 0.0;
//                for(j=0; j<nS; j++)
//                    for(i=0; i<nS; i++)
//                        group_skill_map[g][k][j] += pLe[i] * getA(dt,i,j);//A[i][j];
//            } else { // unknown observation
//                // 2. L = (pL'*A)';
//                for(i=0; i<nS; i++) pLe[i] = group_skill_map[g][k][i]; // copy first;
//                for(i=0; i<nS; i++) group_skill_map[g][k][i] = 0.0; // erase old value
//                for(j=0; j<nS; j++)
//                    for(i=0; i<nS; i++)
//                        group_skill_map[g][k][j] += pLe[i] * getA(dt,i,j);
//            }// ibservations
//        }
//        local_know[0] = 0;
//        for(int l=0; l<n; l++)
//            sprintf(local_know,"%s%s%12.10f",local_know,(strlen(local_know)>0)?",":"",group_skill_map[g][ ar[l] ][0]);
//        if(this->p->predictions>0 && output_this) { // write predictions file if it was opened
//            for(m=0; m<nO; m++)
//                fprintf(fid,"%12.10f%s",local_pred[m],(m<(nO-1))?"\t":"\n");
////            fprintf(fid,"%s\n",local_know);
//            //            for(i=0; i<nS; i++)
//            //                fprintf(fid,"%12.10f%s",local_know[i],(i<(nS-1))?"\t":"\n");
//        }
//        rmse += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
//        rmse_no_null += pow(isTarget-local_pred[this->p->metrics_target_obs],2);
//        accuracy += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5);
//        accuracy_no_null += isTarget == (local_pred[this->p->metrics_target_obs]>=0.5);
//        p = safe01num(local_pred[this->p->metrics_target_obs]);
//        ll -= safelog(  p)*   isTarget  +  safelog(1-p)*(1-isTarget);
//	} // for all data
//    delete(dt);
//	free(local_pred);
//    //	free(local_pred_inner);
//    free3D<NUMBER>(group_skill_map, nG, nK);
//    rmse = sqrt(rmse / this->p->N);
//    rmse_no_null = sqrt(rmse_no_null / (this->p->N - this->p->N_null));
//    if(metrics != NULL) {
//        metrics[0] = ll;
//        metrics[1] = 2*this->n_params + 2*ll;
//        metrics[2] = this->n_params*safelog(this->p->N) + 2*ll;
//        metrics[3] = rmse;
//        metrics[4] = rmse_no_null;
//        metrics[5] = accuracy/this->p->N;
//        metrics[6] = accuracy_no_null/(this->p->N-this->p->N_null);
//    }
//    if(this->p->predictions>0) // close predictions file if it was opened
//        fclose(fid);
//}

// no need computing pL, it's been computed in computeAlpha
void HMMProblem::computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE) {
    NPAR o,m,i, isTarget;
    NPAR nS = this->p->nS, nO = this->p->nO;
    NDAT t;
    NUMBER *local_pred = init1D<NUMBER>(nO);
    NUMBER *pL = init1D<NUMBER>(nS);
    NUMBER *pLe = init1D<NUMBER>(nS);
    NUMBER prob;
    NDAT N = 0;
    for(t=0; t<this->p->N; t++) {
        o = this->p->dat_obs->get(t);
        if(o<0) // if unobserved
            continue;
        isTarget = (this->p->metrics_target_obs == o);
        for(i=0;i<nO; i++) pL[i] = safe01num(getPI(t,i)); // /*safe01num*/(a_PI[i]); // init pL
        for(m=0; m<nO; m++) local_pred[m] = 0.0; // init pCorr
        for(m=0; m<nO; m++)
            for(i=0; i<nS; i++)
                local_pred[m] += pL[i] * getB(t,i,m);//a_B[i][m];
        prob = safe01num(local_pred[this->p->metrics_target_obs]);
        loglik_rmse[0] -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
        loglik_rmse[1] += pow(isTarget - prob, 2);
        loglik_rmse[2] += pow(isTarget - prob, 2); // for RMSE without null skill
        loglik_rmse[4] += this->p->metrics_target_obs == (prob>0.5);
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

void HMMProblem::computeMetrics(NUMBER* metrics) {
    computeLogLikRMSENullSkill(metrics, true /* only SE*/);
    // despite cycling on k-skill, would work for all
    computeLogLikRMSE(metrics, true /* only SE*/);//, this->p->k_nG[k], this->p->k_g_data[k]);
    metrics[3] = metrics[1]; // move Squared Errors from position 2
    metrics[3] = sqrt(metrics[3]/this->p->N);  // convert SE to RMSE
    metrics[4] = metrics[2]; // move Squared Errors from position 2
    metrics[4] = sqrt(metrics[4]/(this->p->N - this->p->N_null));  // convert SE to RMSE
    metrics[5] = 0; // Accuracy all
    metrics[6] = 0; // Accuracy no null
    metrics[1] = 2*this->n_params + 2*metrics[0]/*loglik*/;  // AIC
    metrics[2] = this->n_params*safelog(this->p->N) + 2*metrics[0]/*loglik*/;  // BIC
}

void HMMProblem::FitNullSkill(NUMBER* loglik_rmse, bool keep_SE) {
    if(this->p->nSeqNull==0) {
        this->null_obs_ratio[0] = 1; // set first obs to 1, simplex preserved
        return; // 0 loglik
    }
    NDAT count_all_null_skill = 0;
    struct data *dat; // used as pointer
    // count occurrences
    NCAT g;
    NDAT t, ix;
    NPAR isTarget, o, m;
    for(g=0; g<this->p->nSeqNull; g++) {
        dat = &this->p->all_seq[this->p->nSeq + g]; // null skill seq's are in the end of all skill seqs
        if(dat->tag != 0)
            continue; // observe block
        for(ix=0; ix<dat->n; ix++) {
            t = this->p->dat_tix[ dat->ix1st + ix];
            o = this->p->dat_obs->get( t );
            if(((int)o)>=0) { // -1 we skip \xff in char notation
                this->null_obs_ratio[ o ]++;
                count_all_null_skill++;
            }
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
    // compute RMSE & loglik
    for(g=0; g<this->p->nSeqNull; g++) {
        dat = &this->p->all_seq[this->p->nSeq + g]; // null skill seq's are in the end of all skill seqs
        if(dat->tag != 0)
            continue; // observe block
        for(ix=0; ix<dat->n; ix++) {
            t = this->p->dat_tix[ dat->ix1st + ix];
            o = this->p->dat_obs->get( t );
            if(((int)o)>=0) { // -1 we skip \xff in char notation
                isTarget = o == this->null_skill_obs;
                loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
                loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
            } // if not null observation
        }
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/count_all_null_skill);
}

void HMMProblem::computeLogLikRMSENullSkill(NUMBER* loglik_rmse, bool keep_SE) {
    // compute loglik
    NDAT count_all_null_skill = 0;
    NPAR isTarget;
    NDAT t, ix;
    NPAR o;
    struct data *dat;
    for(NCAT g=0; g<this->p->nSeqNull; g++) {
        dat = &this->p->all_seq[this->p->nSeq + g]; // null skill seq's are in the end of all skill seqs
        if(dat->tag != 0)
            continue; // observe block
        for(ix=0; ix<dat->n; ix++) {
            t = this->p->dat_tix[ dat->ix1st + ix];
            o = this->p->dat_obs->get( t );
            if(((int)o)>=0) { // -1 we skip \xff in char notation
                isTarget = o == this->null_skill_obs;
                loglik_rmse[0] -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1-this->null_skill_obs_prob);
                loglik_rmse[1] += pow(isTarget - this->null_skill_obs_prob, 2);
                count_all_null_skill++;
            }
        }
    }
    if(!keep_SE) loglik_rmse[1] = sqrt(loglik_rmse[1]/count_all_null_skill);
}

FitResult HMMProblem::GradientDescentBit(NCAT ix1, NCAT n, struct data* data, NCAT k_target, NPAR kg_flag, FitBit *fb, bool is1SkillForAll) {
    FitResult fr;
    fr.iter = 1;
    fr.pO0  = 0.0;
    fr.pO   = 0.0;
    fr.conv = 0; // converged
//    
//    
//    while( !fr.conv && fr.iter<=this->p->maxiter ) {
////        clock_t tm = clock();
////        for(int i=0;i<100; i++)
//            computeGradients(ix1, n, data, fb, kg_flag);//a_gradPI, a_gradA, a_gradB);
////        printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
//        
//        if(fr.iter==1) {
//            computeAlphaEtAl(ix1, n, data,k_target, true /*do zero counts*/);
//            fr.pO0 = HMMProblem::getSumLogPOPara(ix1, n, data);
//        }
//        // copy parameter values
//        fb->copy(FBS_PAR, FBS_PARm1);
//        // make step
//        if( fr.iter==1 || this->p->solver!=METHOD_CGD)
//            doLinearStep(ix1, n, data, k_target, fb, (is1SkillForAll)?0:-1/*copy KC0*/); // step for linked skill 0
//        else
//            doConjugateLinearStep(ix1, n, data, k_target, fb, (is1SkillForAll)?0:-1/*copy KC0*/);
//        // if we fit one skill for all
//        if(is1SkillForAll) {
//            NCAT nY = (kg_flag==0/*by skill*/)?this->p->nK:this->p->nG;
//            for(NCAT y=0; y<nY; y++) { // copy the rest
//                if(nY==k_target) continue; /// not sure what we did here :( TODO
//                NUMBER *aPI = this->getPI(y);
//                NUMBER **aA = this->getA(y);
//                NUMBER **aB = this->getB(y);
//                cpy3Params(fb->PI, fb->A, fb->B, aPI, aA, aB, this->p->nS, this->p->nO);
//            }
//        }
//        // converge?
//        fr.conv = fb->checkConvergence();
//        // report if converged
//        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
//            computeAlphaEtAl(ix1, n, data,k_target, true /*do zero counts*/);
//            fr.pO = HMMProblem::getSumLogPOPara(ix1, n, data);
//        } else if (this->p->solver==METHOD_CGD) {
//            fb->copy(FBS_GRAD, FBS_GRADm1);
//        }
//        fr.iter ++;
//    }// single skill loop
////    RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
//    fr.iter--;
    return fr;
}

NCAT HMMProblem::GradientDescentBit(FitResult *frs, FitBit **fbs) {
    NCAT done = 0;
    for(NCAT k=0; k<this->p->nK; k++) // just the non-null skills
        frs[k].iter ++;
    
//    computeGradients(fbs); // compute all that need to be computed
    
    for(NCAT k=0; k<this->p->nK; k++) { // for all skills
        computeGradients(fbs, k);
        if(fabs(frs[k].conv)==1) // already fit, skil it
            continue;
        if(frs[k].iter==1) {
//            computeAlphaEtAl(ix1, n, data,k_target, true /*do zero counts*/); // p_O_params are already well in place
            frs[k].pO0 = HMMProblem::getSumLogPOPara(k);
        }
        // copy parameter values
        fbs[k]->copy(FBS_PAR, FBS_PARm1);
        // make step
        if( frs[k].iter==1 || this->p->solver!=METHOD_CGD)
            doLinearStep(k, fbs[k]); // step for linked skill 0
        else
            doConjugateLinearStep(k, fbs[k]);
        // converge?
        frs[k].conv = 2*fbs[k]->checkConvergence();
        if(frs[k].iter == this->p->maxiter)
            frs[k].conv = -2;
        done += frs[k].conv!=0;
        
        // report if converged
        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (frs[k].conv!=0) )) {
            computeAlphaEtAl(k, true /*do zero counts*/);
            frs[k].pO = HMMProblem::getSumLogPOPara(k);
        } else if (this->p->solver==METHOD_CGD) {
            fbs[k]->copy(FBS_GRAD, FBS_GRADm1);
        }
    } // for all skills
    
    // if we fit one skill for all
//    if(is1SkillForAll) {
//        NCAT nY = (kg_flag==0/*by skill*/)?this->p->nK:this->p->nG;
//        for(NCAT y=0; y<nY; y++) { // copy the rest
//            if(nY==k_target) continue; /// not sure what we did here :( TODO
//            NUMBER *aPI = this->getPI(y);
//            NUMBER **aA = this->getA(y);
//            NUMBER **aB = this->getB(y);
//            cpy3Params(fb->PI, fb->A, fb->B, aPI, aA, aB, this->p->nS, this->p->nO);
//        }
//    }
    return done;
}

NUMBER HMMProblem::GradientDescent0() {
	NCAT x;
    NUMBER loglik = 0;
    NCAT nX, ndat;
    NDAT ix1;
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
    
    // before starting init Alpha and Beta variable arrays and compute them
    this->initAlpha();
    this->initBeta();
    zeroCounts(this->p);
    
    // time testing
//    clock_t tm = clock();
//    for(int i=0;i<100; i++)
        this->computeAlphaEtAl();
//    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
    this->computeBeta();
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill>0) {
        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(0/*start at squence 0*/, this->p->nSeq/*for all sequences*/, this->p->all_seq, 0/*use skill 0*/, 0/* by skill*/, fb, true /*is1SkillForAll*/);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
	}
	//
	// Main fit
	//
	for(x=0; x<nX && this->p->single_skill!=2; x++) { // if not "force single skill" too
        if(this->p->structure==STRUCTURE_SKILL) {
            ix1 = this->p->k_ix1stSeq[x];
            ndat = this->p->k_nG[x];
        } else if(this->p->structure==STRUCTURE_GROUP) {
            fprintf(stderr,"Fitting by Group is not implemented\n");
//            xndat = this->p->g_numk[x];
//            x_data = this->p->g_k_data[x];
        }
        fb->linkPar( this->getPI(x), this->getA(x), this->getB(x));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(ix1, ndat, this->p->all_seq, x, (this->p->structure!=STRUCTURE_SKILL), fb, false /*is1SkillForAll*/);
        if( !this->p->quiet )
            printf("skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", x, fr.iter,fr.pO0,fr.pO,fr.conv);
	} // for all skills
    delete fb; // that takes care of *m1, and *GRAD
    return loglik;
}


NUMBER HMMProblem::GradientDescent() {
    NUMBER loglik = 0, tol = this->p->tol;
    NCAT k, q, nK = this->p->nK, nG = this->p->nG, dones = 0 /*all converged or done counts*/, done = 0 /*done at recent iteration*/;
    NPAR nS = this->p->nS, nO = this->p->nO;
    bool isCGD = this->p->solver==METHOD_CGD; // if we use Conjugate Gradient Descent
    int iter = 1; // iteration count
    
    // before starting init Alpha and Beta variable arrays and compute them
    this->initAlpha();
    this->initBeta();
    zeroCounts(this->p);
//    // time testing
////    clock_t tm = clock();
////    for(int i=0;i<100; i++)
//    this->computeAlphaEtAl(); // << NOT NECESSARY
////    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
////    clock_t tm = clock();
////    for(int i=0;i<100; i++)
//    this->computeBeta();  // << NOT NECESSARY
////    printf("time %12.10f seconds\n",(NUMBER)(clock()-tm)/(100*CLOCKS_PER_SEC));
	
//	//
//	// fit all as 1 skill first
//	//
//	if(this->p->single_skill>0) {
//        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));// link skill 0 (we'll copy fit parameters to others
//    NCAT* original_ks = Calloc(NCAT, this->p->nSeq);
//    for(x=0; x<this->p->nSeq; x++) { original_ks[x] = this->p->all_data[x].k; this->p->all_data[x].k = 0; } // save progonal k's
//    fr = GradientDescentBit(fb);
//    for(x=0; x<this->p->nSeq; x++) { this->p->all_data[x].k = original_ks[x]; } // restore original k's
//    free(original_ks);
//        if( !this->p->quiet )
//            printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
//	}

	//
	// Main fit
	//
    
    // init fit results fit bits and link them to skill parameters
    FitBit   **fbs = Calloc(FitBit*,   (size_t)this->p->nK);
    FitResult *frs = Calloc(FitResult, (size_t)this->p->nK);
    for(k=0; k<this->p->nK; k++) {
        fbs[k] = new FitBit(nS, nO, nK, nG, tol);
        fbs[k]->init(FBS_PARm1);
        fbs[k]->init(FBS_GRAD);
        if(isCGD) {
            fbs[k]->init(FBS_GRADm1);
            fbs[k]->init(FBS_DIRm1);
        }
        fbs[k]->linkPar( this->getPI(k), this->getA(k), this->getB(k));// link skill 0 (we'll copy fit parameters to others
        frs[k].iter = 0;
        frs[k].pO0  = 0.0;
        frs[k].pO   = 0.0;
        frs[k].conv = 0; // converged 2, done -2
    }
    // fit loop
    while(dones < nK) {
        done = GradientDescentBit(frs, fbs); // GDBit runs here, use conv = 2 as `freshly converged`, -2 as `freshly done`, 1/-1 `converged`/`done`
        // FINISHED HERE ^^^^
        for(k=0; k<nK && done>0; k++) { // seek out those done and account for them
            if(fabs(frs[k].conv)==2) { // freshly done
                // turn off fitting on this skill's sequences
                for(q=0; q<this->p->k_nG[k]; q++)
                    this->p->all_seq[ q+this->p->k_ix1stSeq[k] ].tag = 1; // off
                // house-keeping
                loglik += frs[k].pO; // update loglik
                frs[k].conv -= (frs[k].conv>0)?1:-1; // if negative increase, if positive decrease to -1/1
                frs[k].iter = iter;
                if( !this->p->quiet )
                    printf("skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", k, frs[k].iter, frs[k].pO0, frs[k].pO, frs[k].conv);
                done--;  // done with freshly converged
                dones++; // now count it
            } // freshly done
        } // seek out those done and account for them
        iter++;
    } // until we are all done
        
    //
    // Recycle fit results and fit bits
    //
    for(k=0; k<nK; k++)
        delete fbs[k]; // that takes care of all *m1, and *GRAD
    free(fbs);
    free(frs);
    
    return loglik;
}

NUMBER HMMProblem::BaumWelchSkill() {
    fprintf(stderr,"BaumWelch method is not implemented\n");

    
//	NCAT k;
//    NUMBER loglik = 0;
//	
//    FitBit *fb = new FitBit(this->p);
//    fb->init(FBS_PARm1);
//	
//	int conv;
//	int iter; // iteration count
//	NUMBER pO0, pO;
////    bool conv_flags[3] = {true, true, true};
//	
//	//
//	// fit all as 1 skill first
//	//
//    //	if(this->p->fit_single_skill) {
//    //		for(k=0; k<this->p->nK; k++) this->p->mask_skill[k] = true; // mask with k'th==true to be computed
//    //		iter = 1;
//    //		pO0 = 0.0;
//    //		conv = 0; // converged
//    //		while( !conv && iter<=this->p->maxiter ) {
//    //			if(iter>1) {
//    //				toZero1D<NUMBER>(gradPI, this->p->nS);
//    //				toZero2D<NUMBER>(gradA,  this->p->nS, this->p->nS);
//    //				toZero2D<NUMBER>(gradB,  this->p->nS, this->p->nO);
//    //			}
//    //			hmm->computeGradients();
//    //			if(iter==1)
//    //				pO0 = hmm->getSumLogPOParaK(flog);
//    //			// add gradients
//    //			for(k=0; k<this->p->nK; k++) {
//    //				add1DNumbersWeighted(hmm->getGradPI(k), gradPI, this->p->nS, 1.0);
//    //				add2DNumbersWeighted(hmm->getGradA(k),  gradA,  this->p->nS, this->p->nS, 1.0);
//    //				add2DNumbersWeighted(hmm->getGradB(k),  gradB,  this->p->nS, this->p->nO, 1.0);
//    //			}
//    //			// copy old SAVED! values for params, just for skill #0 is enough
//    //			cpy1D<NUMBER>(hmm->getPI(0), PI_m1, this->p->nS);
//    //			cpy2D<NUMBER>(hmm->getA(0),  A_m1,  this->p->nS, this->p->nS);
//    //			cpy2D<NUMBER>(hmm->getB(0),  B_m1,  this->p->nS, this->p->nO);
//    //
//    //			// make step
//    //			for(k=0; k<this->p->nK; k++)
//    //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
//    //			// check convergence, on any skill, e.g. #0
//    //			conv = checkConvergence(hmm->getPI(k), hmm->getA(k), hmm->getB(k), PI_m1, A_m1, B_m1, param);
//    //
//    //			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||/**/ (conv || iter==this->p->maxiter) )) {
//    //				hmm->computeAlphaAndPOParam();
//    //				NUMBER pO = hmm->getSumLogPOParaK(flog);
//    //				printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",iter,pO0,pO,conv);
//    //			}
//    //			iter ++;
//    //		}// single skill loop
//    //
//    //		free(gradPI);				   // from now on, just a pointer
//    //		free2D<NUMBER>(gradA, this->p->nS); // from now on, just a pointer
//    //		free2D<NUMBER>(gradB, this->p->nS); // from now on, just a pointer
//    //	}
//	
//	//
//	// Main fit
//	//
//    //	for(k=0; k<this->p->nK; k++) this->p->mask_skill[k] = false; // mask with k'th==true to be computed
//	for(k=0; k<this->p->nK; k++) {  // for(k=218; k<219; k++) { //
//        NCAT xndat = this->p->k_nG[k];
//        struct data** x_data =  this->p->k_g_data[k];
//        
//		conv = 0; // converged
//		iter = 1; // iteration count
//		pO0 = 0.0;
//        pO = 0.0;
//		
//        fb->linkPar(this->getPI(k), this->getA(k), this->getB(k));
//		
//		while( !conv && iter<=this->p->maxiter ) {
//			if(iter==1) {
//                computeAlphaAndPOParam(xndat, x_data);
//				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
//			}
//			
//			// copy old SAVED! values for params
//            fb->copy(FBS_PAR, FBS_PARm1);
//			
//            //			hmm->zeroLabelsK(k); // reset blocking labels // THIS IS NOT DONE HERE
//			doBaumWelchStep(xndat, this->p->k_g_data[k], fb);// PI, A, B);
//            
//			// check convergence
//            conv = fb->checkConvergence();
//			
//			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
//                computeAlphaAndPOParam(xndat, x_data);
//                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
//                loglik += pO*(pO>0);
//                if(!this->p->quiet)
//                    printf("skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",k,iter,pO0,pO,conv);
//			}
//			
//			iter ++;
//		} // main solver loop
//		// recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
//        RecycleFitData(xndat, x_data, this->p);
//		// recycle
//	} // for all skills
//    delete fb;
//    return loglik;
    return 0;
}

NUMBER HMMProblem::doLinearStep(NCAT ix1, NCAT ndat, struct data* data, NCAT k_target, FitBit *fb, NCAT copy) {//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
//	NPAR i,j,m;
//    NPAR nS = this->p->nS, nO = this->p->nO;
//	// first scale down gradients
//    fb->doLog10ScaleGentle(FBS_GRAD);
//    //	doLog10Scale1DGentle(fb->gradPI, fb->PI, nS);
//    //	doLog10Scale2DGentle(fb->gradA,  fb->A,  nS, nS);
//    //	doLog10Scale2DGentle(fb->gradB,  fb->B,  nS, nO);
//	
//    //    NUMBER *PI_cpy, ** A_cpy, ** B_cpy; // replace with PARcopy in fb
//    //    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
//    fb->init(FBS_PARcopy);
//    
//	NUMBER e = this->p->ArmijoSeed; // step seed
//	bool compliesArmijo = false;
//    //	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
//	NUMBER f_xk = HMMProblem::getSumLogPOPara(ix1, ndat, data);
//	NUMBER f_xkplus1;
//	
//    fb->copy(FBS_PAR, FBS_PARcopy);
//    //    cpy3Params(fb->PI, fb->A, fb->B, PI_cpy, A_cpy, B_cpy, nS, nO);
//    //	cpy1D<NUMBER>(fb->PI, PI_cpy, nS); // save copy
//    //	cpy2D<NUMBER>(fb->A,  A_cpy,  nS, nS); // save copy
//    //	cpy2D<NUMBER>(fb->B,  B_cpy,  nS, nO); // save copy
//	// compute p_k * -p_k
//	NUMBER p_k_by_neg_p_k = 0;
//	for(i=0; i<nS; i++)
//	{
//		if(fb->PI != NULL) p_k_by_neg_p_k -= fb->gradPI[i]*fb->gradPI[i];
//		if(fb->A  != NULL) for(j=0; j<nS; j++) p_k_by_neg_p_k -= fb->gradA[i][j]*fb->gradA[i][j];
//		if(fb->B  != NULL) for(m=0; m<nO; m++) p_k_by_neg_p_k -= fb->gradB[i][m]*fb->gradB[i][m];
//	}
//	int iter = 0; // limit iter steps to 20, via ArmijoMinStep (now 10)
//	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
//		// update
//		for(i=0; i<nS; i++) {
//			if(fb->PI != NULL) fb->PI[i] = fb->PIcopy[i] - e * fb->gradPI[i];
//            if(fb->A  != NULL)
//                for(j=0; j<nS; j++)
//                    fb->A[i][j] = fb->Acopy[i][j] - e * fb->gradA[i][j];
//            if(fb->B  != NULL)
//                for(m=0; m<nO; m++)
//                    fb->B[i][m] = fb->Bcopy[i][m] - e * fb->gradB[i][m];
//		}
//		// scale
//		if( !this->hasNon01Constraints() ) {
//			if(fb->PI != NULL) projectsimplex(fb->PI, nS);
//			for(i=0; i<nS; i++) {
//				if(fb->A  != NULL) projectsimplex(fb->A[i], nS);
//				if(fb->B  != NULL) projectsimplex(fb->B[i], nS);
//			}
//		} else {
//			if(fb->PI != NULL) projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
//			for(i=0; i<nS; i++) {
//				if(fb->A  != NULL) projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
//				if(fb->B  != NULL) projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
//			}
//		}
//        if(copy >= 0) { // copy parameters from position 'copy' to others
//            for(NCAT x=0; (fb->PI != NULL) && x<sizes[0]; x++)
//                if(x!=copy)
//                    cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
//            for(NCAT x=0; (fb->A  != NULL) && x<sizes[1]; x++)
//                if(x!=copy)
//                    cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
//            for(NCAT x=0; (fb->B  != NULL) && x<sizes[2]; x++)
//                if(x!=copy)
//                    cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
//        }
//		// recompute alpha and p(O|param)
//		computeAlphaEtAl(ix1, ndat, data, k_target, true /*do zero counts*/);
//		// compute f(x_{k+1})
//		f_xkplus1 = HMMProblem::getSumLogPOPara(ix1, ndat, data);
//		// compute Armijo compliance
//		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
//		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
//		iter++;
//	} // armijo loop
//    if(!compliesArmijo) { // we couldn't step away from current, copy the inital point back
//        e = 0;
//        fb->copy(FBS_PARcopy, FBS_PAR);
//        f_xkplus1 = f_xk;// do not recompute, use old; HMMProblem::getSumLogPOPara(ix1, ndat, data); //// DELETE
//    }
//    //    RecycleFitData(xndat, x_data, this->p);
//    fb->destroy(FBS_PARcopy);
//    return e;
    return 0;
} // doLinearStep

NUMBER HMMProblem::doLinearStep(NCAT k, FitBit *fb) {//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO;
	// first scale down gradients
    fb->doLog10ScaleGentle(FBS_GRAD);
    fb->init(FBS_PARcopy);

	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
    //	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xk = HMMProblem::getSumLogPOPara(k);
	NUMBER f_xkplus1;
	
    fb->copy(FBS_PAR, FBS_PARcopy);
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
//        if(copy >= 0) { // copy parameters from position 'copy' to others
//            for(NCAT x=0; (fb->PI != NULL) && x<sizes[0]; x++)
//                if(x!=copy)
//                    cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
//            for(NCAT x=0; (fb->A  != NULL) && x<sizes[1]; x++)
//                if(x!=copy)
//                    cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
//            for(NCAT x=0; (fb->B  != NULL) && x<sizes[2]; x++)
//                if(x!=copy)
//                    cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
//        }
		// recompute alpha and p(O|param)
		computeAlphaEtAl(k, true /*do zero counts*/);
		// compute f(x_{k+1})
		f_xkplus1 = HMMProblem::getSumLogPOPara(k);
		// compute Armijo compliance
		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
		iter++;
	} // armijo loop
    if(!compliesArmijo) { // we couldn't step away from current, copy the inital point back
        e = 0;
        fb->copy(FBS_PARcopy, FBS_PAR);
        f_xkplus1 = f_xk;// do not recompute, use old; HMMProblem::getSumLogPOPara(ix1, ndat, data); //// DELETE
    }
    //    RecycleFitData(xndat, x_data, this->p);
    fb->destroy(FBS_PARcopy);
    return e;
} // doLinearStep

NUMBER HMMProblem::doConjugateLinearStep(NCAT ix1, NCAT n, struct data* data, NCAT k_target, FitBit *fb, NCAT copy) {
//	NPAR i,j,m;
//    NPAR nS = this->p->nS, nO = this->p->nO;
//	// first scale down gradients
//    fb->doLog10ScaleGentle(FBS_GRAD);
//    
//    // compute beta_gradient_direction
//    NUMBER beta_grad = 0, beta_grad_num = 0, beta_grad_den = 0;
//    
//    switch (this->p->solver_setting) {
//        case 1: // Fletcher-Reeves
//            for(i=0; i<nS; i++)
//            {
//                if(fb->PI != NULL) {
//                    beta_grad_num = fb->gradPI  [i]*fb->gradPI   [i];
//                    beta_grad_den = fb->gradPIm1[i]*fb->gradPIm1[i];
//                }
//                if(fb->A  != NULL)
//                    for(j=0; j<nS; j++) {
//                        beta_grad_num = fb->gradA  [i][j]*fb->gradA   [i][j];
//                        beta_grad_den = fb->gradAm1[i][j]*fb->gradAm1[i][j];
//                    }
//                if(fb->B  != NULL)
//                    for(m=0; m<nO; m++) {
//                        beta_grad_num = fb->gradB  [i][m]*fb->gradB  [i][m];
//                        beta_grad_den = fb->gradBm1[i][m]*fb->gradBm1[i][m];
//                    }
//            }
//            break;
//        case 2: // Polak–Ribiere
//            for(i=0; i<nS; i++)
//            {
//                if(fb->PI != NULL) {
//                    beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
//                    beta_grad_den =  fb->gradPIm1[i]*fb->gradPIm1[i];
//                }
//                if(fb->A != NULL)
//                    for(j=0; j<nS; j++) {
//                        beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
//                        beta_grad_den =  fb->gradAm1[i][j]*fb->gradAm1[i][j];
//                    }
//                if(fb->B  != NULL)
//                    for(m=0; m<nO; m++) {
//                        beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
//                        beta_grad_den =  fb->gradBm1[i][m]*fb->gradBm1[i][m];
//                    }
//            }
//            break;
//        case 3: // Hestenes-Stiefel
//            for(i=0; i<nS; i++)
//            {
//                if(fb->PI != NULL) {
//                    beta_grad_num = -fb->gradPI[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
//                    beta_grad_den =  fb->dirPIm1[i]*(-fb->gradPI[i] + fb->gradPIm1[i]);
//                }
//                if(fb->A  != NULL)
//                    for(j=0; j<nS; j++) {
//                        beta_grad_num = -fb->gradA[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
//                        beta_grad_den =  fb->dirAm1[i][j]*(-fb->gradA[i][j] + fb->gradAm1[i][j]);
//                    }
//                if(fb->B  != NULL)
//                    for(m=0; m<nO; m++) {
//                        beta_grad_num = -fb->gradB[i][j]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
//                        beta_grad_den =  fb->dirBm1[i][m]*(-fb->gradB[i][j] + fb->gradBm1[i][j]);
//                    }
//            }
//            break;
//        default:
//            fprintf(stderr,"Wrong stepping algorithm (%d) specified for Conjugate Gradient Descent\n",this->p->solver_setting);
//            break;
//    }
//    beta_grad = beta_grad_num / safe0num(beta_grad_den);
//    beta_grad = (beta_grad>=0)?beta_grad:0;
//	
//    // compute new direction (in place of old)
//    fb->toZero(FBS_DIRm1);
//	for(i=0; i<nS; i++)
//	{
//		if(fb->PI != NULL) fb->dirPIm1[i] = -fb->gradPI[i] + beta_grad * fb->dirPIm1[i];
//		if(fb->A  != NULL) for(j=0; j<nS; j++) fb->dirAm1[i][j] = -fb->gradA[i][j] + beta_grad * fb->dirAm1[i][j];
//		if(fb->B  != NULL) for(m=0; m<nO; m++) fb->dirBm1[i][m] = -fb->gradB[i][m] + beta_grad * fb->dirBm1[i][m];
//	}
//	// scale down direction
//    fb->doLog10ScaleGentle(FBS_DIRm1);
//    
//    fb->init(FBS_PARcopy);
//    
//	NUMBER e = this->p->ArmijoSeed; // step seed
//	bool compliesArmijo = false;
//	NUMBER f_xk = HMMProblem::getSumLogPOPara(ix1, n, data);
//	NUMBER f_xkplus1;
//	
//    fb->copy(FBS_PAR, FBS_PARcopy);
//	// compute p_k * -p_k >>>> now current gradient by current direction
//	NUMBER p_k_by_neg_p_k = 0;
//	for(i=0; i<nS; i++)
//	{
//		if(fb->PI != NULL) p_k_by_neg_p_k = fb->gradPI[i]*fb->dirPIm1[i];
//		if(fb->A  != NULL) for(j=0; j<nS; j++) p_k_by_neg_p_k = fb->gradA[i][j]*fb->dirAm1[i][j];
//		if(fb->B  != NULL) for(m=0; m<nO; m++) p_k_by_neg_p_k = fb->gradB[i][m]*fb->dirBm1[i][m];
//	}
//	int iter = 0; // limit iter steps to 20, actually now 10, via ArmijoMinStep
//	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
//		// update
//		for(i=0; i<nS; i++) {
//			if(fb->PI != NULL) fb->PI[i] = fb->PIcopy[i] + e * fb->dirPIm1[i];
//            if(fb->A  != NULL)
//                for(j=0; j<nS; j++)
//                    fb->A[i][j] = fb->Acopy[i][j] + e * fb->dirAm1[i][j];
//            if(fb->B  != NULL)
//                for(m=0; m<nO; m++)
//                    fb->B[i][m] = fb->Bcopy[i][m] + e * fb->dirBm1[i][m];
//		}
//		// scale
//		if( !this->hasNon01Constraints() ) {
//			if(fb->PI != NULL) projectsimplex(fb->PI, nS);
//			for(i=0; i<nS; i++) {
//				if(fb->A != NULL) projectsimplex(fb->A[i], nS);
//				if(fb->B != NULL) projectsimplex(fb->B[i], nS);
//			}
//		} else {
//			if(fb->PI != NULL) projectsimplexbounded(fb->PI, this->getLbPI(), this->getUbPI(), nS);
//			for(i=0; i<nS; i++) {
//				if(fb->A != NULL) projectsimplexbounded(fb->A[i], this->getLbA()[i], this->getUbA()[i], nS);
//				if(fb->B != NULL) projectsimplexbounded(fb->B[i], this->getLbB()[i], this->getUbB()[i], nS);
//			}
//		}
//        // copy parameters from position 'copy' to others
//        if(copy >= 0) {
//            if(fb->PI != NULL)
//                for(NCAT x=0; x<sizes[0]; x++)
//                    if(x!=copy)
//                        cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
//            if(fb->A  != NULL)
//                for(NCAT x=0; x<sizes[1]; x++)
//                    if(x!=copy)
//                        cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
//            if(fb->B  != NULL)
//                for(NCAT x=0; x<sizes[2]; x++)
//                    if(x!=copy)
//                        cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
//        }
//		// recompute alpha and p(O|param)
//		computeAlphaEtAl(ix1, n, data, k_target, true /*do zero counts*/);
//		// compute f(x_{k+1})
//		f_xkplus1 = HMMProblem::getSumLogPOPara(ix1, n, data);
//		// compute Armijo compliance
//		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
//		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
//		iter++;
//	} // armijo loop
//    if(!compliesArmijo) { // failed to step away from initial, reinstate the inital parameters
//        e = 0;
//        fb->copy(FBS_PARcopy, FBS_PAR);
//    }
////    RecycleFitData(xndat, x_data, this->p);
//    fb->destroy(FBS_PARcopy);
//    return e;
    return 0;
} // doLinearStep

NUMBER HMMProblem::doConjugateLinearStep(NCAT k, FitBit *fb) {
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
	NUMBER f_xk = HMMProblem::getSumLogPOPara(k);
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
//        // copy parameters from position 'copy' to others
//        if(copy >= 0) {
//            if(fb->PI != NULL)
//                for(NCAT x=0; x<sizes[0]; x++)
//                    if(x!=copy)
//                        cpy1D<NUMBER>(fb->PI, this->PI[x], nS);
//            if(fb->A  != NULL)
//                for(NCAT x=0; x<sizes[1]; x++)
//                    if(x!=copy)
//                        cpy2D<NUMBER>(fb->A, this->A[x], nS, nS);
//            if(fb->B  != NULL)
//                for(NCAT x=0; x<sizes[2]; x++)
//                    if(x!=copy)
//                        cpy2D<NUMBER>(fb->B, this->B[x], nS, nO);
//        }
		// recompute alpha and p(O|param)
		computeAlphaEtAl(k, true /*do zero counts*/);
		// compute f(x_{k+1})
		f_xkplus1 = HMMProblem::getSumLogPOPara(k);
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

NUMBER HMMProblem::doBarzilaiBorweinStep(NCAT ndat, struct data* data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1) {
    
    fprintf(stderr,"Barzilai-Borwein step is not implemented\n");
    /*
	NPAR i,j,m;
    NPAR nS = this->p->nS, nO = this->p->nO;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
	doLog10Scale2DGentle(a_gradB,  a_B,  nS, nO);
    
    // compute s_k_m1
  	NUMBER *s_k_m1_PI = init1D<NUMBER>(nS);
	NUMBER **s_k_m1_A = init2D<NUMBER>(nS,nS);
	NUMBER **s_k_m1_B = init2D<NUMBER>(nS,nS);
	for(i=0; i<nS; i++)
	{
		s_k_m1_PI[i] = a_PI[i] - a_PI_m1[i];
		for(j=0; j<nS; j++) s_k_m1_A[i][j] = a_A[i][j] - a_A_m1[i][j];
		for(m=0; m<this->   p->nO; m++) s_k_m1_B[i][m] = a_B[i][m] - a_B_m1[i][m];
	}
    // compute alpha_step
    NUMBER alpha_step = 0, alpha_step_num = 0, alpha_step_den = 0;
    // Barzilai Borweig: s' * s / ( s' * (g-g_m1) )
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
    */
    return 0;
}
void HMMProblem::doBaumWelchStep(NCAT ndat, struct data* data, FitBit *fb) {//, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B) {
    fprintf(stderr,"Baum-Welch step is not implemented\n");
    /*
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
	
	NUMBER * b_PI = init1D<NUMBER>(this->p->nS);
	NUMBER ** b_A_num = init2D<NUMBER>(this->p->nS, this->p->nS);
	NUMBER ** b_A_den = init2D<NUMBER>(this->p->nS, this->p->nS);
	NUMBER ** b_B_num = init2D<NUMBER>(this->p->nS, this->p->nO);
	NUMBER ** b_B_den = init2D<NUMBER>(this->p->nS, this->p->nO);
	// compute sums PI
	NUMBER sum_p_O_param = 0;
	for(x=0; x<xndat; x++) {
        if( x_data[x]->tag!=0 ) continue;
		sum_p_O_param += x_data[x]->p_O_param;
    }
	for(x=0; x<xndat; x++) {
        if( x_data[x]->tag!=0 ) continue;
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
    */
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
//	NUMBER *local_pred = init1D<NUMBER>(this->p->nO); // local prediction
//	NUMBER pLe[this->p->nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3D<NUMBER>(this->p->nG, this->p->nK, this->p->nS);
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
//    free3D<NUMBER>(group_skill_map, this->p->nG, this->p->nK);
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
//	NUMBER *local_pred = init1D<NUMBER>(this->p->nO); // local prediction
//	NUMBER pLe[this->p->nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3D<NUMBER>(this->p->nG, this->p->nK, this->p->nS);
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
//    free3D<NUMBER>(group_skill_map, this->p->nG, this->p->nK);
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

void HMMProblem::readNullObsRatio(FILE *fid, NDAT *line_no) {
	NPAR i;
	//
	// read null skill ratios
	//
    fscanf(fid, "Null skill ratios\t");
    this->null_obs_ratio =Calloc(NUMBER, (size_t)this->p->nO);
    this->null_skill_obs      = 0;
    this->null_skill_obs_prob = 0;
	for(i=0; i<this->p->nO; i++) {
        if( i==(this->p->nO-1) ) // end
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


void HMMProblem::readModelBody(FILE *fid, NDAT *line_no) {
	NPAR i,j,m;
	NCAT k = 0;
	string s;
    char col[2048];
    //
    readNullObsRatio(fid, line_no);
    //
    // init param
    //
    this->p->map_group_fwd = new map<string,NCAT>();
    this->p->map_group_bwd = new map<NCAT,string>();
    this->p->map_skill_fwd = new map<string,NCAT>();
    this->p->map_skill_bwd = new map<NCAT,string>();
	//
	// read skills
	//
	for(k=0; k<this->p->nK; k++) {
		// read skill label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;
		this->p->map_skill_fwd->insert(pair<string,NCAT>(s, this->p->map_skill_fwd->size()));
		this->p->map_skill_bwd->insert(pair<NCAT,string>(this->p->map_skill_bwd->size(), s));

        // read PI
        fscanf(fid,"PI\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->PI[k][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->PI[k][i] = atof(col);
        (*line_no)++;
		// read A
        fscanf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++) {
                if(i==(this->p->nS-1) && j==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->A[k][i][j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->A[k][i][j] = atof(col); 
                }
			}
        (*line_no)++;
		// read B
        fscanf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nS; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->B[k][i][m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->B[k][i][m] = atof(col);
                }
			}
        (*line_no)++;
	} // for all k
}

void HMMProblem::cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO) {
    cpy1D<NUMBER>(soursePI, targetPI, nS);
    cpy2D<NUMBER>(sourseA,  targetA,  nS, nS);
    cpy2D<NUMBER>(sourseB,  targetB,  nS, nO);
}

