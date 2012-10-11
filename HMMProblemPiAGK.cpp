//
//  HMMProblemPiAGK.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 9/6/12.
//
//

#include "HMMProblemPiAGK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemPiAGK::HMMProblemPiAGK(struct param *param) {
    init(param);
}

void HMMProblemPiAGK::init(struct param *param) {
    this->n_params = param->nK * 4 + param->nG * 2;
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    
    NUMBER *a_PI = init1DNumber(this->p->nS);// PI[0] = 0.5; PI[1] = 0.5;
	NUMBER **a_A = init2DNumber(this->p->nS,this->p->nS);// A[0][0] = 1.0; A[0][1] = 0.0; A[1][0] = 0.4; A[1][1] = 0.6;
	NUMBER **a_B = init2DNumber(this->p->nS,this->p->nO);// B[0][0] = 0.8; B[0][1] = 0.2; B[1][0] = 0.2; B[1][1] = 0.8;
    
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
	for(i=0; i<((this->p->nS)-1); i++) {
		a_PI[i] = this->p->init_params[i];
		sumPI  += this->p->init_params[i];
	}
	a_PI[this->p->nS-1] = 1 - sumPI;
	// populate A
	offset = this->p->nS-1;
	for(i=0; i<this->p->nS; i++) {
		for(j=0; j<((this->p->nS)-1); j++) {
			idx = offset + i*((this->p->nS)-1) + j;
			a_A[i][j] = this->p->init_params[idx];
			sumA[i]  += this->p->init_params[idx];
		}
		a_A[i][((this->p->nS)-1)]  = 1 - sumA[i];
	}
	// polupale B
	offset = (this->p->nS-1) + this->p->nS*(this->p->nS-1);
	for(i=0; i<this->p->nS; i++) {
		for(j=0; j<((this->p->nO)-1); j++) {
			idx = offset + i*((this->p->nO)-1) + j;
			a_B[i][j] = this->p->init_params[idx];
			sumB[i] += this->p->init_params[idx];
		}
		a_B[i][((this->p->nO)-1)]  = 1 - sumB[i];
	}
    
    // mass produce PI's/PIg's, A's, B's
	if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
		this->PI  = init2DNumber(this->p->nK, this->p->nS);
		this->PIg = init2DNumber(this->p->nG, this->p->nS);
		this->A   = init3DNumber(this->p->nK, this->p->nS, this->p->nS);
		this->Ag  = init3DNumber(this->p->nG, this->p->nS, this->p->nS);
		this->B   = init3DNumber(this->p->nK, this->p->nS, this->p->nO);
        NCAT x;
		for(x=0; x<this->p->nK; x++) {
			cpy1DNumber(a_PI, this->PI[x], this->p->nS);
			cpy2DNumber(a_A,  this->A[x],  this->p->nS, this->p->nS);
			cpy2DNumber(a_B,  this->B[x],  this->p->nS, this->p->nO);
        }
        // PIg start with "no-effect" params,PI[i] = 1/nS
		for(x=0; x<this->p->nG; x++) {
            for(i=0; i<this->p->nS; i++) {
                this->PIg[x][i] = (NUMBER)1/this->p->nS;
                for(j=0; j<this->p->nO; j++)
                    this->Ag[x][i][j] = (NUMBER)1/this->p->nS;
            }
        }
	} else {
		fprintf(stderr,"params do not meet constraints.\n");
		exit(1);
	}
    // destroy setup params
	free(a_PI);
	free2DNumber(a_A, this->p->nS);
	free2DNumber(a_B, this->p->nS);
	
    // populate boundaries
	// populate lb*/ub*
	// *PI
    this->lbPI = init1DNumber(this->p->nS);
    this->ubPI = init1DNumber(this->p->nS);
    this->lbA  = init2DNumber(this->p->nS, this->p->nS);
    this->ubA  = init2DNumber(this->p->nS, this->p->nS);
    this->lbB  = init2DNumber(this->p->nS, this->p->nO);
    this->ubB  = init2DNumber(this->p->nS, this->p->nO);
	for(i=0; i<this->p->nS; i++) {
		lbPI[i] = this->p->param_lo[i];
		ubPI[i] = this->p->param_hi[i];
	}
	// *A
	offset = this->p->nS;
	for(i=0; i<this->p->nS; i++)
		for(j=0; j<this->p->nS; j++) {
			idx = offset + i*this->p->nS + j;
			lbA[i][j] = this->p->param_lo[idx];
			ubA[i][j] = this->p->param_hi[idx];
		}
	// *B
	offset = this->p->nS + this->p->nS*this->p->nS;
	for(i=0; i<this->p->nS; i++)
		for(j=0; j<this->p->nO; j++) {
			idx = offset + i*this->p->nS + j;
			lbB[i][j] = this->p->param_lo[idx];
			ubB[i][j] = this->p->param_hi[idx];
		}
    
	this->gradPIk = NULL;
	this->gradPIg = NULL;
	this->gradA = NULL;
	this->gradB = NULL;
    this->fitK = NULL;
    this->fitG = NULL;
    fitK_countG = NULL;
}

void HMMProblemPiAGK::destroy() {
	// destroy model data
    free(this->null_obs_ratio);
	free2DNumber(this->PIg, this->p->nG);
	free2DNumber(this->PI, this->p->nK);
	free3DNumber(this->A,  this->p->nK, this->p->nS);
	free3DNumber(this->Ag, this->p->nG, this->p->nS);
	free3DNumber(this->B,  this->p->nK, this->p->nS);
	if(this->lbPI!=NULL) free(this->lbPI);
	if(this->ubPI!=NULL) free(this->ubPI);
	if(this->lbA!=NULL) free2DNumber(this->lbA, this->p->nS);
	if(this->ubA!=NULL) free2DNumber(this->ubA, this->p->nS);
	if(this->lbB!=NULL) free2DNumber(this->lbB, this->p->nS);
	if(this->ubB!=NULL) free2DNumber(this->ubB, this->p->nS);
	// destroy fitting data
	// gradPIk,g
	if ( this->gradPIk != NULL) {
		free2DNumber(this->gradPIk, this->p->nK);
		this->gradPIk = NULL;
	}
	if ( this->gradPIg != NULL) {
		free2DNumber(this->gradPIg, this->p->nG);
		this->gradPIg = NULL;
	}
	// gradA
	if ( this->gradAk != NULL) {
		free3DNumber(this->gradAk, this->p->nK, this->p->nS);
		this->gradAk = NULL;
	}
	if ( this->gradAg != NULL) {
		free3DNumber(this->gradAg, this->p->nG, this->p->nS);
		this->gradAg = NULL;
	}
	// gradB
	if ( this->gradB != NULL) {
		free3DNumber(this->gradB, this->p->nK, this->p->nS);
		this->gradB = NULL;
	}
    // free fit flags
    if(this->fitK != NULL)        free(this->fitK);
    if(this->fitG != NULL)        free(this->fitG);
    if(this->fitK_countG != NULL) free(this->fitK_countG);
    
}// ~HMMProblemPiAGK

NUMBER** HMMProblemPiAGK::getPI() { // same as getPIk
	return this->PI;
}

NUMBER** HMMProblemPiAGK::getPIk() {
	return this->PI;
}

NUMBER** HMMProblemPiAGK::getPIg() {
	return this->PIg;
}

NUMBER*** HMMProblemPiAGK::getA() { // same as getPIk
	return this->A;
}

NUMBER*** HMMProblemPiAGK::getAk() {
	return this->A;
}

NUMBER*** HMMProblemPiAGK::getAg() {
	return this->Ag;
}

NUMBER* HMMProblemPiAGK::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER** HMMProblemPiAGK::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemPiAGK::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER* HMMProblemPiAGK::getPIk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER* HMMProblemPiAGK::getPIg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing PI_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->PIg[x];
}

NUMBER** HMMProblemPiAGK::getAk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemPiAGK::getAg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->Ag[x];
}

NUMBER HMMProblemPiAGK::getPI(struct data* dt, NPAR i) {
    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
    return 1/( 1 + (1-p)*(1-q)/(p*q) );
//    return sigmoid( logit( this->PI[dt->k][i] ) + logit( this->PIg[dt->g][i] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiAGK::getA(struct data* dt, NPAR i, NPAR j) {
    NUMBER p = this->A[dt->k][i][j], q = this->Ag[dt->g][i][j];
    return 1/( 1 + (1-p)*(1-q)/(p*q) );
//    return sigmoid( logit( this->A[dt->k][i][j] ) + logit( this->Ag[dt->k][i][j] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiAGK::getB(struct data* dt, NPAR i, NPAR m) {
    return this->B[dt->k][i][m];
}

//NUMBER time_a9 = 0;

void HMMProblemPiAGK::initGrad() {
    // fitK,G - masks
    // PIk
    if( this->gradPIk == NULL )
        this->gradPIk = init2DNumber(this->p->nK, this->p->nS);
    else
        toZero2DNumber(this->gradPIk, this->p->nK, this->p->nS);
    // PIg
    if( this->gradPIg == NULL )
        this->gradPIg = init2DNumber(this->p->nG, this->p->nS);
    else
        toZero2DNumber(this->gradPIg, this->p->nG, this->p->nS);
    // Ak
    if( this->gradAk == NULL ) {
        this->gradAk  = init3DNumber(this->p->nK, this->p->nS, this->p->nS);
    }
    else
        toZero3DNumber(this->gradAk, this->p->nK, this->p->nS, this->p->nS);
    // Ag
    if( this->gradAg == NULL ) {
        this->gradAg  = init3DNumber(this->p->nG, this->p->nS, this->p->nS);
    }
    else
        toZero3DNumber(this->gradAg, this->p->nG, this->p->nS, this->p->nS);
    // B
    if( this->gradB == NULL ) {
        this->gradB  = init3DNumber(this->p->nK, this->p->nS, this->p->nO);
    }
    else
        toZero3DNumber(this->gradB, this->p->nK, this->p->nS, this->p->nO);
}

void HMMProblemPiAGK::computeGradients() {
	initGrad();
	NCAT x, k, g;
	NDAT t, xndat;
	NPAR i, j, o;
    struct data ** x_data = NULL;
    NUMBER update, combined_pLo;
    for(k=0; k<this->p->nK; k++) { // for all sKills
        xndat  = this->p->k_numg[k];
        x_data = this->p->k_g_data[k];
        if( !this->fitK[k] && this->fitK_countG[k]==0) // if fitK flag NOT raised, and ALL of relevant fitG flags are NOT raised
            continue;
        computeAlphaAndPOParam(xndat, x_data);
        computeBeta(xndat, x_data);
        for(x=0; x<xndat; x++) { // for all Groups within this sKill
            if( x_data[x]->cnt!=0 ) continue; // this is an older superceding mechanism
            t = 0;
            o = x_data[x]->obs[t];
            g = x_data[x]->g;
            if(this->fitK[k] || this->fitG[g]) { // if fitK or fitG are raised
                for(i=0; i<this->p->nS; i++) {
                    // logit combination
                    combined_pLo = sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
                    update =  combined_pLo * (1-combined_pLo) * x_data[x]->beta[t][i] * B[k][i][o] / safe0num(x_data[x]->p_O_param);
                    this->gradPIk[k][i] -= update * this->fitK[k];
                    this->gradPIg[g][i] -= update * this->fitG[g];
                    // Corbett combination - todo
                }
            }
            // if fitK not raised - continue
            if( !this->fitK[k] )
                continue;
            for(t=0; t<x_data[x]->ndat; t++) {
                o = x_data[x]->obs[t];
                // Gradient with respect to A
                // \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
                if( t>0 ) {
                    for(i=0; i<this->p->nS /*&& param->fitparam[1]>0*/; i++)
                        for(j=0; j<this->p->nS; j++) {
                            this->gradAk[k][i][j] -= x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                            this->gradAg[g][i][j] -= x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                        }
                }// if not first obs in sequence
                // Gradient with respect to B
                for(i=0; i<this->p->nS /*&& param->fitparam[2]>0*/; i++)
                    this->gradB[k][i][o] -= x_data[g]->alpha[t][i] * x_data[g]->beta[t][i] / safe0num(x_data[g]->p_O_param * B[k][i][o]);
            } // for all observations within skill-group
        } // for all Groups
        RecycleFitData(xndat, x_data, this->p);
    }// for all sKills
    
} // computeGradients()

void HMMProblemPiAGK::computeGradientsK(NCAT k, NUMBER* a_gradPI, NUMBER** a_gradA, NUMBER** a_gradB) {
	toZero1DNumber(a_gradPI, this->p->nS);
	toZero2DNumber(a_gradA , this->p->nS, this->p->nS);
	toZero2DNumber(a_gradB , this->p->nS, this->p->nO);
    
	NCAT x, g;
	NDAT t, xndat;
	NPAR i, j, o;
    struct data ** x_data = NULL;
    NUMBER update, combined_PI, combined_A;
    xndat  = this->p->k_numg[k];
    x_data = this->p->k_g_data[k];
    
    computeAlphaAndPOParam(xndat, x_data);
    computeBeta(xndat, x_data);
    for(x=0; x<xndat; x++) { // for all Groups within this sKill
        if( x_data[x]->cnt!=0 ) continue; // this is an older superceding mechanism
        t = 0;
        o = x_data[x]->obs[t];
        g = x_data[x]->g;
        for(i=0; i<this->p->nS; i++) {
            // logit combination
            combined_PI = sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            update =  combined_PI * (1-combined_PI) * x_data[x]->beta[t][i] * this->B[k][i][o] / safe0num(x_data[x]->p_O_param);
            a_gradPI[i] -= update;
            //this->gradPIg[gidx][i] -= update * this->fitG[g];
            // Corbett combination - todo
        }
        for(t=0; t<x_data[x]->ndat; t++) {
            o = x_data[x]->obs[t];
            // Gradient with respect to A
            // \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
            if( t>0 ) {
                for(i=0; i<this->p->nS /*&& param->fitparam[1]>0*/; i++)
                    for(j=0; j<this->p->nS; j++) {
                        combined_A = sigmoid( logit(this->A[k][i][j] + logit(this->Ag[g][i][j] ) ) );
                        update =  combined_A * (1-combined_A) * x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                        a_gradA[i][j] -= update;
                    }
            }// if not first obs in sequence
            // Gradient with respect to B
            for(i=0; i<this->p->nS /*&& param->fitparam[2]>0*/; i++)
                a_gradB[i][o] -= x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * this->B[k][i][o]);
        } // for all observations within skill-group
    } // for all Groups
    RecycleFitData(xndat, x_data, this->p);
    
} // computeGradients()

void HMMProblemPiAGK::computeGradientsG(NCAT g, NUMBER* a_gradPI, NUMBER** a_gradA) {
	toZero1DNumber(a_gradPI, this->p->nS);
	toZero2DNumber(a_gradA , this->p->nS, this->p->nS);
    
	NCAT x, k;
	NDAT t, xndat;
	NPAR i, j, o;
    struct data ** x_data = NULL;
    NUMBER update, combined_PI, combined_A;
    xndat  = this->p->g_numk[g];
    x_data = this->p->g_k_data[g];
    computeAlphaAndPOParam(xndat, x_data);
    computeBeta(xndat, x_data);
    for(x=0; x<xndat; x++) { // for all Groups within this sKill
        if( x_data[x]->cnt!=0 ) continue; // this is an older superceding mechanism
        t = 0;
        o = x_data[x]->obs[t];
        k = x_data[x]->k;
        for(i=0; i<this->p->nS; i++) {
            // logit combination
            combined_PI = sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            update =  combined_PI * (1-combined_PI) * x_data[x]->beta[t][i] * B[k][i][o] / safe0num(x_data[x]->p_O_param);
            a_gradPI[i] -= update;
            // Corbett combination - todo
        }
        for(t=0; t<x_data[x]->ndat; t++) {
            o = x_data[x]->obs[t];
            // Gradient with respect to A
            // \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
            if( t>0 ) {
                for(i=0; i<this->p->nS /*&& param->fitparam[1]>0*/; i++)
                    for(j=0; j<this->p->nS; j++) {
                        combined_A = sigmoid( logit(this->Ag[g][i][j] ) + logit(this->A[k][i][j]) );
                        update =  combined_A * (1-combined_A) * x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
                        a_gradA[i][j] -= update;
                    }
            }// if not first obs in sequence
        } // for all observations within skill-group
    } // for all Groups
    RecycleFitData(xndat, x_data, this->p);
} // computeGradients()

void HMMProblemPiAGK::toFile(const char *filename) {
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
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
		fprintf(fid,"Ag\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%10.8f%s",this->Ag[g][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PI[k][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"Ak\t");
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

void HMMProblemPiAGK::fit() {
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    switch(this->p->solver)
    {
        case BKT_GD_APIgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
            loglik_rmse[0] += GradientDescentPLoSKillGroupOtherSkill();
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemPiAGK::GradientDescentPLoSKillGroupOtherSkill() {
	NCAT k, g;
    NUMBER loglik = 0;
    bool conv_flagsK[3] = {true, true,  true};
    bool conv_flagsG[3] = {true, true, false};
    
	bool conv;
	int iter; // iteration count
    NCAT x;
    if(this->fitK == NULL) {
        this->fitK = Malloc(bool, this->p->nK);
        this->fitK_countG = Malloc(NCAT, this->p->nK);
        for(x=0; x<this->p->nK; x++) this->fitK[x] = true;
        for(x=0; x<this->p->nK; x++) this->fitK_countG[x] = this->p->k_numg[k];
    }
    if(this->fitG == NULL) {
        this->fitG = Malloc(bool, this->p->nG);
        for(x=0; x<this->p->nG; x++) this->fitG[x] = true; //
    }
	
	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill==1) {
        NUMBER *gradPI = init1DNumber(this->p->nS);
        NUMBER **gradA = init2DNumber(this->p->nS,this->p->nS);
        NUMBER **gradB = init2DNumber(this->p->nS,this->p->nS);
        NUMBER *PI_m1 = init1DNumber(this->p->nS);			// value on previous iteration
        NUMBER **A_m1 = init2DNumber(this->p->nS,this->p->nS);
        NUMBER **B_m1 = init2DNumber(this->p->nS,this->p->nO);
        NUMBER pO0, pO;
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
        for(x=0; x<this->p->nG; x++) this->fitG[x] = false; // do not fit group gradients
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1) {
				toZero1DNumber(gradPI, this->p->nS);
				toZero2DNumber(gradA,  this->p->nS, this->p->nS);
				toZero2DNumber(gradB,  this->p->nS, this->p->nO);
			}
            // all gradients are computed in one call
            HMMProblemPiAGK::computeGradients();
			// add gradients
			for(k=0; k<this->p->nK; k++) {
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(this->gradPIk[k], gradPI, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradAk[k],  gradA,  this->p->nS, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradB[k],  gradB,  this->p->nS, this->p->nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
			cpy1DNumber(this->getPI(0), PI_m1, this->p->nS);
			cpy2DNumber(this->getA(0),  A_m1,  this->p->nS, this->p->nS);
			cpy2DNumber(this->getB(0),  B_m1,  this->p->nS, this->p->nO);
			
			// make step only on K
			for(k=0; k<this->p->nK; k++) {
                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], this->gradPIk[k], this->gradAk[k], this->gradB[k], gradPI, gradA, gradB);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(this->getPI(k), this->getA(k), this->getB(k), PI_m1, A_m1, B_m1, conv_flagsK);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                for(k=0; k<this->p->nK; k++) {
                    computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
        for(x=0; x<this->p->nG; x++) this->fitG[x] = true; // from now DO fit group gradients
        free(PI_m1);
        free2DNumber(A_m1, this->p->nS);
        free2DNumber(B_m1, this->p->nS);
        free(gradPI);
        free2DNumber(gradA, this->p->nS);
        free2DNumber(gradB, this->p->nS);
	}
	
	//
	// Main fit
	//
    int first_iteration_qualify = 0; // at what iteration, qualification for skill/group convergence should start
    int iterations_to_qualify = 2; // how many concecutive iterations necessary for skill/group to qualify as converged
    NPAR* iter_qual_skill = Calloc(NPAR, this->p->nK);
    NPAR* iter_qual_group = Calloc(NPAR, this->p->nG);
    int skip_k = 0, skip_g = 0;
    iter = 1; // iteration count
    NUMBER pO, pO0;
    
    NUMBER* PI;// pointer
    NUMBER** A;// pointer
    NUMBER** B;// pointer
    NUMBER* PI_m1 = init1DNumber(this->p->nS);
    NUMBER** A_m1 = init2DNumber(this->p->nS, this->p->nS);
    NUMBER** B_m1 = init2DNumber(this->p->nS, this->p->nO);
    NUMBER* a_gradPI = init1DNumber(this->p->nS);
    NUMBER** a_gradA = init2DNumber(this->p->nS, this->p->nS);
    NUMBER** a_gradB = init2DNumber(this->p->nS, this->p->nO);
    int i = 0;
    NUMBER e;// step
    while(skip_k<this->p->nK || skip_g<this->p->nG) {
        //
        // Skills first
        //
        for(k=0; k<this->p->nK && skip_k<this->p->nK; k++) { // for all A,B-by-skill
            if(iter_qual_skill[k]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->k_numg[k];
            struct data** x_data = this->p->k_g_data[k];
            conv = 0; // converged
            iter = 1; // iteration count
            pO0 = 0.0;
            pO = 0.0;
            PI = this->PI[k];// pointer stays same through fitting
            A  = this->A[k]; // pointer stays same through fitting
            B  = this->B[k]; // pointer stays same through fitting
            toZero1DNumber(a_gradPI, this->p->nS);
            toZero2DNumber(a_gradA,  this->p->nS, this->p->nS);
            toZero2DNumber(a_gradB,  this->p->nS, this->p->nO);
            
            while( !conv && iter<=this->p->maxiter ) {
                computeGradientsK(k, a_gradPI, a_gradA, a_gradB);
                
                if(iter==1) {
                    pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
//                    if(k==4) {
//                        int z = 0;
//                        NUMBER *metrics =init1DNumber(5);
//                        for(i=0; i<this->p->nS; i++)
//                            printf("%s%6.4f%s",(i==0)?"PI: ":"",PI[i],(i==0)?"  ":"\n");
//                        for(i=0; i<this->p->nS; i++)
//                            for(NPAR j=0; j<this->p->nS; j++)
//                                printf("%s%6.4f%s",(i==0 && j==0)?" A: ":"",A[i][j],(i==1 && j==1)?"\n":"  ");
//                        for(i=0; i<this->p->nS; i++)
//                            for(NPAR j=0; j<this->p->nO; j++)
//                                printf("%s%6.4f%s",(i==0 && j==0)?" B: ":"",B[i][j],(i==1 && j==1)?"\n":"  ");
//                        for(x=0; x<xndat; x++) {
//                            metrics[0] = 0.0;
//                            HMMProblemPiAGK::computeLogLikRMSE(metrics, false, 1, &x_data[x], this->PI[k], this->PIg, this->A[k], this->Ag, this->B[k], this->p);
//                            printf("predict_i LL=%15.7f, alfa_i LL=%15.7f, (pOscaled=%28.22f, ndat=%4d)\n",metrics[0],-safelog(x_data[x]->p_O_param),x_data[x]->p_O_param*1000000000, x_data[x]->ndat);
//                        }
//                        HMMProblemPiAGK::computeLogLikRMSE(metrics, false, this->p->k_numg[k], this->p->k_g_data[k], this->PI[k], this->PIg, this->A[k], this->Ag, this->B[k], this->p);
//                        free(metrics);
//                    }
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy1DNumber(PI, PI_m1, this->p->nS);
                cpy2DNumber(A,  A_m1,  this->p->nS, this->p->nS);
                cpy2DNumber(B,  B_m1,  this->p->nS, this->p->nO);
                
                e = doLinearStepSkill(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
                
                // check convergence
                conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flagsK);
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1/*e<=this->p->tol*/ || skip_g==this->p->nG) { // converged quick, or don't care (others all converged
                    iter_qual_skill[k]++;
                    if(iter_qual_skill[k]==iterations_to_qualify || skip_g==this->p->nG) {// criterion met, or don't care (others all converged)
                        if(skip_g==this->p->nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                        skip_k++;
//                        report
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                            printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,iter,pO0,pO,conv);
                        }
                    }
                }
                else
                    iter_qual_skill[k]=0;
            }
            // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
            RecycleFitData(xndat, x_data, this->p);
        } // for all A,B-by-skill
        
        //
        // PIg second
        //
        for(g=0; g<this->p->nG && skip_g<this->p->nG; g++) { // for all PI-by-user
            if(iter_qual_group[g]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->g_numk[g];
            struct data** x_data = this->p->g_k_data[g];
            conv = 0; // converged
            iter = 1; // iteration count
            pO0 = 0.0;
            pO = 0.0;
            PI = this->PIg[g]; // pointer stays same through fitting
            A = this->Ag[g]; // pointer stays same through fitting
            toZero1DNumber(a_gradPI, this->p->nS); // zero these, for the sake of selective fitting
            toZero2DNumber(a_gradA,  this->p->nS, this->p->nS);
            while( !conv && iter<=this->p->maxiter ) {
                computeGradientsG(g, a_gradPI, a_gradA);
                if(iter==1) {
                    pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy1DNumber(PI, PI_m1, this->p->nS);
                cpy2DNumber(A,  A_m1,  this->p->nS, this->p->nS);
                
                e = doLinearStepPLoGroup(xndat, x_data, PI, A, a_gradPI, a_gradA);
                
                // check convergence
                conv = checkConvergence(PI, A, NULL, PI_m1, A_m1, NULL, conv_flagsG);
                
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1/*e<=this->p->tol*/ || skip_k==this->p->nK) { // converged quick, or don't care (others all converged
                    iter_qual_group[g]++;
                    if(iter_qual_group[g]==iterations_to_qualify || skip_k==this->p->nK) {// criterion met, or don't care (others all converged)
                        if(skip_k==this->p->nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                        skip_g++;
                        //                        report
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                            printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,iter,pO0,pO,conv);
                        }
                    }
                }
                else
                    iter_qual_group[g]=0;
            }
            // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
            RecycleFitData(xndat, x_data, this->p);
        } // for all PI-by-user
        i++; // run counter
    }// external loop for alternating PI and A,B
    free(PI_m1);
    free2DNumber(A_m1, this->p->nS);
    free2DNumber(B_m1, this->p->nS);
    free(a_gradPI);
    free2DNumber(a_gradA, this->p->nS);
    free2DNumber(a_gradB, this->p->nS);
    
    // compute loglik
//    printf("\n\n");
//    /*for(k=0; k<this->p->nK; k++) */ k = 4; { // for all A,B-by-skill
//        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
//        pO = HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
//        NUMBER *metrics =init1DNumber(5);
//        loglik +=pO*(pO>0);
//        computeLogLikRMSE(metrics, false, this->p->k_numg[k], this->p->k_g_data[k]);
//        printf("predict LL=%15.7f, alfa LL=%15.7f\n",metrics[0],pO);
//        RecycleFitData(this->p->k_numg[k], this->p->k_g_data[k], this->p);
//        free(metrics);
//    }
//    printf("fit loglik=%15.7f\n",loglik);
    return loglik;
}

NUMBER HMMProblemPiAGK::doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
                                       NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
	NPAR i,j,m;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, this->p->nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  this->p->nS, this->p->nS);
	doLog10Scale2DGentle(a_gradB,  a_B,  this->p->nS, this->p->nO);
	
	NUMBER *PI_cpy = init1DNumber(this->p->nS); // safe copy
	NUMBER **A_cpy = init2DNumber(this->p->nS, this->p->nS); // safe copy
	NUMBER **B_cpy = init2DNumber(this->p->nS, this->p->nO); // safe copy
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy1DNumber(a_PI, PI_cpy, this->p->nS); // save copy
	cpy2DNumber(a_A,  A_cpy,  this->p->nS, this->p->nS); // save copy
	cpy2DNumber(a_B,  B_cpy,  this->p->nS, this->p->nO); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<this->p->nS; i++)
	{
		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
		for(j=0; j<this->p->nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
		for(m=0; m<this->p->nO; m++) p_k_by_neg_p_k -= a_gradB[i][m]*a_gradB[i][m];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<this->p->nS; i++) {
			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
			for(j=0; j<this->p->nS; j++)
				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
			for(m=0; m<this->p->nO; m++)
				a_B[i][m] = B_cpy[i][m] - e * a_gradB[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			projectsimplex(a_PI, this->p->nS);
			for(i=0; i<this->p->nS; i++) {
				projectsimplex(a_A[i], this->p->nS);
				projectsimplex(a_B[i], this->p->nS);
			}
		} else {
			projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), this->p->nS);
			for(i=0; i<this->p->nS; i++) {
				projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], this->p->nS);
				projectsimplexbounded(a_B[i], this->getLbB()[i], this->getUbB()[i], this->p->nS);
			}
		}
		// recompute alpha and p(O|param)
        //		zeroLabels(xndat, x_data); // THIS IS NOT DONE HERE
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
        cpy1DNumber(PI_cpy, a_PI, this->p->nS); // save copy
        cpy2DNumber(A_cpy,  a_A,  this->p->nS, this->p->nS); // save copy
        cpy2DNumber(B_cpy,  a_B,  this->p->nS, this->p->nO); // save copy
    }
    RecycleFitData(xndat, x_data, this->p);
	free(PI_cpy);
	free2DNumber(A_cpy, this->p->nS);
	free2DNumber(B_cpy, this->p->nS);
    return e;
} // doLinearStep

NUMBER HMMProblemPiAGK::doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER *a_gradPI, NUMBER **a_gradA) {
	NPAR i, j;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, this->p->nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  this->p->nS, this->p->nS);
	
	NUMBER *PI_cpy = init1DNumber(this->p->nS); // safe copy
	NUMBER **A_cpy = init2DNumber(this->p->nS, this->p->nS); // safe copy
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy1DNumber(a_PI, PI_cpy, this->p->nS); // save copy
	cpy2DNumber(a_A,  A_cpy,  this->p->nS, this->p->nS); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<this->p->nS; i++) {
		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
		for(j=0; j<this->p->nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
    }
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<this->p->nS; i++) {
			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
			for(j=0; j<this->p->nS; j++)
				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			projectsimplex(a_PI, this->p->nS);
			for(i=0; i<this->p->nS; i++)
				projectsimplex(a_A[i], this->p->nS);
		} else {
			projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), this->p->nS);
			for(i=0; i<this->p->nS; i++)
				projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], this->p->nS);
		}
		// recompute alpha and p(O|param)
        //		zeroLabels(xndat, x_data); // THIS IS NOT DONE HERE
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
        cpy1DNumber(PI_cpy, a_PI, this->p->nS); // save copy
        cpy2DNumber(A_cpy,  a_A,  this->p->nS, this->p->nS); // save copy
    }
	free(PI_cpy);
	free2DNumber(A_cpy, this->p->nS);
    RecycleFitData(xndat, x_data, this->p);
    return e;
} // doLinearStep
