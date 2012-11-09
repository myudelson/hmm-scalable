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
    this->sizes[0] = param->nK;
    this->sizes[1] = param->nK;
    this->sizes[2] = param->nK;
    this->n_params = param->nK * 4 + param->nG * 2;
    init(param);
}

void HMMProblemPiAGK::init(struct param *param) {
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;

	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, nO);
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
    
    // mass produce PI's/PIg's, A's, B's
	if( true /*checkPIABConstraints(a_PI, a_A, a_B)*/ ) {
		this->PI  = init2DNumber(nK, nS);
		this->A   = init3DNumber(nK, nS, nS);
		this->B   = init3DNumber(nK, nS, nO);
		this->PIg = init2DNumber(nG, nS);
		this->Ag  = init3DNumber(nG, nS, nS);
        NCAT x;
		for(x=0; x<nK; x++) {
			cpy1DNumber(a_PI, this->PI[x], nS);
			cpy2DNumber(a_A,  this->A[x],  nS, nS);
			cpy2DNumber(a_B,  this->B[x],  nS, nO);
        }
        // PIg start with "no-effect" params,PI[i] = 1/nS
		for(x=0; x<nG; x++) {
            for(i=0; i<nS; i++) {
                this->PIg[x][i] = (NUMBER)1/nS;
                for(j=0; j<nO; j++)
                    this->Ag[x][i][j] = (NUMBER)1/nS;
            }
        }
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
//	this->gradPIg = NULL;
//	this->gradA = NULL;
//	this->gradAg = NULL;
//	this->gradB = NULL;
    this->fitK = NULL;
    this->fitG = NULL;
    fitK_countG = NULL;
}

HMMProblemPiAGK::~HMMProblemPiAGK() {
    destroy();
}

void HMMProblemPiAGK::destroy() {
	// destroy model data
	free2DNumber(this->PIg, this->p->nG);
	free3DNumber(this->Ag, this->p->nG, this->p->nS);
	// destroy fitting data
//	// gradPI,g
//	if ( this->gradPIg != NULL) {
//		free2DNumber(this->gradPIg, this->p->nG);
//		this->gradPIg = NULL;
//	}
//	// gradA
//	if ( this->gradAg != NULL) {
//		free3DNumber(this->gradAg, this->p->nG, this->p->nS);
//		this->gradAg = NULL;
//	}
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

void HMMProblemPiAGK::computeGradientsK(NCAT k, NUMBER* a_gradPI, NUMBER** a_gradA, NUMBER** a_gradB) {
	toZero1DNumber(a_gradPI, this->p->nS);
	toZero2DNumber(a_gradA , this->p->nS, this->p->nS);
	toZero2DNumber(a_gradB , this->p->nS, this->p->nO);
    
	NCAT x, g;
	NDAT t, xndat;
	NPAR i, j, o;
    struct data ** x_data = NULL;
    NUMBER update, combined_PI, combined_A, deriv_logit_K;
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
            deriv_logit_K = 1 / safe0num( this->PI[k][i] * (1-this-> PI[k][i]) );
            update =  combined_PI * (1-combined_PI) * deriv_logit_K * x_data[x]->beta[t][i] * this->B[k][i][o] / safe0num(x_data[x]->p_O_param);
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
                        deriv_logit_K = 1 / safe0num( this->A[k][i][j] * (1-this->A[k][i][j]) );
                        update =  combined_A * (1-combined_A) * deriv_logit_K * x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
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
    NUMBER update, combined_PI, combined_A, deriv_logit_G;
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
            deriv_logit_G = 1 / safe0num( this->PIg[g][i] * (1-this-> PIg[g][i]) );
            update =  combined_PI * (1-combined_PI) * deriv_logit_G * x_data[x]->beta[t][i] * B[k][i][o] / safe0num(x_data[x]->p_O_param);
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
                        deriv_logit_G = 1 / safe0num( this->Ag[g][i][j] * (1-this->Ag[g][i][j]) );
                        update =  combined_A * (1-combined_A) * deriv_logit_G * x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
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
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    NUMBER loglik = 0;
    bool conv_flagsK[3] = {true, true,  true};
    bool conv_flagsG[3] = {true, true, false};
    
	bool conv;
	int iter; // iteration count
    NCAT x;
    if(this->fitK == NULL) {
        this->fitK = Malloc(bool, nK);
        this->fitK_countG = Malloc(NCAT, nK);
        for(x=0; x<nK; x++) this->fitK[x] = true;
        for(x=0; x<nK; x++) this->fitK_countG[x] = this->p->k_numg[k];
    }
    if(this->fitG == NULL) {
        this->fitG = Malloc(bool, nG);
        for(x=0; x<nG; x++) this->fitG[x] = true; //
    }
	
    NUMBER* PI;// pointer
    NUMBER** A;// pointer
    NUMBER** B;// pointer
    NUMBER *PI_m1, ** A_m1, ** B_m1;
    init3Params(PI_m1, A_m1, B_m1, nS, nO);
    NUMBER *a_gradPI, ** a_gradA, ** a_gradB;
    init3Params(a_gradPI, a_gradA, a_gradB, nS, nO);
    NUMBER pO0, pO;
	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill==1) {
		iter = 1;
		pO0 = 0.0;
        pO = 0.0;
        PI = this->getPI(0);
        A =  this->getA(0);
        B =  this->getB(0);
		conv = 0; // converged
        NUMBER *gradPI_sum, ** gradA_sum, ** gradB_sum;
        init3Params(gradPI_sum, gradA_sum, gradB_sum, nS, nO);
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1)
                toZero3Params(gradPI_sum, gradA_sum, gradB_sum, nS, nO);
			// add gradients
			for(k=0; k<nK; k++) {
                // all gradients are computed in one call
                computeGradientsK(k, a_gradPI, a_gradA, a_gradB);
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(a_gradPI, gradPI_sum, nS, 1.0);
				add2DNumbersWeighted(a_gradA,  gradA_sum,  nS, nS, 1.0);
				add2DNumbersWeighted(a_gradB,  gradB_sum,  nS, nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
            cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
                
			// make step only on K
			for(k=0; k<nK; k++) {
                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], this->PI[k], this->A[k], this->B[k], a_gradPI, a_gradA, a_gradB);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(this->getPI(k), this->getA(k), this->getB(k), PI_m1, A_m1, B_m1, conv_flagsK);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                for(k=0; k<nK; k++) {
                    computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
        toZero3Params(PI_m1, A_m1, B_m1, nS, nO);
        free(gradPI_sum);
        free2DNumber(gradA_sum, nS);
        free2DNumber(gradB_sum, nS);
	}
	
	//
	// Main fit
	//
    int first_iteration_qualify = 0; // at what iteration, qualification for skill/group convergence should start
    int iterations_to_qualify = 2; // how many concecutive iterations necessary for skill/group to qualify as converged
    NPAR* iter_qual_skill = Calloc(NPAR, nK);
    NPAR* iter_qual_group = Calloc(NPAR, nG);
    int skip_k = 0, skip_g = 0;
    iter = 1; // iteration count
    
    int i = 0;
    NUMBER e;// step
    while(skip_k<nK || skip_g<nG) {
        //
        // Skills first
        //
        for(k=0; k<nK && skip_k<nK; k++) { // for all A,B-by-skill
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

            while( !conv && iter<=this->p->maxiter ) {
                computeGradientsK(k, a_gradPI, a_gradA, a_gradB);
                
                if(iter==1) {
                    pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                }
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
                
                e = doLinearStepSkill(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
                // check convergence
                conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flagsK);
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1/*e<=this->p->tol*/ || skip_g==nG) { // converged quick, or don't care (others all converged
                    iter_qual_skill[k]++;
                    if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
                        if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
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
        for(g=0; g<nG && skip_g<nG; g++) { // for all PI-by-user
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
            while( !conv && iter<=this->p->maxiter ) {
                computeGradientsG(g, a_gradPI, a_gradA);
                if(iter==1) {
                    pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy1DNumber(PI, PI_m1, nS);
                cpy2DNumber(A,  A_m1,  nS, nS);
                
                e = doLinearStepPLoGroup(xndat, x_data, PI, A, a_gradPI, a_gradA);
                
                // check convergence
                conv = checkConvergence(PI, A, NULL, PI_m1, A_m1, NULL, conv_flagsG);
                
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1/*e<=this->p->tol*/ || skip_k==nK) { // converged quick, or don't care (others all converged
                    iter_qual_group[g]++;
                    if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                        if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
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
    free2DNumber(A_m1, nS);
    free2DNumber(B_m1, nS);
    free(a_gradPI);
    free2DNumber(a_gradA, nS);
    free2DNumber(a_gradB, nS);
    
    // compute loglik
    //    printf("\n\n");
    //    /*for(k=0; k<nK; k++) */ k = 4; { // for all A,B-by-skill
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
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
	doLog10Scale2DGentle(a_gradB,  a_B,  nS, nO);
	
    NUMBER *PI_cpy, ** A_cpy, ** B_cpy;
    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy1DNumber(a_PI, PI_cpy, nS); // save copy
	cpy2DNumber(a_A,  A_cpy,  nS, nS); // save copy
	cpy2DNumber(a_B,  B_cpy,  nS, nO); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++)
	{
		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
		for(j=0; j<nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
		for(m=0; m<nO; m++) p_k_by_neg_p_k -= a_gradB[i][m]*a_gradB[i][m];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
			for(j=0; j<nS; j++)
				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
			for(m=0; m<nO; m++)
				a_B[i][m] = B_cpy[i][m] - e * a_gradB[i][m];
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
        cpy1DNumber(PI_cpy, a_PI, nS); // save copy
        cpy2DNumber(A_cpy,  a_A,  nS, nS); // save copy
        cpy2DNumber(B_cpy,  a_B,  nS, nO); // save copy
    }
    RecycleFitData(xndat, x_data, this->p);
	free(PI_cpy);
	free2DNumber(A_cpy, nS);
	free2DNumber(B_cpy, nS);
    return e;
} // doLinearStep

NUMBER HMMProblemPiAGK::doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER *a_gradPI, NUMBER **a_gradA) {
	NPAR i, j;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
	// first scale down gradients
	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
	
	NUMBER *PI_cpy = init1DNumber(nS); // safe copy
	NUMBER **A_cpy = init2DNumber(nS, nS); // safe copy
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy1DNumber(a_PI, PI_cpy, nS); // save copy
	cpy2DNumber(a_A,  A_cpy,  nS, nS); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<nS; i++) {
		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
		for(j=0; j<nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
    }
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<nS; i++) {
			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
			for(j=0; j<nS; j++)
				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			projectsimplex(a_PI, nS);
			for(i=0; i<nS; i++)
				projectsimplex(a_A[i], nS);
		} else {
			projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), nS);
			for(i=0; i<nS; i++)
				projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], nS);
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
        cpy1DNumber(PI_cpy, a_PI, nS); // save copy
        cpy2DNumber(A_cpy,  a_A,  nS, nS); // save copy
    }
	free(PI_cpy);
	free2DNumber(A_cpy, nS);
    RecycleFitData(xndat, x_data, this->p);
    return e;
} // doLinearStep
