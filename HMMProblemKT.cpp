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

//
//  HMMProblemKT.cpp - K-skill Transfers
//

#include "HMMProblemKT.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>

#include "liblinear/linear.h"

HMMProblemKT::HMMProblemKT(struct param *param) {
    for(NPAR i=0; i<3; i++) this->sizes[i] = param->nK;
    this->n_params = param->nK * 4 + param->nK*param->nK + param->nK;
    init(param);
}

void HMMProblemKT::init(struct param *param) {
    HMMProblem::init(param);
    this->liblinear_weights = Calloc(NUMBER, ((NDAT)param->nK)*((NDAT)param->nK) + (NDAT)param->nK);
}

HMMProblemKT::~HMMProblemKT() {
    free(this->liblinear_weights);
}// ~HMMProblemKT

// getters for computing alpha, beta, gamma
NUMBER HMMProblemKT::getPI(struct data* dt, NPAR i) {
    return this->PI[dt->k][i];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemKT::getA (struct data* dt, NPAR i, NPAR j) {
    return this->A[dt->k][i][j];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemKT::getB (struct data* dt, NPAR i, NPAR m) {
    return this->B[dt->k][i][m];
}

void HMMProblemKT::toFile(const char *filename) {
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    
    // write solved id
    writeSolvedId(fid, this->p);
    
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
	fprintf(fid,"Transfer matrix");
    for(NDAT r=0; r<( ((NDAT)this->p->nK)*((NDAT)this->p->nK) + (NDAT)this->p->nK); r++)
        fprintf(fid,"%s%10.8f",(r%(this->p->nK + 1)==0)?"\n":"\t",this->liblinear_weights[r]);
	fprintf(fid,"\n");
    
	fclose(fid);
}

void HMMProblemKT::fit() {
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*do RMSE*/);
    switch(this->p->solver)
    {
        case BKT_GD_T: // Gradient Descent, transfer K-K
            loglik_rmse[0] += GradientDescentKT2();
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

void HMMProblemKT::producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt) {
    if(this->liblinear_weights==NULL) {
        fprintf(stderr,"Transfer weights are NULL\n");
        exit(1);
    }
    NPAR m, i;
    NCAT k, kk;
    NUMBER *local_pred_inner = init1DNumber(this->p->nO);
    NUMBER *val = init1DNumber(this->p->nO);
    for(m=0; m<this->p->nO; m++) local_pred[m] = 0.0;
    for(int l=0; l<nks; l++) {
        for(m=0; m<this->p->nO; m++) local_pred_inner[m] = 0.0;
        k = ks[l];
        for(kk = 0; kk<this->p->nK; kk++) {
            dt->k = kk;
            for(m=0; m<this->p->nO; m++)
                val[m] = 0;
            for(m=0; m<this->p->nO; m++)
                for(i=0; i<this->p->nS; i++)
                    val[m] += group_skill_map[dt->g][kk][i] * getB(dt,i,m);//B[i][m];
            NUMBER weight = this->liblinear_weights[(NDAT)k*((NDAT)this->p->nK+(NDAT)1)+(NDAT)kk];
            for(m=0; m<this->p->nO; m++)
                local_pred_inner[m] += logit(val[m]) * weight;
        }
        double bias = this->liblinear_weights[ ((NDAT)k+(NDAT)1)*((NDAT)this->p->nK+(NDAT)1) - (NDAT)1 ];
        for(m=0; m<this->p->nO; m++)
            local_pred_inner[m] = sigmoid( local_pred_inner[m] + NUMBER(bias) );
        for(m=0; m<this->p->nO; m++)
            local_pred[m] += local_pred_inner[m]; // local_pred[m] = 0.0;
    }
    if(nks>1) {
        for(m=0; m<this->p->nO; m++)
            local_pred[m] /= nks;
        //            projectsimplex(local_pred, this->p->nO);
    }
    free(local_pred_inner);
    free(val);
}

NUMBER HMMProblemKT::GradientDescentKT1() {
	NCAT k;
    NUMBER loglik = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    
	NUMBER *PI = NULL; // just pointer
	NUMBER **A = NULL; //
	NUMBER **B = NULL; //
    NUMBER *PI_m1, ** A_m1, ** B_m1;
    init3Params(PI_m1, A_m1, B_m1, nS, nO);
    NUMBER *a_gradPI, ** a_gradA, ** a_gradB;
    init3Params(a_gradPI, a_gradA, a_gradB, nS, nO);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
        NUMBER *a_gradPI_sum, ** a_gradA_sum, ** a_gradB_sum;
        init3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS, nO);
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1) {
				toZero1DNumber(a_gradPI, nS);
				toZero2DNumber(a_gradA,  nS, nS);
				toZero2DNumber(a_gradB,  nS, nO);
				toZero1DNumber(a_gradPI_sum, nS);
				toZero2DNumber(a_gradA_sum,  nS, nS);
				toZero2DNumber(a_gradB_sum,  nS, nO);
			}
			// add gradients
			for(k=0; k<nK; k++) {
                computeGradients(this->p->k_numg[k], this->p->k_g_data[k], a_gradPI, a_gradA, a_gradB);
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(a_gradPI, a_gradPI_sum, nS, 1.0);
				add2DNumbersWeighted(a_gradA,  a_gradA_sum,  nS, nS, 1.0);
				add2DNumbersWeighted(a_gradB,  a_gradB_sum,  nS, nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
			cpy1DNumber(HMMProblem::getPI(0), PI_m1, nS);
			cpy2DNumber(HMMProblem::getA(0),  A_m1,  nS, nS);
			cpy2DNumber(HMMProblem::getB(0),  B_m1,  nS, nO);
			
			// make step
			for(k=0; k<nK; k++) {
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, a_gradPI_sum, a_gradA_sum, a_gradB_sum);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), PI_m1, A_m1, B_m1, conv_flags);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
                for(k=0; k<nK; k++) {
                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
        free(a_gradPI_sum);
        free2DNumber(a_gradA_sum, nS);
        free2DNumber(a_gradB_sum, nS);
        toZero1DNumber(a_gradPI, nS);
        toZero2DNumber(a_gradA,  nS, nS);
        toZero2DNumber(a_gradB,  nS, nO);
	}
	//
	// Main fit
	//
	for(k=0; k<nK; k++) {  // for(k=218; k<219; k++) { //
        NCAT xndat = this->p->k_numg[k];
        struct data** x_data = this->p->k_g_data[k];
		
		conv = 0; // converged
		iter = 1; // iteration count
		pO0 = 0.0;
        pO= 0.0;
		
		PI = HMMProblem::getPI(k);         // pointers stay same through fitting
		A  = HMMProblem::getA(k);
		B  = HMMProblem::getB(k);
		
		while( !conv && iter<=this->p->maxiter ) {
			computeGradients(xndat, x_data, a_gradPI, a_gradA, a_gradB);
			if(iter==1) {
				toZero1DNumber(a_gradPI, nS);
				toZero2DNumber(a_gradA,  nS, nS);
				toZero2DNumber(a_gradB,  nS, nO);
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
			cpy1DNumber(PI, PI_m1, nS);
			cpy2DNumber(A,  A_m1,  nS, nS);
			cpy2DNumber(B,  B_m1,  nS, nO);
            
			doLinearStep(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
            
			// check convergence
			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, nS);
                computeAlphaAndPOParam(xndat, x_data);
                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                //                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                loglik += pO*(pO>0);
                if(!this->p->quiet)
                    printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",k,iter,pO0,pO,conv);
			}
			iter ++;
		} // main solver loop
        RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
		// recycle
	} // for all skills
	free(PI_m1);
	free2DNumber(A_m1, nS);
	free2DNumber(B_m1, nS);
	free(a_gradPI);
	free2DNumber(a_gradA, nS);
	free2DNumber(a_gradB, nS);
    //
    // fit Transfer weights
    //
    
    for(NCAT k=0; k<(nK); k++) {
        NDAT idx = (NDAT)k*(((NDAT)nK)+1) + (NDAT)k;
        this->liblinear_weights[idx] = 1;
    }
    
    struct problem ll_prob;
    struct parameter ll_param;
    struct feature_node *ll_x_space;
    
    // vvvv for all skils
    // grab max mem beforehand
    NDAT l = this->p->N - this->p->N_null;
    fprintf(stdout,"Grab memory\n");
	ll_prob.y = Malloc(double,l);
	ll_prob.x = Malloc(struct feature_node *,l);
    long long elements = (long long)l * ((long long)nK + 1); // +1 more is for ?, but it's there
	ll_x_space = Malloc(struct feature_node,elements+l); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
    
    fprintf(stdout,"Create problem\n");
    createLiblinearProblem(ll_prob, ll_param, ll_x_space);
    fprintf(stdout,"Train model\n");
    struct model* model_=train(&ll_prob, &ll_param);
    if(model_->nr_feature != ((long long)nK*nK)) {
        fprintf(stderr,"Number of features in Logistic Regression is not correct\n");
        exit(1);
    }
    NUMBER multiplier = 1;
    if(model_->label[0] != 1) // inverse coefficients if first class label is not 1 (documented in LL)
        multiplier = -1;
    for(NDAT r=0; r<((long long)nK*nK + 1); r++)
        if(model_->w[r] != 0) {
            this->liblinear_weights[r] = multiplier * model_->w[r];
        }
    fprintf(stdout,"Destroy model\n");
    free_and_destroy_model(&model_);
    
    // free all afterwards
    fprintf(stdout,"Free LL data\n");
    recycleLiblinearProblem(ll_prob, ll_param, ll_x_space);
    // ^^^^ for all skils
    // ^^^^^^^^^^^^^^^^^^^^
    //
    /**/
    // report pO per group
    //    for(NCAT g=0; g<nG; g++) {
    //        //        computeAlphaAndPOParamG(g);
    //        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
    //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        printf("group %3d p(O|param)= %15.7f\n",g,pO);
    //    }
    
    return loglik;
}

NUMBER HMMProblemKT::GradientDescentKT2() {
	NCAT k;
    NUMBER loglik = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    
	NUMBER *PI = NULL; // just pointer
	NUMBER **A = NULL; //
	NUMBER **B = NULL; //
    NUMBER *PI_m1, ** A_m1, ** B_m1;
    init3Params(PI_m1, A_m1, B_m1, nS, nO);    
    NUMBER *a_gradPI, ** a_gradA, ** a_gradB;
    init3Params(a_gradPI, a_gradA, a_gradB, nS, nO);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
        NUMBER *a_gradPI_sum, ** a_gradA_sum, ** a_gradB_sum;
        init3Params(a_gradPI_sum, a_gradA_sum, a_gradB_sum, nS, nO);
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1) {
				toZero1DNumber(a_gradPI, nS);
				toZero2DNumber(a_gradA,  nS, nS);
				toZero2DNumber(a_gradB,  nS, nO);
				toZero1DNumber(a_gradPI_sum, nS);
				toZero2DNumber(a_gradA_sum,  nS, nS);
				toZero2DNumber(a_gradB_sum,  nS, nO);
			}
			// add gradients
			for(k=0; k<nK; k++) {
                computeGradients(this->p->k_numg[k], this->p->k_g_data[k], a_gradPI, a_gradA, a_gradB);
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(a_gradPI, a_gradPI_sum, nS, 1.0);
				add2DNumbersWeighted(a_gradA,  a_gradA_sum,  nS, nS, 1.0);
				add2DNumbersWeighted(a_gradB,  a_gradB_sum,  nS, nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
			cpy1DNumber(HMMProblem::getPI(0), PI_m1, nS);
			cpy2DNumber(HMMProblem::getA(0),  A_m1,  nS, nS);
			cpy2DNumber(HMMProblem::getB(0),  B_m1,  nS, nO);
			
			// make step
			for(k=0; k<nK; k++) {
                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, a_gradPI_sum, a_gradA_sum, a_gradB_sum);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), PI_m1, A_m1, B_m1, conv_flags);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
                for(k=0; k<nK; k++) {
                    //                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, nS);
                    //                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
        free(a_gradPI_sum);
        free2DNumber(a_gradA_sum, nS);
        free2DNumber(a_gradB_sum, nS);
	}
	//
	// Main fit
	//
	for(k=0; k<nK; k++) {  // for(k=218; k<219; k++) { //
        NCAT xndat = this->p->k_numg[k];
        struct data** x_data = this->p->k_g_data[k];
		
		conv = 0; // converged
		iter = 1; // iteration count
		pO0 = 0.0;
        pO= 0.0;
		
		PI = HMMProblem::getPI(k);         // pointers stay same through fitting
		A  = HMMProblem::getA(k);
		B  = HMMProblem::getB(k);
		
		while( !conv && iter<=this->p->maxiter ) {
			computeGradients(xndat, x_data, a_gradPI, a_gradA, a_gradB);
			if(iter==1) {
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
			cpy1DNumber(PI, PI_m1, nS);
			cpy2DNumber(A,  A_m1,  nS, nS);
			cpy2DNumber(B,  B_m1,  nS, nO);
            
			doLinearStep(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
            
			// check convergence
			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, nS);
                computeAlphaAndPOParam(xndat, x_data);
                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                //                pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                loglik += pO*(pO>0);
                if(!this->p->quiet)
                    printf("skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",k,iter,pO0,pO,conv);
			}
			iter ++;
		} // main solver loop
        RecycleFitData(xndat, x_data, this->p); // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
		// recycle
	} // for all skills
	free(PI_m1);
	free2DNumber(A_m1, nS);
	free2DNumber(B_m1, nS);
	free(a_gradPI);
	free2DNumber(a_gradA, nS);
	free2DNumber(a_gradB, nS);
    
    //
    // fit Transfer weights
    //
    for(NCAT k=0; k<(nK); k++) {
        NDAT idx = (NDAT)k*(((NDAT)nK)+1) + (NDAT)k;
        //        fprintf(stdout, "t.idx=%5d / %d (nK=%d)\n",idx,((NDAT)nK)*((NDAT)nK)+(NDAT)nK,nK);
        this->liblinear_weights[idx] = 1;
    }
    
    struct problem ll_prob;
    struct parameter ll_param;
    struct feature_node *ll_x_space;
    
    // vvvv for all skils
    NDAT counts[nK];
    for(NCAT k=0; k<nK; k++) {
        counts[k] = 0;
        for(NCAT gidx=0; gidx<this->p->k_numg[k]; gidx++)
            counts[k] += this->p->k_g_data[k][gidx]->n;
    }
    NDAT n_max = 0;
    for(NCAT k=0; k<nK; k++)
        n_max = (n_max<counts[k])?counts[k]:n_max;
    // grab max mem beforehand
	ll_prob.y = Malloc(double,n_max);
	ll_prob.x = Malloc(struct feature_node *,n_max);
    long long elements = (long long)n_max * ((long long)nK + 1); // +1 more is for ?, but it's there
	ll_x_space = Malloc(struct feature_node,elements+n_max); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
    
    for(NCAT k=0; k<nK; k++) { //
        ll_prob.l = counts[k];
        createLiblinearProblem(ll_prob, ll_param, ll_x_space, k);
        struct model* model_=train(&ll_prob, &ll_param);
        if(model_->nr_feature != (nK)) {
            fprintf(stderr,"Number of features in Logistic Regression is not correct\n");
            exit(1);
        }
        NUMBER multiplier = 1;
        if(model_->label[0] != 1) // inverse coefficients if first class label is not 1 (documented in LL)
            multiplier = -1;
        for(NCAT r=0; r<(nK + 1); r++)
            if(model_->w[r] != 0) {
                this->liblinear_weights[ (NDAT)k*((NDAT)nK+(NDAT)1) + (NDAT)r ] = multiplier * model_->w[r];
            }
        free_and_destroy_model(&model_);
    }
    // free all afterwards
    recycleLiblinearProblem(ll_prob, ll_param, ll_x_space);
    // ^^^^ for all skils
    // ^^^^^^^^^^^^^^^^^^^^
    //
    /**/
    // report pO per group
    //    for(NCAT g=0; g<nG; g++) {
    //        //        computeAlphaAndPOParamG(g);
    //        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
    //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        printf("group %3d p(O|param)= %15.7f\n",g,pO);
    //    }
    
    return loglik;
}


