//
//  HMMProblemKT.cpp - K-skill Transfers
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
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
    this->n_params = param->nK * 4 + param->nK*param->nK;
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
    
	NUMBER *PI = NULL; // just pointer
	NUMBER **A = NULL; //
	NUMBER **B = NULL; //
	NUMBER *PI_m1 = init1DNumber(this->p->nS);			// value on previous iteration
	NUMBER **A_m1 = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **B_m1 = init2DNumber(this->p->nS,this->p->nO);
	NUMBER *a_gradPI = init1DNumber(this->p->nS);
	NUMBER **a_gradA = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **a_gradB = init2DNumber(this->p->nS,this->p->nS);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1) {
				toZero1DNumber(a_gradPI, this->p->nS);
				toZero2DNumber(a_gradA,  this->p->nS, this->p->nS);
				toZero2DNumber(a_gradB,  this->p->nS, this->p->nO);
			}
			// add gradients
			for(k=0; k<this->p->nK; k++) {
                HMMProblem::computeGradients(this->p->k_numg[k], this->p->k_g_data[k], a_gradPI, a_gradA, a_gradB, this->p);
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(this->gradPI, a_gradPI, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradA,  a_gradA,  this->p->nS, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradB,  a_gradB,  this->p->nS, this->p->nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
			cpy1DNumber(HMMProblem::getPI(0), PI_m1, this->p->nS);
			cpy2DNumber(HMMProblem::getA(0),  A_m1,  this->p->nS, this->p->nS);
			cpy2DNumber(HMMProblem::getB(0),  B_m1,  this->p->nS, this->p->nO);
			
			// make step
			for(k=0; k<this->p->nK; k++) {
                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, a_gradPI, a_gradA, a_gradB);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), PI_m1, A_m1, B_m1, conv_flags);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
                for(k=0; k<this->p->nK; k++) {
                    //                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, this->p->nS);
                    //                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
	}
	//
	// Main fit
	//
	for(k=0; k<this->p->nK; k++) {  // for(k=218; k<219; k++) { //
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
			computeGradients(xndat, x_data, a_gradPI, a_gradA, a_gradB, this->p);
			if(iter==1) {
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
			cpy1DNumber(PI, PI_m1, this->p->nS);
			cpy2DNumber(A,  A_m1,  this->p->nS, this->p->nS);
			cpy2DNumber(B,  B_m1,  this->p->nS, this->p->nO);
            
			doLinearStep(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
            
			// check convergence
			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, this->p->nS);
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
	free2DNumber(A_m1, this->p->nS);
	free2DNumber(B_m1, this->p->nS);
	free(a_gradPI);
	free2DNumber(a_gradA, this->p->nS);
	free2DNumber(a_gradB, this->p->nS);
    //
    // fit Transfer weights
    //
    
    for(NCAT k=0; k<(this->p->nK); k++) {
        NDAT idx = (NDAT)k*(((NDAT)this->p->nK)+1) + (NDAT)k;
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
    long long elements = (long long)l * ((long long)this->p->nK + 1); // +1 more is for ?, but it's there
	ll_x_space = Malloc(struct feature_node,elements+l); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
    
    fprintf(stdout,"Create problem\n");
    createLiblinearProblem(ll_prob, ll_param, ll_x_space);
    fprintf(stdout,"Train model\n");
    struct model* model_=train(&ll_prob, &ll_param);
    if(model_->nr_feature != ((long long)this->p->nK*this->p->nK)) {
        fprintf(stderr,"Number of features in Logistic Regression is not correct\n");
        exit(1);
    }
    NUMBER multiplier = 1;
    if(model_->label[0] != 1) // inverse coefficients if first class label is not 1 (documented in LL)
        multiplier = -1;
    for(NDAT r=0; r<((long long)this->p->nK*this->p->nK + 1); r++)
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
//    for(NCAT g=0; g<this->p->nG; g++) {
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
    
	NUMBER *PI = NULL; // just pointer
	NUMBER **A = NULL; //
	NUMBER **B = NULL; //
	NUMBER *PI_m1 = init1DNumber(this->p->nS);			// value on previous iteration
	NUMBER **A_m1 = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **B_m1 = init2DNumber(this->p->nS,this->p->nO);
	NUMBER *a_gradPI = init1DNumber(this->p->nS);
	NUMBER **a_gradA = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **a_gradB = init2DNumber(this->p->nS,this->p->nS);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flags[3] = {true, true, true};
	
	//
	// fit all as 1 skill first
	//
	if(this->p->single_skill==1) {
		iter = 1;
		pO0 = 0.0;
		conv = 0; // converged
		while( !conv && iter<=this->p->maxiter ) {
			if(iter>1) {
				toZero1DNumber(a_gradPI, this->p->nS);
				toZero2DNumber(a_gradA,  this->p->nS, this->p->nS);
				toZero2DNumber(a_gradB,  this->p->nS, this->p->nO);
			}
			// add gradients
			for(k=0; k<this->p->nK; k++) {
                HMMProblem::computeGradients(this->p->k_numg[k], this->p->k_g_data[k], a_gradPI, a_gradA, a_gradB, this->p);
                if(iter==1)
                    pO0 += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
				add1DNumbersWeighted(this->gradPI, a_gradPI, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradA,  a_gradA,  this->p->nS, this->p->nS, 1.0);
				add2DNumbersWeighted(this->gradB,  a_gradB,  this->p->nS, this->p->nO, 1.0);
			}
			// copy old SAVED! values for params, just for skill #0 is enough
			cpy1DNumber(HMMProblem::getPI(0), PI_m1, this->p->nS);
			cpy2DNumber(HMMProblem::getA(0),  A_m1,  this->p->nS, this->p->nS);
			cpy2DNumber(HMMProblem::getB(0),  B_m1,  this->p->nS, this->p->nO);
			
			// make step
			for(k=0; k<this->p->nK; k++) {
                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
                doLinearStep(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, a_gradPI, a_gradA, a_gradB);
            }
			// check convergence, on any skill, e.g. #0
			conv = checkConvergence(HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), PI_m1, A_m1, B_m1, conv_flags);
			
			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                NUMBER pO = 0.0;
                for(k=0; k<this->p->nK; k++) {
                    //                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, this->p->nS);
                    //                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                    HMMProblem::computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
                    pO += HMMProblem::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
                }
				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
			}
			iter ++;
		}// single skill loop
	}
	//
	// Main fit
	//
	for(k=0; k<this->p->nK; k++) {  // for(k=218; k<219; k++) { //
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
			computeGradients(xndat, x_data, a_gradPI, a_gradA, a_gradB, this->p);
			if(iter==1) {
				pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
			}
			
			// copy old SAVED! values for params
			cpy1DNumber(PI, PI_m1, this->p->nS);
			cpy2DNumber(A,  A_m1,  this->p->nS, this->p->nS);
			cpy2DNumber(B,  B_m1,  this->p->nS, this->p->nO);
            
			doLinearStep(xndat, x_data, PI, A, B, a_gradPI, a_gradA, a_gradB);
            
			// check convergence
			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
			
			if( ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                //                HMMProblem::computeAlphaAndPOParam(xndat, x_data, PI, A, B, this->p->nS);
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
	free2DNumber(A_m1, this->p->nS);
	free2DNumber(B_m1, this->p->nS);
	free(a_gradPI);
	free2DNumber(a_gradA, this->p->nS);
	free2DNumber(a_gradB, this->p->nS);

    //
    // fit Transfer weights
    //
    for(NCAT k=0; k<(this->p->nK); k++) {
        NDAT idx = (NDAT)k*(((NDAT)this->p->nK)+1) + (NDAT)k;
        //        fprintf(stdout, "t.idx=%5d / %d (nK=%d)\n",idx,((NDAT)this->p->nK)*((NDAT)this->p->nK)+(NDAT)this->p->nK,this->p->nK);
        this->liblinear_weights[idx] = 1;
    }
    
    struct problem ll_prob;
    struct parameter ll_param;
    struct feature_node *ll_x_space;
    
    // vvvv for all skils
    NDAT counts[this->p->nK];
    for(NCAT k=0; k<this->p->nK; k++) {
        counts[k] = 0;
        for(NCAT gidx=0; gidx<this->p->k_numg[k]; gidx++)
            counts[k] += this->p->k_g_data[k][gidx]->ndat;
    }
    NDAT n_max = 0;
    for(NCAT k=0; k<this->p->nK; k++)
        n_max = (n_max<counts[k])?counts[k]:n_max;
    // grab max mem beforehand
	ll_prob.y = Malloc(double,n_max);
	ll_prob.x = Malloc(struct feature_node *,n_max);
    long long elements = (long long)n_max * ((long long)this->p->nK + 1); // +1 more is for ?, but it's there
	ll_x_space = Malloc(struct feature_node,elements+n_max); // here's the '+prob.l' grab for bias, what the hell is 'elemets++' earlier
    
    for(NCAT k=0; k<this->p->nK; k++) { //
        ll_prob.l = counts[k];
        createLiblinearProblem(ll_prob, ll_param, ll_x_space, k);
        struct model* model_=train(&ll_prob, &ll_param);
        if(model_->nr_feature != (this->p->nK)) {
            fprintf(stderr,"Number of features in Logistic Regression is not correct\n");
            exit(1);
        }
        NUMBER multiplier = 1;
        if(model_->label[0] != 1) // inverse coefficients if first class label is not 1 (documented in LL)
            multiplier = -1;
        for(NCAT r=0; r<(this->p->nK + 1); r++)
            if(model_->w[r] != 0) {
                this->liblinear_weights[ (NDAT)k*((NDAT)this->p->nK+(NDAT)1) + (NDAT)r ] = multiplier * model_->w[r];
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
    //    for(NCAT g=0; g<this->p->nG; g++) {
    //        //        computeAlphaAndPOParamG(g);
    //        computeAlphaAndPOParam(this->p->k_numg[k], this->p->k_g_data[k]);
    //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        //        pO = HMMProblem::getSumLogPOPara(this->p->g_numk[g], this->p->g_k_data[g]);
    //        printf("group %3d p(O|param)= %15.7f\n",g,pO);
    //    }
    
    return loglik;
}


