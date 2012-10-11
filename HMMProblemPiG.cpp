//
//  HMMProblemPiG.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 8/24/12.
//
//

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include "HMMProblemPiG.h"
#include <map>

HMMProblemPiG::HMMProblemPiG(struct param *param) {
    switch (param->solver) {
        case BKT_GD_PIg: // Gradient Descent: PI by group, A,B by skill
            this->sizes[0] = param->nG;
            this->sizes[1] = param->nK;
            this->sizes[2] = param->nK;
            this->n_params = param->nG + 3 * param->nK;
            break;
        default:
            fprintf(stderr,"Method specified is not supported and should have been caught earlier\n");
            break;
    }
    init(param);
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiG::getPI(struct data* dt, NPAR i) {
    return this->PI[dt->g][i];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiG::getA(struct data* dt, NPAR i, NPAR j) {
    return this->A[dt->k][i][j];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiG::getB(struct data* dt, NPAR i, NPAR m) {
    return this->B[dt->k][i][m];
}


void HMMProblemPiG::computeGradientsPI(NCAT xndat, struct data** x_data, NUMBER *a_gradPI) {
	toZero1DNumber(a_gradPI, this->p->nS);

	computeAlphaAndPOParam(xndat, x_data);
	computeBeta(xndat, x_data);
	
	NCAT x, k;
	NDAT t;
	NPAR i, o;
	// pre-compute smallest p(O|param)
	NUMBER min_p_O_param = 1;
    //	for(x=0; x<xndat && scale_p_O_param; x++)
    //		if( x_data[x]->p_O_param < min_p_O_param)
    //			min_p_O_param = x_data[x]->p_O_param;
    
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
        k = x_data[x]->k;
		// Mask for gradPI is handled differently
		// Gradient with respect to PI
		t = 0;
		o = x_data[x]->obs[t];
		for(i=0; i<this->p->nS /*&& param->fitparam[0]>0*/; i++)
			a_gradPI[i] -= min_p_O_param * x_data[x]->beta[t][i] * getB(x_data[x],i,o) / safe0num(x_data[x]->p_O_param);
        // skip gradients for A,B
	} // for all groups in skill
	
} // computeGradients()

void HMMProblemPiG::computeGradientsAB(NCAT xndat, struct data** x_data, NUMBER** a_gradA, NUMBER **a_gradB) {
	toZero2DNumber(a_gradA , this->p->nS, this->p->nS);
	toZero2DNumber(a_gradB , this->p->nS, this->p->nO);

	computeAlphaAndPOParam(xndat, x_data);
	computeBeta(xndat, x_data);
	
	NCAT x;
	NDAT t;
	NPAR i, j, o;
	// pre-compute smallest p(O|param)
	NUMBER min_p_O_param = 1;
    //	for(x=0; x<xndat && scale_p_O_param; x++)
    //		if( x_data[x]->p_O_param < min_p_O_param)
    //			min_p_O_param = x_data[x]->p_O_param;
    
	for(x=0; x<xndat; x++) {
		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
		// Mask for gradPI is handled differently
		// Gradient with respect to PI
		t = 0;
        // skip gradient for PI
		for(t=0; t<x_data[x]->ndat; t++) {
			o = x_data[x]->obs[t];
			// Gradient with respect to A
			// \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
			if( t>0 ) {
				for(i=0; i<this->p->nS /*&& param->fitparam[1]>0*/; i++)
					for(j=0; j<this->p->nS; j++)
						a_gradA[i][j] -= min_p_O_param * x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
				
			}// if not first obs in sequence
			// Gradient with respect to B
			for(i=0; i<this->p->nS /*&& param->fitparam[2]>0*/; i++)
				a_gradB[i][o] -= min_p_O_param * x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * getB(x_data[x],i,o));
		} // for all observations within skill-group
		
	} // for all groups in skill
	
} // computeGradients()

void HMMProblemPiG::toFile(const char *filename) {
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
    for(g=0; g<this->p->nG; g++) {
        it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		fprintf(fid,"PI\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PI[g][i],(i==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
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

void HMMProblemPiG::fit() {
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*compute RMSE*/);
    switch(this->p->solver)
    {
        case BKT_GD_PIg: // pLo - per stident, pT,pS,pG - per skill - alternating
            loglik_rmse[0] += GradientDescentPLoGroupOtherSkill();
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemPiG::GradientDescentPLoGroupOtherSkill(){
	NCAT k, g;
    NUMBER loglik;
    //    this->p->fitparam[0]=0;
    //    this->p->fitparam[1]=0;
    
	NUMBER *PI = NULL; // just pointer
	NUMBER **A = NULL; //
	NUMBER **B = NULL; //
	NUMBER *PI_m1 = init1DNumber(this->p->nS);			// value on previous iteration
	NUMBER **A_m1 = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **B_m1 = init2DNumber(this->p->nS,this->p->nO);
	NUMBER *gradPI = init1DNumber(this->p->nS);
	NUMBER **gradA = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **gradB = init2DNumber(this->p->nS,this->p->nS);
	NUMBER *gradPI_m1 = init1DNumber(this->p->nS);
	NUMBER **gradA_m1 = init2DNumber(this->p->nS,this->p->nS);
	NUMBER **gradB_m1 = init2DNumber(this->p->nS,this->p->nS);
	
	int conv;
	int iter; // iteration count
	NUMBER pO0, pO;
    bool conv_flagsK[3] = {false, true,  true};
    bool conv_flagsG[3] = {true,  false, false};
	
    //	//
    //	// fit all as 1 skill first
    //	//
    //	if(this->p->single_skill==1) {
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
    //				pO0 = HMMProblem::getSumLogPOPara(param);
    //			// add gradients
    //			for(g=0; g<this->p->nG; g++)
    //				add1DNumbersWeighted(hmm->getGradPI(g), gradPI, this->p->nS, 1.0);
    //			for(k=0; k<this->p->nK; k++) {
    //				add2DNumbersWeighted(hmm->getGradA(k),  gradA,  this->p->nS, this->p->nS, 1.0);
    //				add2DNumbersWeighted(hmm->getGradB(k),  gradB,  this->p->nS, this->p->nO, 1.0);
    //			}
    //			// copy old SAVED! values for params, just for skill #0 is enough
    //            // no selectivity here, these should change
    //			cpy1DNumber(hmm->getPI(0), PI_m1, this->p->nS);
    //			cpy2DNumber(hmm->getA(0),  A_m1,  this->p->nS, this->p->nS);
    //			cpy2DNumber(hmm->getB(0),  B_m1,  this->p->nS, this->p->nO);
    //
    //			// make step
    //			for(k=0; k<this->p->nK; k++) {
    //                //				doLevinsonRabinerSondhi1982Step(hmm->getPI(k), hmm->getA(k), hmm->getB(k), gradPI, gradA, gradB, param);
    //                doLinearStepPLoGroupOther(this->p->k_numg[k], this->p->k_g_data[k], PI, A, B, gradPI, gradA, gradB, param, hmm);
    //            }
    //			// check convergence, on any skill, e.g. #0
    //			conv = checkConvergence(hmm->getPI(k), hmm->getA(k), hmm->getB(k), PI_m1, A_m1, B_m1, param);
    //
    //			if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
    //				hmm->computeAlphaAndPOParam();
    //				NUMBER pO = HMMProblem::getSumLogPOPara(param);
    //				printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",iter,pO0,pO,conv);
    //			}
    //			iter ++;
    //		}// single skill loop
    //	}
    //
    
	//
	// Main fit
	//
    // this is a loop of repeated fits of PI(user), and A,B(skill)
    NPAR* mask_skill = Calloc(NPAR, this->p->nK);
    NPAR* mask_group = Calloc(NPAR, this->p->nG);
    int first_iteration_qualify = 0; // at what iteration, qualification for skill/group convergence should start
    int iterations_to_qualify = 2; // how many concecutive iterations necessary for skill/group to qualify as converged
    int skip_k = 0, skip_g = 0;
    for(int i=0; i<10; i++) {
        //
        // PI first
        //
        for(g=0; g<this->p->nG && skip_k<this->p->nK && skip_g<this->p->nG; g++) { // for all PI-by-user
            if(mask_group[g]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->g_numk[g];
            struct data** x_data = this->p->g_k_data[g];
            conv = 0; // converged
            iter = 1; // iteration count
            pO0 = 0.0;
            pO = 0.0;
            PI = HMMProblem::getPI(g); // pointer stays same through fitting
            toZero1DNumber(gradPI_m1, this->p->nS); // zero these, for the sake of selective fitting
            //
            // *** FINISHED HERE
            //
            while( !conv && iter<=this->p->maxiter ) {
                if(iter>1) cpy1DNumber(gradPI, gradPI_m1, this->p->nS);
                
                computeGradientsPI(xndat,x_data, gradPI);
                if(iter==1) {
                    pO0 = HMMProblemPiG::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy1DNumber(PI, PI_m1, this->p->nS);
                
                doLinearStepPLoGroupOtherPI(xndat, x_data, PI, gradPI);
                
                // check convergence
                conv = checkConvergence(PI, NULL, NULL, PI_m1, NULL, NULL, conv_flagsG);
                
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1) {
                    mask_group[g]++;
                    if(mask_group[g]==iterations_to_qualify) {
                        skip_g++;
                        // report
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                            printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,iter,pO0,pO,conv);
                        }
                    }
                }
                else
                    mask_group[g]=0;
            }
            // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
            RecycleFitData(xndat, x_data, this->p);
        } // for all PI-by-user
        
        //
        // A,B next
        //
        for(k=0; k<this->p->nK; k++) { // for all A,B-by-skill
            if(mask_skill[k]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->k_numg[k];
            struct data** x_data = this->p->k_g_data[k];
            conv = 0; // converged
            iter = 1; // iteration count
            pO0 = 0.0;
            pO = 0.0;
            A = HMMProblem::getA(k); // pointer stays same through fitting
            B = HMMProblem::getB(k); // pointer stays same through fitting
            toZero2DNumber(gradA_m1,  this->p->nS, this->p->nS);
            toZero2DNumber(gradB_m1,  this->p->nS, this->p->nO);
            //
            // *** FINISHED HERE
            //
            while( !conv && iter<=this->p->maxiter ) {
                if(iter>1) {
                    cpy2DNumber(gradA,  gradA_m1,  this->p->nS, this->p->nS);
                    cpy2DNumber(gradB,  gradB_m1,  this->p->nS, this->p->nO);
                }
                
                computeGradientsAB(xndat, x_data, gradA, gradB);
                if(iter==1) {
                    pO0 = HMMProblemPiG::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
                cpy2DNumber(A,  A_m1,  this->p->nS, this->p->nS);
                cpy2DNumber(B,  B_m1,  this->p->nS, this->p->nO);
                
                doLinearStepPLoGroupOtherAB(xndat, x_data,  A, B, gradA, gradB);
                
                // check convergence
                conv = checkConvergence(NULL, A, B, NULL, A_m1, B_m1, conv_flagsK);
                
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1) {
                    mask_skill[k]++;
                    if(mask_skill[k]==iterations_to_qualify) {
                        skip_k++;
                        // report
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (conv || iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            pO = HMMProblem::getSumLogPOPara(xndat, x_data);
                            printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,iter,pO0,pO,conv);
                        }
                    }
                }
                else
                    mask_skill[k]=0;
            }
            // recycle memory (Alpha, Beta, p_O_param, Xi, Gamma)
            RecycleFitData(xndat, x_data, this->p);
        } // for all A,B-by-skill
        
    }// external loop for alternating PI and A,B
    
    // compute loglik
    for(k=0; k<this->p->nK; k++) { // for all A,B-by-skill
        pO = HMMProblemPiG::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
        loglik +=pO*(pO>0);
    }
    
    free(mask_skill);
    free(mask_group);
	free(PI_m1);
	free2DNumber(A_m1, this->p->nS);
	free2DNumber(B_m1, this->p->nS);
	free(gradPI_m1);
	free2DNumber(gradA_m1, this->p->nS);
	free2DNumber(gradB_m1, this->p->nS);
    return loglik;
}

NUMBER HMMProblemPiG::doLinearStepPLoGroupOtherPI(NCAT xndat, struct data** x_data, NUMBER *a_PI,
                                                 NUMBER *a_gradPI) {
	NPAR i;
	// first scale down gradients
	/*if(this->p->fitparam[0]>0) */doLog10Scale1D(a_gradPI, this->p->nS);
	
	NUMBER *PI_cpy = init1DNumber(this->p->nS); // safe copy
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy1DNumber(a_PI, PI_cpy, this->p->nS); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<this->p->nS; i++)
	{
		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<this->p->nS; i++) {
			/*if(this->p->fitparam[0]>0) */a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			/*if(this->p->fitparam[0]>0) */projectsimplex(a_PI, this->p->nS);
		} else {
			/*if(this->p->fitparam[0]>0) */projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), this->p->nS);
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
    }
	free(PI_cpy);
    RecycleFitData(xndat, x_data, this->p);
    return e;
} // doLinearStep

NUMBER HMMProblemPiG::doLinearStepPLoGroupOtherAB(NCAT xndat, struct data** x_data, NUMBER **a_A, NUMBER **a_B,
                                                 NUMBER **a_gradA, NUMBER **a_gradB) {
	NPAR i,j,m;
	// first scale down gradients
	doLog10Scale2D(a_gradA, this->p->nS, this->p->nS);
	doLog10Scale2D(a_gradB, this->p->nS, this->p->nO);
	
	NUMBER **A_cpy = init2DNumber(this->p->nS, this->p->nS); // safe copy
	NUMBER **B_cpy = init2DNumber(this->p->nS, this->p->nO); // safe copy
	NUMBER e = this->p->ArmijoSeed; // step seed
	bool compliesArmijo = false;
	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
	NUMBER f_xkplus1;
	
	cpy2DNumber(a_A,  A_cpy,  this->p->nS, this->p->nS); // save copy
	cpy2DNumber(a_B,  B_cpy,  this->p->nS, this->p->nO); // save copy
	// compute p_k * -p_k
	NUMBER p_k_by_neg_p_k = 0;
	for(i=0; i<this->p->nS; i++)
	{
		for(j=0; j<this->p->nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
		for(m=0; m<this->p->nO; m++) p_k_by_neg_p_k -= a_gradB[i][m]*a_gradB[i][m];
	}
	int iter = 0; // limit iter steps to 20
	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
		// update
		for(i=0; i<this->p->nS; i++) {
			for(j=0; j<this->p->nS /*&& (this->p->fitparam[1]>0)*/ ; j++)
				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
			for(m=0; m<this->p->nO /*&& (this->p->fitparam[2]>0)*/ ; m++)
				a_B[i][m] = B_cpy[i][m] - e * a_gradB[i][m];
		}
		// scale
		if( !this->hasNon01Constraints() ) {
			for(i=0; i<this->p->nS /*&& ( (this->p->fitparam[1]>0) || (this->p->fitparam[2]>0) )*/; i++) {
				/*if(this->p->fitparam[1]>0) */projectsimplex(a_A[i], this->p->nS);
				/*if(this->p->fitparam[2]>0) */projectsimplex(a_B[i], this->p->nS);
			}
		} else {
			for(i=0; i<this->p->nS /*&& ( (this->p->fitparam[1]>0) || (this->p->fitparam[2]>0) )*/; i++) {
				/*if(this->p->fitparam[1]>0) */projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], this->p->nS);
				/*if(this->p->fitparam[2]>0) */projectsimplexbounded(a_B[i], this->getLbB()[i], this->getUbB()[i], this->p->nS);
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
        cpy2DNumber(A_cpy,  a_A,  this->p->nS, this->p->nS); // save copy
        cpy2DNumber(B_cpy,  a_B,  this->p->nS, this->p->nO); // save copy
    }
	free2DNumber(A_cpy, this->p->nS);
	free2DNumber(B_cpy, this->p->nS);
    RecycleFitData(xndat, x_data, this->p);
    return e;
} // doLinearStep