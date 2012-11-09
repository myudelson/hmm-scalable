//
//  HMMProblemPiGKPloGK.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 8/31/12.
//
//

#include "HMMProblemPiGK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemPiGK::HMMProblemPiGK(struct param *param) {
    for(NPAR i=0; i<3; i++) this->sizes[i] = param->nK;
    this->n_params = param->nK * 4 + param->nG;
    init(param);
}

void HMMProblemPiGK::init(struct param *param) {
	this->p = param;
    NPAR nS = this->p->nS, nO = this->p->nO; //NCAT nK = this->p->nK, nG = this->p->nG;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, this->p->nO);
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
		this->A   = init3DNumber(this->p->nK, this->p->nS, this->p->nS);
		this->B   = init3DNumber(this->p->nK, this->p->nS, this->p->nO);
		this->PIg = init2DNumber(this->p->nG, this->p->nS);
        NCAT x;
		for(x=0; x<this->p->nK; x++) {
			cpy1DNumber(a_PI, this->PI[x], this->p->nS);
			cpy2DNumber(a_A,  this->A[x],  this->p->nS, this->p->nS);
			cpy2DNumber(a_B,  this->B[x],  this->p->nS, this->p->nO);
        }
        // PIg start with same params
		for(x=0; x<this->p->nG; x++)
			cpy1DNumber(a_PI, this->PIg[x], this->p->nS);
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
    init3Params(this->lbPI, this->lbA, this->lbB, nS, nO);
    init3Params(this->ubPI, this->ubA, this->ubB, nS, nO);
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
//	this->gradPI = NULL;
//	this->gradPIg = NULL;
//	this->gradA = NULL;
//	this->gradB = NULL;
    this->fitK = NULL;
    this->fitG = NULL;
    fitK_countG = NULL;
}

HMMProblemPiGK::~HMMProblemPiGK() {
    destroy();
}

void HMMProblemPiGK::destroy() {
	// destroy additional model data
	free2DNumber(this->PIg, this->p->nG);
	// destroy fitting data
//	// gradPI,g
//	if ( this->gradPIg != NULL) {
//		free2DNumber(this->gradPIg, this->p->nG);
//		this->gradPIg = NULL;
//	}
    // free fit flags
    if(this->fitK != NULL)        free(this->fitK);
    if(this->fitG != NULL)        free(this->fitG);
    if(this->fitK_countG != NULL) free(this->fitK_countG);    
    
}// ~HMMProblemPiGK

NUMBER** HMMProblemPiGK::getPI() { // same as getPIk
	return this->PI;
}

NUMBER** HMMProblemPiGK::getPIk() {
	return this->PI;
}

NUMBER** HMMProblemPiGK::getPIg() {
	return this->PIg;
}

NUMBER* HMMProblemPiGK::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER** HMMProblemPiGK::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemPiGK::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER* HMMProblemPiGK::getPIk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER* HMMProblemPiGK::getPIg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing PI_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->PIg[x];
}

NUMBER HMMProblemPiGK::getPI(struct data* dt, NPAR i) {
    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
    return 1/( 1 + (1-p)*(1-q)/(p*q) );
//    return sigmoid( logit( this->PI[dt->k][i] ) + logit( this->PIg[dt->g][i] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiGK::getA(struct data* dt, NPAR i, NPAR j) {
    return this->A[dt->k][i][j];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemPiGK::getB(struct data* dt, NPAR i, NPAR m) {
    return this->B[dt->k][i][m];
}

void HMMProblemPiGK::setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag){
    NDAT t = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    o = dt->obs[t];
    if(kg_flag == 0) { // k
        for(i=0; i<this->p->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( this->PI[ dt->k ][i] * (1-this->PI[ dt->k ][i]) );
			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * getB(dt,i,o) / safe0num(dt->p_O_param);
        }
    }
    else
        for(i=0; i<this->p->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( this->PIg[ dt->g ][i] * (1-this->PIg[ dt->g ][i]) );
			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * getB(dt,i,o) / safe0num(dt->p_O_param);
        }
}

void HMMProblemPiGK::toFile(const char *filename) {
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
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
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

void HMMProblemPiGK::fit() {
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    switch(this->p->solver)
    {
        case BKT_GD_PIgk: // Gradient Descent, pLo=f(K,G), other by K
            loglik_rmse[0] += GradientDescentPLoSKillGroupOtherSkill(0/*K*/);
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemPiGK::GradientDescentPLoSKillGroupOtherSkill(NPAR kg_flag) {
	NCAT k, g;
    // NPAR nS = this->p->nS, nO = this->p->nO;
    NCAT nK = this->p->nK, nG = this->p->nG;
    NUMBER loglik = 0;
//    bool conv_flagsK[3] = {true, true,  true};
//    bool conv_flagsG[3] = {true, false, false};
    
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
    
    FitBit *fb = new FitBit(this->p);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    fb->init(FBS_GRADm1);
    fb->init(FBS_DIRm1);
//    NUMBER *PI, **A, **B; // just pointer
//    NUMBER *PI_m1, ** A_m1, ** B_m1;
//    init3Params(PI_m1, A_m1, B_m1, nS, nO);
//    NUMBER *a_gradPI, ** a_gradA, ** a_gradB;
//    init3Params(a_gradPI, a_gradA, a_gradB, nS, nO);

	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill==1) {
        GradientDescent1Skill(fb);
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
    NUMBER pO, pO0;
//    NUMBER e; // step
    
    int i = 0;
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
            fb->linkPar(this->getPI(k), this->getA(k), this->getB(k));
//            PI = this->PI[k];// pointer stays same through fitting
//            A  = this->A[k]; // pointer stays same through fitting
//            B  = this->B[k]; // pointer stays same through fitting
//            toZero1DNumber(a_gradPI, nS);
//            toZero2DNumber(a_gradA,  nS, nS);
//            toZero2DNumber(a_gradB,  nS, nO);
            
            while( !conv && iter<=this->p->maxiter ) {
                computeGradients(xndat, x_data, fb, 0/*K here*/);// a_gradPI, a_gradA, a_gradB);
                if(iter==1) {
                    pO0 = HMMProblem::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                //            cpy3Params(PI, A, B, PI_m1, A_m1, B_m1, nS, nO);
                fb->copy(FBS_PAR, FBS_PARm1);
                
                doLinearStep(xndat, x_data, fb, -1/*co copy*/);//PI, A, B, a_gradPI, a_gradA, a_gradB);
                
                // check convergence
                //			conv = checkConvergence(PI, A, B, PI_m1, A_m1, B_m1, conv_flags);
                conv = fb->checkConvergence();
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1 /*e<=this->p->tol*/ || skip_g==nG) { // converged quick, or don't care (others all converged
                    iter_qual_skill[k]++;
                    if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
                        if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                        skip_k++;
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
        //        tm0 = clock();
        for(g=0; g<nG && skip_g<nG; g++) { // for all PI-by-user
            if(iter_qual_group[g]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->g_numk[g];
            struct data** x_data = this->p->g_k_data[g];
            conv = 0; // converged
            iter = 1; // iteration count
            pO0 = 0.0;
            pO = 0.0;
            fb->linkPar(this->getPIg(g), NULL, NULL);
//            PI = this->PIg[g]; // pointer stays same through fitting
//            toZero1DNumber(a_gradPI, nS); // zero these, for the sake of selective fitting
            while( !conv && iter<=this->p->maxiter ) {
                computeGradients(xndat, x_data, fb, 1/*G here*/);
                if(iter==1) {
                    pO0 = HMMProblemPiGK::getSumLogPOPara(xndat, x_data);
                }
                
                // copy old SAVED! values for params
                // no selectivity here, these should change
//                cpy1DNumber(PI, PI_m1, nS);
                fb->copy(FBS_PAR, FBS_PARm1);
                
//                e = doLinearStepPLoGroup(xndat, x_data, PI, a_gradPI);
                doLinearStep(xndat, x_data, fb, -1/*co copy*/);//PI, A, B, a_gradPI, a_gradA, a_gradB);
                
                // check convergence
//                conv = checkConvergence(PI, NULL, NULL, PI_m1, NULL, NULL, conv_flagsG);
                conv = fb->checkConvergence();
                
                iter ++;
            } // main solver loop
            iter--; // to turn in right
            // count convergence
            if(i>=first_iteration_qualify) {
                if(iter==1 /*e<=this->p->tol*/ || skip_k==nK) { // converged quick, or don't care (others all converged
                    iter_qual_group[g]++;
                    if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                        if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                        skip_g++;
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
        //		printf("time G elapsed is %8.6f seconds\n",(NUMBER)(clock()-tm0)/CLOCKS_PER_SEC);
        i++; // run counter
    }// external loop for alternating PI and A,B
//    free(PI_m1);
//    free2DNumber(A_m1, nS);
//    free2DNumber(B_m1, nS);
//    free(a_gradPI);
//    free2DNumber(a_gradA, nS);
//    free2DNumber(a_gradB, nS);
    delete fb;
    
    // compute loglik
    for(k=0; k<nK; k++) { // for all A,B-by-skill
        pO = HMMProblemPiGK::getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
        loglik +=pO*(pO>0);
    }
    return loglik;
}

//NUMBER HMMProblemPiGK::doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
//                                         NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
//	NPAR i,j,m;
//    NPAR nS = this->p->nS, nO = this->p->nO; //NCAT nK = this->p->nK, nG = this->p->nG;
//	// first scale down gradients
//	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
//	doLog10Scale2DGentle(a_gradA,  a_A,  nS, nS);
//	doLog10Scale2DGentle(a_gradB,  a_B,  nS, nO);
//	
//    NUMBER *PI_cpy, ** A_cpy, ** B_cpy;
//    init3Params(PI_cpy, A_cpy, B_cpy, nS, nO);
//	NUMBER e = this->p->ArmijoSeed; // step seed
//	bool compliesArmijo = false;
//	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
//	NUMBER f_xkplus1;
//	
//	cpy1DNumber(a_PI, PI_cpy, nS); // save copy
//	cpy2DNumber(a_A,  A_cpy,  nS, nS); // save copy
//	cpy2DNumber(a_B,  B_cpy,  nS, nO); // save copy
//	// compute p_k * -p_k
//	NUMBER p_k_by_neg_p_k = 0;
//	for(i=0; i<nS; i++)
//	{
//		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
//		for(j=0; j<nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
//		for(m=0; m<nO; m++) p_k_by_neg_p_k -= a_gradB[i][m]*a_gradB[i][m];
//	}
//	int iter = 0; // limit iter steps to 20
//	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
//		// update
//		for(i=0; i<nS; i++) {
//			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
//			for(j=0; j<nS; j++)
//				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
//			for(m=0; m<nO; m++)
//				a_B[i][m] = B_cpy[i][m] - e * a_gradB[i][m];
//		}
//		// scale
//		if( !this->hasNon01Constraints() ) {
//			projectsimplex(a_PI, nS);
//			for(i=0; i<nS; i++) {
//				projectsimplex(a_A[i], nS);
//				projectsimplex(a_B[i], nS);
//			}
//		} else {
//			projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), nS);
//			for(i=0; i<nS; i++) {
//				projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], nS);
//				projectsimplexbounded(a_B[i], this->getLbB()[i], this->getUbB()[i], nS);
//			}
//		}
//		// recompute alpha and p(O|param)
//        //		zeroLabels(xndat, x_data); // THIS IS NOT DONE HERE
//		computeAlphaAndPOParam(xndat, x_data);
//		// compute f(x_{k+1})
//		f_xkplus1 = HMMProblem::getSumLogPOPara(xndat, x_data);
//		// compute Armijo compliance
//		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
//		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
//		iter++;
//	} // armijo loop
//    if(!compliesArmijo) {
//        e = 0;
//        cpy1DNumber(PI_cpy, a_PI, nS); // save copy
//        cpy2DNumber(A_cpy,  a_A,  nS, nS); // save copy
//        cpy2DNumber(B_cpy,  a_B,  nS, nO); // save copy
//    }
//    RecycleFitData(xndat, x_data, this->p);
//	free(PI_cpy);
//	free2DNumber(A_cpy, nS);
//	free2DNumber(B_cpy, nS);
//    return e;
//} // doLinearStep
//
//NUMBER HMMProblemPiGK::doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER *a_gradPI) {
//	NPAR i;
//    NPAR nS = this->p->nS, //nO = this->p->nO; //NCAT nK = this->p->nK, nG = this->p->nG;	// first scale down gradients
//	doLog10Scale1DGentle(a_gradPI, a_PI, nS);
//	
//	NUMBER *PI_cpy = init1DNumber(nS); // safe copy
//	NUMBER e = this->p->ArmijoSeed; // step seed
//	bool compliesArmijo = false;
//	NUMBER f_xk = HMMProblemPiGK::getSumLogPOPara(xndat, x_data);
//	NUMBER f_xkplus1;
//	
//	cpy1DNumber(a_PI, PI_cpy, nS); // save copy
//	// compute p_k * -p_k
//	NUMBER p_k_by_neg_p_k = 0;
//	for(i=0; i<nS; i++)
//		p_k_by_neg_p_k -= a_gradPI[i]*a_gradPI[i];
//	int iter = 0; // limit iter steps to 20
//	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
//		// update
//		for(i=0; i<nS; i++)
//			a_PI[i] = PI_cpy[i] - e * a_gradPI[i];
//		// scale
//		if( !this->hasNon01Constraints() )
//			projectsimplex(a_PI, nS);
//		else
//			projectsimplexbounded(a_PI, this->getLbPI(), this->getUbPI(), nS);
//		// recompute alpha and p(O|param)
//        //		zeroLabels(xndat, x_data); // THIS IS NOT DONE HERE
//        computeAlphaAndPOParam(xndat, x_data);
//		// compute f(x_{k+1})
//		f_xkplus1 = HMMProblemPiGK::getSumLogPOPara(xndat, x_data);
//		// compute Armijo compliance
//		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
//		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
//		iter++;
//	} // armijo loop
//    if(!compliesArmijo) {
//        e = 0;
//        cpy1DNumber(PI_cpy, a_PI, nS); // save copy
//    }
//	free(PI_cpy);
//    RecycleFitData(xndat, x_data, this->p);
//    return e;
//} // doLinearStep
