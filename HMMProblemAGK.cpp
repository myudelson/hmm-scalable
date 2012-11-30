//
//  HMMProblemAGK.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
//

#include "HMMProblemAGK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemAGK::HMMProblemAGK(struct param *param) {
    this->sizes[0] = param->nK;
    this->sizes[1] = param->nK;
    this->sizes[2] = param->nK;
    this->n_params = param->nK * 4 + param->nG;
    init(param);
}

void HMMProblemAGK::init(struct param *param) {
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    
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
		this->Ag  = init3DNumber(nG, nS, nS);
		this->B   = init3DNumber(nK, nS, nO);
        NCAT x;
		for(x=0; x<nK; x++) {
			cpy1DNumber(a_PI, this->PI[x], nS);
			cpy2DNumber(a_A,  this->A[x],  nS, nS);
			cpy2DNumber(a_B,  this->B[x],  nS, nO);
        }
        // PIg start with "no-effect" params,PI[i] = 1/nS
		for(x=0; x<nG; x++) {
            for(i=0; i<nS; i++) {
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
    //	this->gradA = NULL;
    //	this->gradB = NULL;
//    this->fitK = NULL;
//    this->fitG = NULL;
//    fitK_countG = NULL;
}

HMMProblemAGK::~HMMProblemAGK() {
    destroy();
}

void HMMProblemAGK::destroy() {
	// destroy model data
	free3DNumber(this->Ag, this->p->nG, this->p->nS);
//	// destroy fitting data
//	if ( this->gradAg != NULL) {
//		free3DNumber(this->gradAg, this->p->nG, this->p->nS);
//		this->gradAg = NULL;
//	}
    // free fit flags
//    if(this->fitK != NULL)        free(this->fitK);
//    if(this->fitG != NULL)        free(this->fitG);
//    if(this->fitK_countG != NULL) free(this->fitK_countG);
    
}// ~HMMProblemAGK

NUMBER** HMMProblemAGK::getPI() { // same as getPIk
	return this->PI;
}

NUMBER*** HMMProblemAGK::getA() { // same as getPIk
	return this->A;
}

NUMBER*** HMMProblemAGK::getAk() {
	return this->A;
}

NUMBER*** HMMProblemAGK::getAg() {
	return this->Ag;
}

NUMBER* HMMProblemAGK::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->PI[x];
}

NUMBER** HMMProblemAGK::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemAGK::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER** HMMProblemAGK::getAk(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemAGK::getAg(NCAT x) {
	if( x > (this->p->nG-1) ) {
		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
		exit(1);
	}
	return this->Ag[x];
}

NUMBER HMMProblemAGK::getPI(struct data* dt, NPAR i) {
    return this->PI[dt->k][i];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemAGK::getA(struct data* dt, NPAR i, NPAR j) {
    NUMBER p = this->A[dt->k][i][j], q = this->Ag[dt->g][i][j];
    return 1/( 1 + (1-p)*(1-q)/(p*q) );
    //    return sigmoid( logit( this->A[dt->k][i][j] ) + logit( this->Ag[dt->k][i][j] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemAGK::getB(struct data* dt, NPAR i, NPAR m) {
    return this->B[dt->k][i][m];
}

//NUMBER time_a9 = 0;

//void HMMProblemAGK::initGrad() {
//    // fitK,G - masks
//    // PI
//    if( this->gradPI == NULL )
//        this->gradPI = init2DNumber(this->p->nK, this->p->nS);
//    else
//        toZero2DNumber(this->gradPI, this->p->nK, this->p->nS);
//    // Ak
//    if( this->gradA == NULL ) {
//        this->gradA  = init3DNumber(this->p->nK, this->p->nS, this->p->nS);
//    }
//    else
//        toZero3DNumber(this->gradA, this->p->nK, this->p->nS, this->p->nS);
//    // Ag
//    if( this->gradAg == NULL ) {
//        this->gradAg  = init3DNumber(this->p->nG, this->p->nS, this->p->nS);
//    }
//    else
//        toZero3DNumber(this->gradAg, this->p->nG, this->p->nS, this->p->nS);
//    // B
//    if( this->gradB == NULL ) {
//        this->gradB  = init3DNumber(this->p->nK, this->p->nS, this->p->nO);
//    }
//    else
//        toZero3DNumber(this->gradB, this->p->nK, this->p->nS, this->p->nO);
//}

//void HMMProblemAGK::computeGradients() {
//	initGrad();
//	NCAT x, k, g;
//	NDAT t, xndat;
//	NPAR i, j, o;
//    struct data ** x_data = NULL;
//    NUMBER update;
//    for(k=0; k<this->p->nK; k++) { // for all sKills
//        xndat  = this->p->k_numg[k];
//        x_data = this->p->k_g_data[k];
//        if( !this->fitK[k] && this->fitK_countG[k]==0) // if fitK flag NOT raised, and ALL of relevant fitG flags are NOT raised
//            continue;
//        computeAlphaAndPOParam(xndat, x_data);
//        computeBeta(xndat, x_data);
//        for(x=0; x<xndat; x++) { // for all Groups within this sKill
//            if( x_data[x]->cnt!=0 ) continue; // this is an older superceding mechanism
//            t = 0;
//            o = x_data[x]->obs[t];
//            g = x_data[x]->g;
//            if(this->fitK[k] || this->fitG[g]) { // if fitK or fitG are raised
//                for(i=0; i<this->p->nS; i++) {
//                    // logit combination
//                    update =  getPI(x_data[x],i) * (1-getPI(x_data[x],i)) * x_data[x]->beta[t][i] * getB(x_data[x],i,o) / safe0num(x_data[x]->p_O_param);
//                    // Corbett combination - todo
//                }
//            }
//            // if fitK not raised - continue
//            if( !this->fitK[k] )
//                continue;
//            for(t=0; t<x_data[x]->ndat; t++) {
//                o = x_data[x]->obs[t];
//                // Gradient with respect to A
//                // \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
//                if( t>0 ) {
//                    for(i=0; i<this->p->nS /*&& param->fitparam[1]>0*/; i++)
//                        for(j=0; j<this->p->nS; j++) {
//                            this->gradA[k][i][j] -= x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
//                            this->gradAg[g][i][j] -= x_data[x]->beta[t][j] * this->B[k][j][o] * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
//                        }
//                }// if not first obs in sequence
//                // Gradient with respect to B
//                for(i=0; i<this->p->nS /*&& param->fitparam[2]>0*/; i++)
//                    this->gradB[k][i][o] -= x_data[g]->alpha[t][i] * x_data[g]->beta[t][i] / safe0num(x_data[g]->p_O_param * B[k][i][o]);
//            } // for all observations within skill-group
//        } // for all Groups
//        RecycleFitData(xndat, x_data, this->p);
//    }// for all sKills
//    
//} // computeGradients()

void HMMProblemAGK::setGradA (struct data* dt, FitBit *fb, NPAR kg_flag){
    NDAT t;
    NPAR o, i, j;
    NUMBER combined, deriv_logit;
    if(kg_flag == 0) { // k
        for(t=1; t<dt->ndat; t++) {
            o = dt->obs[t];
            for(i=0; i<this->p->nS /*&& fitparam[1]>0*/; i++)
                for(j=0; j<this->p->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = 1 / safe0num( this->A[ dt->k ][i][j] * (1-this->A[ dt->k ][i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * getB(dt,j,o) * dt->alpha[t-1][i] / safe0num(dt->p_O_param);
                }
        }
    }
    else
        for(t=1; t<dt->ndat; t++) {
            o = dt->obs[t];
            for(i=0; i<this->p->nS /*&& fitparam[1]>0*/; i++)
                for(j=0; j<this->p->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = 1 / safe0num( this->Ag[ dt->g ][i][j] * (1-this->Ag[ dt->g ][i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * getB(dt,j,o) * dt->alpha[t-1][i] / safe0num(dt->p_O_param);
                }
        }
}

//void HMMProblemAGK::computeGradients(NCAT xndat, struct data** x_data, FitBit *fb) {
//    NPAR nS = this->p->nS;
//    fb->toZero(FBS_GRAD);
//    
//	computeAlphaAndPOParam(xndat, x_data);
//	computeBeta(xndat, x_data);
//	NCAT x;
//	NDAT t;
//	NPAR i, j, o;
//    NUMBER combined, deriv_logit, update;
//    
//	for(x=0; x<xndat; x++) {
//		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
//		// Mask for gradPI is handled differently
//		// Gradient with respect to PI
//		t = 0;
//		o = x_data[x]->obs[t];
//		for(i=0; i<nS /*&& fitparam[0]>0*/; i++) {
//			fb->gradPI[i] -= x_data[x]->beta[t][i] * getB(x_data[x],i,o) / safe0num(x_data[x]->p_O_param);
//        }
//        
//		for(t=0; t<x_data[x]->ndat; t++) {
//			o = x_data[x]->obs[t];
//			// Gradient with respect to A
//			// \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
//			if( t>0 ) {
//				for(i=0; i<nS /*&& fitparam[1]>0*/; i++)
//					for(j=0; j<nS; j++) {
//                        combined = getA(x_data[x],i,j);
//                        deriv_logit = 1 / safe0num( this->A[ x_data[x]->k ][i][j] * (1-this->A[ x_data[x]->k ][i][j]) );
//                        update =  combined * (1-combined) * deriv_logit * x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
//                        fb->gradA[i][j] -= update;
////						fb->gradA[i][j] -= x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
//                    }
//				
//			}// if not first obs in sequence
//			// Gradient with respect to B
//			for(i=0; i<nS /*&& fitparam[2]>0*/; i++) {
//				fb->gradB[i][o] -= x_data[x]->alpha[t][i] * x_data[x]->beta[t][i] / safe0num(x_data[x]->p_O_param * getB(x_data[x],i,o));
//            }
//		} // for all observations within skill-group
//	} // for all groups in skill
//} // computeGradients()
//
//void HMMProblemAGK::computeGradientsG(NCAT xndat, struct data** x_data, FitBit *fb) {
//    NPAR nS = this->p->nS;
//    fb->toZero(FBS_GRAD);
//    
//	computeAlphaAndPOParam(xndat, x_data);
//	computeBeta(xndat, x_data);
//	NCAT x;
//	NDAT t;
//	NPAR i, j, o;
//    NUMBER combined, deriv_logit, update;
//    
//	for(x=0; x<xndat; x++) {
//		if( x_data[x]->cnt!=0 ) continue; // ... and the thing has not been computed yet (e.g. from group to skill)
//		// Mask for gradPI is handled differently
//		t = 0;
//		o = x_data[x]->obs[t];
//		for(t=0; t<x_data[x]->ndat; t++) {
//			o = x_data[x]->obs[t];
//			// Gradient with respect to A
//			// \frac{\partial J}{\partial a_{ij}} = - \frac{1}{L_{tot}} \sum_{t=2}^T \beta_t(j) b_j(o_t) \alpha_{t-1}(i)
//			if( t>0 ) {
//				for(i=0; i<nS; i++)
//					for(j=0; j<nS; j++) {
//                        combined = getA(x_data[x],i,j);
//                        deriv_logit = 1 / safe0num( this->Ag[ x_data[x]->g ][i][j] * (1-this->Ag[ x_data[x]->g ][i][j]) );
//                        update =  combined * (1-combined) * deriv_logit * x_data[x]->beta[t][j] * getB(x_data[x],j,o) * x_data[x]->alpha[t-1][i] / safe0num(x_data[x]->p_O_param);
//                        fb->gradA[i][j] -= update;
//                    }
//			}// if not first obs in sequence
//		} // for all observations within skill-group
//	} // for all groups in skill
//} // computeGradients()

void HMMProblemAGK::toFile(const char *filename) {
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
		fprintf(fid,"Ag\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%10.8f%s",this->Ag[g][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PI\t");
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

void HMMProblemAGK::fit() {
    NUMBER* loglik_rmse = init1DNumber(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure==STRUCTURE_Agk)
        loglik_rmse[0] += GradientDescent();
    else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemAGK::GradientDescent() {
	NCAT k, g;
    /*NPAR nS = this->p->nS, nO = this->p->nO;*/ NCAT nK = this->p->nK, nG = this->p->nG;
    NUMBER loglik = 0;
    FitResult fr;
    FitBit *fb = new FitBit(this->p);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    if(this->p->solver==METHOD_CGD) {
        fb->init(FBS_GRADm1);
        fb->init(FBS_DIRm1);
    }
	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill==1) {
        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(0/*use skill 0*/, this->p->ndata, this->p->k_data, 0/* by skill*/, fb, true /*is1SkillForAll*/);
        printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
    }
	
	//
	// Main fit
	//
    int first_iteration_qualify = 0; // at what iteration, qualification for skill/group convergence should start
    int iterations_to_qualify = 2; // how many concecutive iterations necessary for skill/group to qualify as converged
    NPAR* iter_qual_skill = Calloc(NPAR, nK);
    NPAR* iter_qual_group = Calloc(NPAR, nG);
    int skip_k = 0, skip_g = 0;
    
    int i = 0; // count runs
    while(skip_k<nK || skip_g<nG) {
        //
        // Skills first
        //
        for(k=0; k<nK && skip_k<nK; k++) { // for all A,B-by-skill
            if(iter_qual_skill[k]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->k_numg[k];
            struct data** x_data = this->p->k_g_data[k];
            // link and fit
            fb->linkPar( this->getPI(k), this->getA(k), this->getB(k));// link skill 0 (we'll copy fit parameters to others
            fr = GradientDescentBit(k/*use skill x*/, xndat, x_data, 0/*skill*/, fb, false /*is1SkillForAll*/);
            // decide on convergence
            if(i>=first_iteration_qualify) {
                if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==nG) { // converged quick, or don't care (others all converged
                    iter_qual_skill[k]++;
                    if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
                        if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                        skip_k++;
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                        }
                    }
                }
                else
                    iter_qual_skill[k]=0;
            } // decide on convergence
        } // for all skills
        //
        // PIg second
        //
        for(g=0; g<nG && skip_g<nG; g++) { // for all PI-by-user
            if(iter_qual_group[g]==iterations_to_qualify)
                continue;
            NCAT xndat = this->p->g_numk[g];
            struct data** x_data = this->p->g_k_data[g];
            // vvvvvvvvvvvvvvvvvvvvv ONLY PART THAT IS DIFFERENT FROM HMMProblemPiGK
            fb->linkPar(NULL, this->getAg(g), NULL);
            // ^^^^^^^^^^^^^^^^^^^^^
            // decide on convergence
            fr = GradientDescentBit(g/*use skill x*/, xndat, x_data, 1/*group*/, fb, false /*is1SkillForAll*/);
            if(i>=first_iteration_qualify) {
                if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==nK) { // converged quick, or don't care (others all converged
                    iter_qual_group[g]++;
                    if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                        if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                        skip_g++;
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                            computeAlphaAndPOParam(xndat, x_data);
                            printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                        }
                    }
                }
                else
                    iter_qual_group[g]=0;
            } // decide on convergence
        } // for all groups
        i++;
    }
    delete fb;
    // compute loglik
    fr.pO = 0.0;
    for(k=0; k<nK; k++) { // for all A,B-by-skill
        fr.pO = getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
        loglik +=fr.pO*(fr.pO>0);
    }
    return loglik;
}

//NUMBER HMMProblemAGK::doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
//                                          NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB) {
//	NPAR i,j,m;
//    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
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
//NUMBER HMMProblemAGK::doLinearStepPLoGroup(NCAT xndat, struct data** x_data,  NUMBER **a_A, NUMBER **a_gradA) {
//	NPAR i, j;
//	// first scale down gradients
//	doLog10Scale2DGentle(a_gradA,  a_A,  this->p->nS, this->p->nS);
//	
//	NUMBER *PI_cpy = init1DNumber(this->p->nS); // safe copy
//	NUMBER **A_cpy = init2DNumber(this->p->nS, this->p->nS); // safe copy
//	NUMBER e = this->p->ArmijoSeed; // step seed
//	bool compliesArmijo = false;
//	NUMBER f_xk = HMMProblem::getSumLogPOPara(xndat, x_data);
//	NUMBER f_xkplus1;
//	
//	cpy2DNumber(a_A,  A_cpy,  this->p->nS, this->p->nS); // save copy
//	// compute p_k * -p_k
//	NUMBER p_k_by_neg_p_k = 0;
//	for(i=0; i<this->p->nS; i++) {
//		for(j=0; j<this->p->nS; j++) p_k_by_neg_p_k -= a_gradA[i][j]*a_gradA[i][j];
//    }
//	int iter = 0; // limit iter steps to 20
//	while( !compliesArmijo && e > this->p->ArmijoMinStep) {
//		// update
//		for(i=0; i<this->p->nS; i++) {
//			for(j=0; j<this->p->nS; j++)
//				a_A[i][j] = A_cpy[i][j] - e * a_gradA[i][j];
//		}
//		// scale
//		if( !this->hasNon01Constraints() ) {
//			for(i=0; i<this->p->nS; i++)
//				projectsimplex(a_A[i], this->p->nS);
//		} else {
//			for(i=0; i<this->p->nS; i++)
//				projectsimplexbounded(a_A[i], this->getLbA()[i], this->getUbA()[i], this->p->nS);
//		}
//		// recompute alpha and p(O|param)
//        //		zeroLabels(xndat, x_data); // THIS IS NOT DONE HERE
//        computeAlphaAndPOParam(xndat, x_data);
//		// compute f(x_{k+1})
//		f_xkplus1 = HMMProblem::getSumLogPOPara(xndat, x_data);
//		// compute Armijo compliance
//		compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
//		e /= (compliesArmijo)?1:this->p->ArmijoReduceFactor;
//		iter++;
//	} // armijo loop
//    if(!compliesArmijo) {
//        e = 0;
//        cpy2DNumber(A_cpy,  a_A,  this->p->nS, this->p->nS); // save copy
//    }
//	free(PI_cpy);
//	free2DNumber(A_cpy, this->p->nS);
//    RecycleFitData(xndat, x_data, this->p);
//    return e;
//} // doLinearStep
