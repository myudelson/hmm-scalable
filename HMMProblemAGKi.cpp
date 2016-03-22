//
//  HMMProblemAGKi.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
//

#include "HMMProblemAGKi.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemAGKi::HMMProblemAGKi(struct param *param) : HMMProblemAGK(param) {
}

HMMProblemAGKi::~HMMProblemAGKi() {
    destroy();
}

void HMMProblemAGKi::init(struct param *param) {
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


void HMMProblemAGKi::destroy() {
	if( this->Ag != NULL) {
        free3D<NUMBER>(this->Ag, this->p->nG, this->p->nS);
        this->Ag = NULL;
    }
}// ~HMMProblemAGKi


void HMMProblemAGKi::setGradA (struct data* dt, FitBit *fb, NPAR kg_flag){
    NDAT t;
    NPAR o, i, j;
    NUMBER combined, deriv_logit;
    if(kg_flag == 0) { // k
        for(t=1; t<dt->n; t++) {
//            o = dt->obs[t];
            o = this->p->dat_obs->get( dt->ix[t] );
            for(i=0; i<this->p->nS /*&& fitparam[1]>0*/; i++)
                for(j=0; j<this->p->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = (1 + logit(this->Ag[ dt->g ][i][j])/*ADDITION*/) / safe0num( this->A[ dt->k ][i][j] * (1-this->A[ dt->k ][i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * ((o<0)?1:getB(dt,j,o)) * dt->alpha[t-1][i] /
                        safe0num(dt->p_O_param) + L2penalty(this->p,this->A[ dt->k ][i][j], 0.5); // PENALTY
                }
        }
    }
    else
        for(t=1; t<dt->n; t++) {
//            o = dt->obs[t];
            o = this->p->dat_obs->get( dt->ix[t] );
            for(i=0; i<this->p->nS /*&& fitparam[1]>0*/; i++)
                for(j=0; j<this->p->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = ( 1 + logit(this->A[ dt->k ][i][j])/*ADDITION*/) / safe0num( this->Ag[ dt->g ][i][j] * (1-this->Ag[ dt->g ][i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * ((o<0)?1:getB(dt,j,o)) * dt->alpha[t-1][i]
                        / safe0num(dt->p_O_param) + L2penalty(this->p,this->Ag[ dt->g ][i][j], 0.5); // PENALTY
                }
        }
}

void HMMProblemAGKi::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure==STRUCTURE_Agki)
        loglik_rmse[0] += GradientDescent();
    else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NUMBER HMMProblemAGKi::GradientDescent() {
	NCAT k, g;
    /*NPAR nS = this->p->nS, nO = this->p->nO;*/ NCAT nK = this->p->nK, nG = this->p->nG;
    NUMBER loglik = 0;
    FitResult fr;
    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    if(this->p->solver==METHOD_CGD) {
        fb->init(FBS_GRADm1);
        fb->init(FBS_DIRm1);
    }
	//
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill>0) {
        fb->linkPar( this->getPI(0), this->getA(0), this->getB(0));// link skill 0 (we'll copy fit parameters to others
        fr = GradientDescentBit(0/*use skill 0*/, this->p->nSeq, this->p->k_data, 0/* by skill*/, fb, true /*is1SkillForAll*/);
        printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
    }
	
	//
	// Main fit
	//
    if( this->p->single_skill!=2 ) { // if not "force single skill"
        int first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
        int iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
        NPAR* iter_qual_skill = Calloc(NPAR, (size_t)nK);
        NPAR* iter_qual_group = Calloc(NPAR, (size_t)nG);
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
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_qual_group != NULL) free(iter_qual_group);
    } // if not "force single skill"
    
    delete fb;
    // compute loglik
    fr.pO = 0.0;
    for(k=0; k<nK; k++) { // for all A,B-by-skill
        fr.pO = getSumLogPOPara(this->p->k_numg[k], this->p->k_g_data[k]);
        loglik +=fr.pO*(fr.pO>0);
    }
    return loglik;
}

