/*
 
 Copyright (c) 2012-2017, Michael (Mikhail) Yudelson
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
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utilsSt.h"
#include "FitBitSt.h"
#include <math.h>
#include "HMMProblemSt.h"
#include <map>

HMMProblemSt::HMMProblemSt() {
}

HMMProblemSt::HMMProblemSt(struct task *task) {
//    printf("HMMProblemSt::HMMProblemSt\n");
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    this->model_param_n = task->nK * 4;
    init(task);
}// HMMProblemSt

void HMMProblemSt::init(struct task *task) {
//    printf("HMMProblemSt::init\n");
   // init
    this->task = task;
    this->non01constraints = true;
    this->null_obs_ratios = Calloc(NUMBER, (size_t)this->task->nO);
//    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    if( this->task->solver == METHOD_CGD && this->task->solver_setting == -1)
        this->task->solver_setting = 1; // default Fletcher-Reeves
    
    NPAR nS = this->task->nS, nO = this->task->nO;
    NCAT nK = this->task->nK, nG = this->task->nG;
    
    //
    // setup params
    //
    this->skill1_n = nS * (1 + nS + nO); // number of params per 1 skill slot
    this->param_skill = init1D<NUMBER>( (NDAT)(this->skill1_n*nK) );
    this->gradient_skill = init1D<NUMBER>( (NDAT)(this->skill1_n*nK) );
    this->active_set_block = initToValue1D<NPAR>((NDAT)nK, 1);
    this->unblocked_simplex = initToValue1D<NPAR>((NDAT)(nK * (1 + 2 * nS)), 1);
    
    NPAR i, j, m;
    int offset, idx;

    // write default parameters first
    NUMBER sum = 0;
    
    // vv populate skill parameters
    for(NCAT k=0; k<nK; k++) {
        // vv populate PI
        for(i=0; i<((nS)-1); i++) {
            this->param_skill[ k*this->skill1_n + i ] = this->task->init_param_values[i];
            sum  += this->task->init_param_values[i];
        } // ^^ populate PI
        this->param_skill[ k*this->skill1_n + nS - 1 ] = 1 - sum;
        // vv populate A
        offset = (int)(nS-1);
        for(i=0; i<nS; i++) {
            sum = 0;
            for(j=0; j<((nS)-1); j++) {
                idx = (int)(offset + i*((nS)-1) + j);
                this->param_skill[ k*this->skill1_n + nS + i*nS + j ] = this->task->init_param_values[idx];
                sum  += this->task->init_param_values[idx];
            }
            this->param_skill[ k*this->skill1_n + nS + i*nS + nS-1 ] = 1 - sum;
        } // ^^ populate A
        // vv populate B
        offset = (int)((nS-1) + nS*(nS-1));
        for(i=0; i<nS; i++) {
            sum = 0;
            for(m=0; m<((nO)-1); m++) {
                idx = (int)(offset + i*((nO)-1) + m);
                this->param_skill[ k*this->skill1_n + nS + nS*nS + i*nO + m ] = this->task->init_param_values[idx];
                sum += this->task->init_param_values[idx];
            }
            this->param_skill[ k*this->skill1_n + nS + nS*nS + i*nO + nO - 1 ] = 1 - sum;
        } // ^^ populate B
    } // ^^ populate skill parameters
    
    // check only first mass-produced skill
    if( this->task->do_not_check_constraints==0 && !checkSkillConstraints(this->param_skill)) {
        fprintf(stderr,"params do not meet constraints.\n");
        exit(1);
    }

    if(task->initfile[0]!=0) { // from setup file
        // if needs be -- read in init params from a file
        this->readModel(task->initfile, false /* read and upload but not overwrite*/);
    }
    
    // init and populate boundaries
    this->lb_param_skill = init1D<NUMBER>(this->skill1_n);
    cpy1D(this->task->param_values_lb, this->lb_param_skill, this->skill1_n);
    this->ub_param_skill = init1D<NUMBER>(this->skill1_n);
    cpy1D(this->task->param_values_ub, this->ub_param_skill, this->skill1_n);
    
    // init supportive data structures
    NDAT Nst = this->task->Nst;
    this->alpha = init2D<NUMBER>(Nst, nS); // forward variables, 2D array of Nst rows nS values in each ((1..nS)){Nst}
    this->po_param_gk = init2D<NUMBER>(nG, nK); //
    this->loglik_k = init1D<NUMBER>(nK); //
    this->ndat_k = init1D<NDAT>(nK);
    this->beta = init2D<NUMBER>(Nst, nS); // backward variables, 2D array of Nst rows nS values in each ((1..nS)){Nst}
    this->scale = (this->task->is_scaled==1) ? initToValue1D<NUMBER>(Nst, 0) : NULL;
    this->gamma = NULL;
    this->xi = NULL;
    if(this->task->solver == METHOD_BW) {
        this->gamma = init2D<NUMBER>(Nst, nS);
        this->xi = init3D<NUMBER>(Nst, nS, nS);
    }
    this->backward_ix = init1D<NDAT>(Nst); // Nst sized array of pointing to previous student-skill repetition, -1 means none before
    this->forward_ix = init1D<NDAT>(Nst); // Nst sized array of pointing to next student-skill repetition, -1 means no more
    this->is_fwd_bwd_built=false; // flag for noting whether backward_ix & forward_ix are already built
    
    this->group_skill_map = init3D<NUMBER>(nG, nK, nS);//UNBOOST

}// HMMProblemSt::init

HMMProblemSt::~HMMProblemSt() {
//    printf("HMMProblemSt::~HMMProblemSt\n");
    NDAT Nst = this->task->Nst;
//    printf("before null_obs_ratios\n");
    if(this->null_obs_ratios != NULL) free(this->null_obs_ratios);
//    printf("before param_skill\n");
    if(this->param_skill != NULL) free(this->param_skill);
//    printf("before gradient_skill\n");
    if(this->gradient_skill != NULL) free(this->gradient_skill);
//    printf("before lb_param_skill\n");
    if(this->lb_param_skill != NULL) free(this->lb_param_skill);
//    printf("before ub_param_skill\n");
    if(this->ub_param_skill != NULL) free(this->ub_param_skill);
//    printf("before alpha %d\n",this->alpha != NULL);
    if(this->alpha != NULL) free2D(this->alpha, Nst);
//    printf("before po_param_gk\n");
    if(this->po_param_gk != NULL) free2D(this->po_param_gk, this->task->nG);
//    printf("before loglik_k\n");
    if(this->loglik_k != NULL) free(this->loglik_k);
//    printf("before ndat_k\n");
    if(this->ndat_k != NULL) free(this->ndat_k);
//    printf("before beta\n");
    if(this->beta != NULL) free2D(this->beta, Nst);
//    printf("before scale\n");
    if(this->scale != NULL) free(this->scale);
//    printf("before gamma\n");
    if(this->gamma != NULL) free2D(this->gamma, Nst);
//    printf("before xi\n");
    if(this->xi != NULL) free3D(this->xi, Nst, this->task->nS);
//    printf("before backward_ix\n");
    if(this->backward_ix != NULL) free(this->backward_ix);
//    printf("before forward_ix\n");
    if(this->forward_ix != NULL) free(this->forward_ix);
    
    free3D(this->group_skill_map, this->task->nG, this->task->nK); //UNBOOST
    
//    printf("before active_set_block\n");
    free(this->active_set_block);
//    printf("before unblocked_simplex\n");
    free(this->unblocked_simplex);
    
}// ~HMMProblemSt

NDAT HMMProblemSt::PI(NCAT k, NPAR i) {
    return this->skill1_n*k + i;
}
NDAT HMMProblemSt::A(NCAT k, NPAR i, NPAR j) {
    NPAR nS = this->task->nS;
    return this->skill1_n*k + nS + i * nS + j;
}
NDAT HMMProblemSt::B(NCAT k, NPAR i, NPAR o){
    NPAR nS = this->task->nS, nO = this->task->nO;
    return this->skill1_n*k + (1+nS)*nS + i * nO + o;
}

NUMBER* HMMProblemSt::getPI(NCAT x) {
	if( x > (this->task->nK-1) ) {
		fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, this->task->nK-1);
		exit(1);
	}
	return &this->param_skill[this->skill1_n*x + 0];
}

NUMBER* HMMProblemSt::getA(NCAT x) {
    NPAR nS = this->task->nS;
    NCAT nK = this->task->nK;
	if( x > (nK-1) ) {
        fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, nK-1);
		exit(1);
	}
	return &this->param_skill[this->skill1_n*x + nS];
}

NUMBER* HMMProblemSt::getB(NCAT x) {
    NPAR nS = this->task->nS;
    NCAT nK = this->task->nK;
	if( x > (nK-1) ) {
        fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, nK-1);
		exit(1);
	}
	return &this->param_skill[ this->skill1_n*x + nS*(nS+1) ];
}

NUMBER HMMProblemSt::getPI(NCAT x, NPAR i) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    NCAT nK = this->task->nK;
    if( x > (this->task->nK-1) ) {
        fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, nK-1);
        exit(1);
    }
    // checks done
    return this->param_skill[this->skill1_n*x + i];
}

NUMBER HMMProblemSt::getA(NCAT x, NPAR i, NPAR j) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    NPAR nS = this->task->nS;
    NCAT nK = this->task->nK;
    if( x > (nK-1) ) {
        fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, nK-1);
        exit(1);
    }
    // checks done
    return this->param_skill[this->skill1_n*x + nS + i*nS + j];
}

NUMBER HMMProblemSt::getB(NCAT x, NPAR i, NPAR m) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    NPAR nS = this->task->nS;
    NPAR nO = this->task->nO;
    NCAT nK = this->task->nK;
    if( x > (nK-1) ) {
        fprintf(stderr,"While accessing PI, skill index %d exceeded the last skill index %d.\n", x, nK-1);
        exit(1);
    }
    if(m<0)
        return 1;
    // checks done
    return this->param_skill[this->skill1_n*x + nS * (1 + nS )+ i*nO + m];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemSt::getPI(struct context* ctx, NPAR i) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    // checks done
    return this->param_skill[this->skill1_n*ctx->k + i];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemSt::getA (struct context* ctx, NPAR i, NPAR j) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    // checks done
    NPAR nS = this->task->nS;
    return this->param_skill[this->skill1_n*ctx->k + nS + i*nS + j];
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemSt::getB (struct context* ctx, NPAR i, NPAR m) {
//    if(task->structure != STRUCTURE_SKILL) {
//        fprintf(stderr,"Model structure specified is not supported and should have been caught earlier\n");
//        exit(1);
//    }
    // checks done
    NPAR nS = this->task->nS;
    NPAR nO = this->task->nO;
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;
    return this->param_skill[this->skill1_n*ctx->k + nS * (1 + nS )+ i*nO + m];
}

bool HMMProblemSt::checkSkillConstraints(NUMBER* param_skill) {
    NPAR nS = this->task->nS, nO = this->task->nO;
	NCAT skill1_n = nS * (1 + nS + nO); // number of params per 1 skill slot
    NUMBER sum = 0.0;
    
    for(NPAR l=0; l<skill1_n; l++) {
        if( param_skill[l]>1.0 || param_skill[l]<0.0)
            return false;
        sum += param_skill[l];
    }
    if( sum/(skill1_n/2) != 1.0 )
        return false;
    return true;
}

void HMMProblemSt::initAlphaEtAl() {
    NPAR nS = this->task->nS;
    NCAT nK = this->task->nK, nG = this->task->nG;
    NDAT Nst = this->task->Nst;
    NPAR f_is_scaled = this->task->is_scaled==1;
    toZero2D(this->alpha, Nst, nS);
    // if scaled, new parameter values are multiplied, if not scaled – added
    if(f_is_scaled)
        toValue2D(this->po_param_gk, nG, nK, 1.0);
    else
        toZero2D(this->po_param_gk, nG, nK);
    toZero1D(this->loglik_k, nK);
    
    toValue3D(this->group_skill_map, (NDAT)nG, (NDAT)nK, (NDAT)nS, -1.0); // -1 means not set (use pLo)
}

NUMBER HMMProblemSt::computeAlphaAndPOParam(NUMBER *metrics_res) {
    NPAR nS = this->task->nS, nO = this->task->nO, o = -1, isTarget;
    NDAT Nst = this->task->Nst, N = this->task->N;
    NCAT g,k, nK = this->task->nK, nG = this->task->nG;
    NPAR f_metrics_target_obs = this->task->metrics_target_obs;
    NPAR f_is_scaled = this->task->is_scaled==1;
    // prepare values
    NUMBER loglik = 0.0, loglik_nonull = 0.0, sse = 0.0, sse_nonull = 0.0, ncorr = 0.0, ncorr_nonull = 0.0;
    NUMBER p = 0.0;///, corr = 0.0;
    
    // initialize reused parameters before run
//    NUMBER ***group_skill_map; // traced state values for student (group) * skill //UNBOOST
//    group_skill_map = init3D<NUMBER>(nG, nK, nS);//UNBOOST
//    toValue3D(group_skill_map, (NDAT)nG, (NDAT)nK, (NDAT)nS, -1.0); // -1 means not set (use pLo)
    initAlphaEtAl();

    NDAT **alpha_gk_ix = NULL;
    
    NUMBER **predict = this->task->dat_predict; // running value of prediction, do we save it?
    NUMBER *predict_k = this->task->dat_predict_k; // running value of prediction, do we save it?
    if(metrics_res!=NULL & this->task->predictions==2) {
        if(predict_k==NULL) predict_k = initToValue1D<NUMBER>(this->task->Nst,0.0);
        else toZero1D<NUMBER>(predict_k, this->task->Nst);
    }
    // toZero2D(predict,N,nO); // don't blast everything, the predictions are inited to 0 line by line as the overall data is updated according to the active set

    struct context* ctx = new context;

    if(!this->is_fwd_bwd_built) { // if helper structures were not built
        toValue1D(this->forward_ix, Nst, -1);
        toValue1D(this->backward_ix, Nst, -1);
        toValue1D(this->ndat_k, nK, 0);
        alpha_gk_ix = initToValue2D<NDAT>(nG, nK, -1);
    }
    
    for(NDAT t=0; t<N; t++) {
        o = this->task->dat_obs[t]; // observation y in the data is 1-right, 0-wrong; math assumes HMM-Scalable 1-right, 2-wrong, with -1 taken out, so 0-right, 1-wrong
        isTarget = this->task->metrics_target_obs == o;
        //corr = 1-o; // corr 1 - right, 0 - wrong
        g = this->task->dat_group[t]; // -1, because in data they were 1-starting
        ctx->g=g;
        ctx->o=o;
        ctx->t=t;
        // grab skill array (if exists)
        NCAT *ar;
        NPAR n, i;
        getSkillsAtRow(this->task, t, &ar, &n);
        // deal with null skill
        if(ar[0]<0 && metrics_res!=NULL) { // account for no skill label only if we need metrics (likely called from predict)
            isTarget = this->null_skill_obs==o;
            // old: rmse += pow(isTarget - hmm->null_skill_obs_prob,2);
            sse += pow(isTarget - this->null_obs_ratios[f_metrics_target_obs],2);
            // old: accuracy += isTarget == (hmm->null_skill_obs_prob>=0.5);
            ncorr += isTarget == (this->null_obs_ratios[f_metrics_target_obs]==maxn(this->null_obs_ratios,nO) && this->null_obs_ratios[f_metrics_target_obs] > 1/nO);
            loglik -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1 - this->null_skill_obs_prob);
            
            for(NPAR m=0; m<nO; m++) {
                predict[t][m] = this->null_obs_ratios[m];
            }
        }
        if(ar[0]<0) {
            continue;
        }
        
        // produce corrects first
        producePCorrect(group_skill_map, ar, n, predict[t], ctx); //UNBOOST
        
        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);
        
        for(int l=0; l<n && n_active==n/**/; l++) {
            k = ar[l];
            ctx->k=k;
            NDAT tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array
            if(f_is_scaled) this->scale[tt] = 0;
            
            NDAT pre_tt = (this->is_fwd_bwd_built) ? this->backward_ix[tt]: alpha_gk_ix[g][k]; // previous tt index of forward variable alpha, alpha_gk_ix[g][k] is running (prior) tt index of this group/skill (gk)
            
            if(pre_tt==-1) { // it's alpha(1,i)
                // compute \alpha_1(i) = \pi_i b_i(o_1)
                for(i=0; i<nS; i++) {
                    this->alpha[tt][i] = getPI(ctx,i) * ((o<0)?1:getB(ctx,i,o)); // if observatiob unknown use 1
                    if(f_is_scaled) {
                        this->scale[tt] += this->alpha[tt][i];
                    }
                }
            } else { // it's alpha(t,i)
                // compute \alpha_{t}(i) = b_j(o_{t})\sum_{j=1}^N{\alpha_{t-1}(j) a_{ji}}
                for(NPAR i=0; i<nS; i++) {
                    for(NPAR j=0; j<nS; j++) {
                        alpha[tt][i] += alpha[pre_tt][j] * getA(ctx,j,i);
                    }
                    this->alpha[tt][i] *= ((o<0)?1:getB(ctx,i,o)); // if observatiob unknown use 1
//                    if( this->alpha[t][i] < 0 || this->alpha[t][i] > 1)
//                        fprintf(stderr, "ERROR! alpha value is not within [0, 1] range!\n");
                    if(f_is_scaled) this->scale[tt] += alpha[tt][i];
                }
            }
            // update backward_ix, forward_ix if necessary
            if(!this->is_fwd_bwd_built) {
                if( pre_tt != -1 ) { // there exists a prior datapoint for this student and skill
                  this->forward_ix[pre_tt] = tt; // point from that previous position to current tt index
                  this->backward_ix[tt] = pre_tt; // point from current position to the previous
                }
                // now update the last tt index
                alpha_gk_ix[g][k] = tt;
                // count datapoints per k
                this->ndat_k[k]++;
            }// update backward_ix, forward_ix if necessary
        } // all skills in a row
        
        // update per row values, skill, group, Elo, etc
        updateValuesLocal(group_skill_map, ar, n, predict[t], ctx); //UNBOOST
        
        // update obj val
        p = safe01num( predict[t][f_metrics_target_obs]);
        loglik -= safelog(  p)*   isTarget  +  safelog(1-p)*(1-isTarget);
        if(metrics_res!=NULL) { // we are predicting
            loglik_nonull -= safelog(  p)*   isTarget  +  safelog(1-p)*(1-isTarget);
            sse += pow(isTarget-predict[t][f_metrics_target_obs],2);
            sse_nonull += pow(isTarget-predict[t][f_metrics_target_obs],2);
            ncorr += isTarget == (predict[t][f_metrics_target_obs]==maxn(predict[t],nO) && predict[t][f_metrics_target_obs]>1/nO);
            ncorr_nonull += isTarget == (predict[t][f_metrics_target_obs]==maxn(predict[t],nO) && predict[t][f_metrics_target_obs]>1/nO);
            if(this->task->predictions==2) {
                for(int l=0; l<n; l++) {
                    k = ar[l];
                    NDAT tt = l + this->task->dat_skill_rix[t];
                    predict_k[tt] = group_skill_map[g][k][0];
                }
            }
        }



    } // all rows in the data
    
    // compute this->po_param_gk and this->loglik_k
    for(int t=0; t<N; t++) {
        g = this->task->dat_group[t]; // -1, because in data they were 1-starting
        NCAT *ar;
        NPAR n, i;
        getSkillsAtRow(this->task, t, &ar, &n);
        
        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);
        
        for(int l=0; l<n && (n_active==n)/**/; l++) {
            k = ar[l];
            ctx->k=k;
            NDAT tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array
            if(f_is_scaled) // if scaled – multiply all scales
                this->po_param_gk[g][k] *= this->scale[tt];
            else { // if not scaled, add alpha's for the last row in g/k sequence
                if(this->forward_ix[tt]==-1) {
                    for(i=0; i<nS; i++) this->po_param_gk[g][k] += this->alpha[tt][i];
                    this->loglik_k[k] -= safelog(this->po_param_gk[g][k]);
                }
            }
        }
    } // compute this->po_param_gk & compute this->loglik_k
    // compute this->loglik_k

    // recycle
    this->is_fwd_bwd_built = true; // set in any way
    if(alpha_gk_ix!=NULL) free2D(alpha_gk_ix, nG);
    if(metrics_res != NULL) {
        metrics_res[0] = loglik;
        metrics_res[1] = loglik_nonull;
        metrics_res[2] = sse;
        metrics_res[3] = sse_nonull;
        metrics_res[4] = ncorr;
        metrics_res[5] = ncorr_nonull;
        if(this->task->predictions==2) {
            this->task->dat_predict_k = predict_k;
        }
    }
    
//    free3D(group_skill_map, this->task->nG, this->task->nK); //UNBOOST
    delete ctx;
    return loglik; //TODO, figure out a diff way to sum it, and not multiple times
}

void HMMProblemSt::computeBeta() {
    NPAR nS = this->task->nS, o = -1;//, isTarget;
    NDAT Nst = this->task->Nst, N = this->task->N;
    NCAT g,k;
    NPAR f_is_scaled = this->task->is_scaled==1;
    struct context* ctx = new context;
    toZero2D(this->beta, Nst, nS);

    // Backwards pass
    for(int t=(N-1); t>=0; t--) {
        if(t==96) {
            int a =0;
        }
        o = this->task->dat_obs[t]; // observation y in the data is 1-right, 0-wrong; math assumes HMM-Scalable 1-right, 2-wrong, with -1 taken out, so 0-right, 1-wrong
        // isTarget = this->task->metrics_target_obs == o;
        g = this->task->dat_group[t]; // -1, because in data they were 1-starting
        ctx->g=g;
        ctx->o=o;
        ctx->t=t;
        // grab skill array (if exists)
        NCAT *ar;
        NPAR n, i, j;
        getSkillsAtRow(this->task, t, &ar, &n);
        // deal with null skill
        if(ar[0]<0) { // if no skill label
            // fprintf(stderr, "WARNING! We are not dealing with skill-less observation\n");
            continue;
        }
        
        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);

        for(NPAR l=0; l<n && (n_active==n)/**/; l++) {
            k = ar[l];
            ctx->k=k;
            NDAT tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array

            NDAT fwd_ix = this->forward_ix[tt];
            if(fwd_ix==-1) { // this is the end of student,skill sequence, no forward index (-1)
                for(i=0;i<nS;i++)
                    beta[tt][i] = (f_is_scaled)?this->scale[tt]:1.0;
            } else { // not the end of student,skill sequence
                // o here is o from next step in oo_tt
                NPAR o_next = this->task->dat_obs_stacked[fwd_ix];
                for(i=0; i<nS; i++) {
                    for(j=0; j<nS; j++)
                        this->beta[tt][i] += this->beta[fwd_ix][j] * getA(ctx,i,j) * ((o_next<0)?1:getB(ctx,j,o_next)); // if observatiob unknown use 1
                    if(f_is_scaled) this->beta[tt][i] *= this->scale[tt];
                }
            }
        } // all skills in a row
        // tt_off -= fd->skill_n[t]; // tt offset
    } // all rows in the data
    delete ctx;
}

void HMMProblemSt::computeXiGamma(){
    NPAR nS = this->task->nS, o_tp1 = -1, i, j;
    NDAT Nst = this->task->Nst, N = this->task->N;
    NCAT g,k;
    // prepare values
    toZero2D(this->gamma, Nst, nS);
    toZero3D(this->xi, Nst, nS, nS);
    struct context* ctx = new context;
    NUMBER denom = 0.0;

    for(NDAT t=0; t<N; t++) {
        g = this->task->dat_group[t]; 
        ctx->g=g;
        ctx->t=t;
        // grab skill array (if exists)
        NCAT *ar;
        NPAR n;
        getSkillsAtRow(this->task, t, &ar, &n);
        // deal with null skill
        if(ar[0]<0) { // if no skill label
            // fprintf(stderr, "WARNING! We are not dealing with skill-less observation\n");
            continue;
        }
        
        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);

        for(int l=0; l<n && (n_active==n)/**/; l++) {
            k = ar[l];
            ctx->k=k;
            NDAT tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array
            NDAT next_tt = this->forward_ix[tt];
            if(next_tt==-1) continue;
            o_tp1 = this->task->dat_obs_stacked[next_tt];
            
            denom = 0.0;
            for(i=0; i<nS; i++) {
                for(j=0; j<nS; j++) {
                    denom += this->alpha[tt][i] * getA(ctx,i,j) * this->beta[next_tt][j] * ((o_tp1<0)?1:getB(ctx,j,o_tp1));
                }
            }
            for(i=0; i<nS; i++) {
                for(j=0; j<nS; j++) {
                    this->xi[tt][i][j] = this->alpha[tt][i] * getA(ctx,i,j) * this->beta[next_tt][j] * ((o_tp1<0)?1:getB(ctx,j,o_tp1)) / ((denom>0)?denom:1); //
                    this->gamma[tt][i] += this->xi[tt][i][j];
                }
            }
        } // all row's skills
    } // all N rows
    delete ctx;
}

NUMBER HMMProblemSt::computeGradients() {
    NPAR nS = this->task->nS, o = -1;//, isTarget;
    NDAT Nst = this->task->Nst, N = this->task->N, c = 0;
    NCAT g,k, nK = this->task->nK;
    // prepare values
    toZero2D(this->alpha, Nst, nS);
    struct context* ctx = new context;
    
    toZero1D(this->gradient_skill, this->skill1_n*nK);
    
    NUMBER loglik = computeAlphaAndPOParam(NULL); // just alpha and p(O|param)
    computeBeta();
    // this->toFileAlphaBetaObs("trainhmmst_alphabetaobs0.txt", 0); // printing out alpha and beta arrays
    
    // main gradient computation
    for(int t=0; t<N; t++) {
        o = this->task->dat_obs[t]; // observation y in the data is 1-right, 0-wrong; math assumes HMM-Scalable 1-right, 2-wrong, with -1 taken out, so 0-right, 1-wrong
        //isTarget = this->task->metrics_target_obs == o;
        g = this->task->dat_group[t]; // -1, because in data they were 1-starting
        ctx->g=g;
        ctx->o=o;
        ctx->t=t;
        // grab skill array (if exists)
        NCAT *ar;
        NPAR n, i, j;
        getSkillsAtRow(this->task, t, &ar, &n);
        // deal with null skill
        if(ar[0]<0) { // if no skill label
//            fprintf(stderr, "WARNING! We are not dealing with skill-less observation\n");
            continue;
        }
        
        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);

        for(int l=0; l<n && (n_active==n)/**/; l++) {
            k = ar[l];
            ctx->k=k;
            NDAT tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array
            NDAT tt_m1 = this->backward_ix[tt];

            // grad PI
            if (tt_m1==-1) {
                c += (k==0);
                for(i=0; i<nS; i++) {
                    this->gradient_skill[PI(k,i)] -= this->beta[tt][i] * ((o<0)?1:this->param_skill[B(k,i,o)]) / safe0num(this->po_param_gk[g][k]);
                    //fb->gradPI[i] -= dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
                }
            } else {
            // grad A
                for(i=0; i<nS; i++)
                    for(j=0; j<nS; j++)
                        this->gradient_skill[A(k,i,j)] -= this->beta[tt][j] * ((o<0)?1:this->param_skill[B(k,j,o)]) * this->alpha[tt_m1][i] / safe0num(this->po_param_gk[g][k]);
            }
            // grad B
            for(i=0; i<nS; i++)
                this->gradient_skill[B(k,i,o)] -= this->alpha[tt][i] * this->beta[tt][i] / safe0num(this->po_param_gk[g][k] * this->param_skill[B(k,i,o)]);
        }// for all skills in a row
    }// for all rows
//    if( this->task->Cslices>0 ) { // penalty
//        fb->addL2Penalty(FBV_PI, this->p, (NUMBER)ndat);
//        fb->addL2Penalty(FBV_A, this->p, (NUMBER)ndat);
//        fb->addL2Penalty(FBV_B, this->p, (NUMBER)ndat);
//    }
    
//    NUMBER sum = 0;
//    for(int l=0;l<(this->skill1_n*nK);l++) sum += this->gradient_skill[l];
    delete ctx;
    
    return loglik;
} // computeGradients()

void HMMProblemSt::toFile(const char *filename) {
//    switch(this->task->structure)
//    {
//        case STRUCTURE_SKILL:
            toFileSkill(filename);
//            break;
//        case STRUCTURE_GROUP:
//            toFileGroup(filename);
//            break;
//        default:
//            fprintf(stderr,"Solver specified is not supported.\n");
//            break;
//    }
}

void HMMProblemSt::toFileSkill(const char *filename) {
    NPAR nS = this->task->nS, nO = this->task->nO;
    NCAT nK = this->task->nK;
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    // write solved id
    writeSolverInfo(fid, this->task);
    
	fprintf(fid,"Null skill ratios\t");
	for(NPAR m=0; m<nO; m++)
		fprintf(fid," %10.7f%s",this->null_obs_ratios[m],(m==(this->task->nO-1))?"\n":"\t");
    
	NCAT k;
	std::map<NCAT,std::string>::iterator it;
	for(k=0;k<nK;k++) {
		it = this->task->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		NPAR i,j,m;
		fprintf(fid,"PI\t");
		for(i=0; i<nS; i++)
			fprintf(fid,"%12.10f%s",this->param_skill[PI(k,i)],(i==(nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<nS; i++)
			for(j=0; j<nS; j++)
				fprintf(fid,"%12.10f%s",this->param_skill[A(k,i,j)],(i==(nS-1) && j==(nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<nS; i++)
			for(m=0; m<nO; m++)
				fprintf(fid,"%12.10f%s",this->param_skill[B(k,i,m)],(i==(nS-1) && m==(nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblemSt::toFileGroup(const char *filename) {
    NPAR nS = this->task->nS, nO = this->task->nO;
    NCAT nG = this->task->nG;
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    
    // write solved id
    writeSolverInfo(fid, this->task);
    
	fprintf(fid,"Null skill ratios\t");
	for(NPAR m=0; m<nO; m++)
		fprintf(fid," %10.7f%s",this->null_obs_ratios[m],(m==(nO-1))?"\n":"\t");
	NCAT g;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<nG;g++) {
		it = this->task->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		NPAR i,j,m;
        fprintf(fid,"PI\t");
        for(i=0; i<=nS; i++)
            fprintf(fid,"%12.10f%s",this->param_skill[PI(g,i)],(i==(nS-1))?"\n":"\t");
        fprintf(fid,"A\t");
        for(i=0; i<=nS; i++)
            for(j=0; j<=nS; j++)
                fprintf(fid,"%12.10f%s",this->param_skill[A(g,i,j)],(i==(nS-1) && j==(nS-1))?"\n":"\t");
        fprintf(fid,"B\t");
        for(i=0; i<=nS; i++)
            for(m=0; m<=nO; m++)
                fprintf(fid,"%12.10f%s",this->param_skill[B(g,i,m)],(i==(nS-1) && m==(nO-1))?"\n":"\t");
	}
	fclose(fid);
}

//void HMMProblemSt::producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt) {
void HMMProblemSt::producePCorrect(NUMBER*** group_skill_map, NCAT* skills, NPAR n_skills, NUMBER* local_pred, struct context *ctx) {
    NPAR m, i, nO = this->task->nO, nS = this->task->nS;
    NCAT k;
    NUMBER *local_pred_inner = init1D<NUMBER>(nO);
    for(m=0; m<nO; m++) local_pred[m] = 0.0;
    
    for(int l=0; l<n_skills; l++) {
        for(m=0; m<nO; m++) local_pred_inner[m] = 0.0;
        k = skills[l];
        ctx->k=k;
        NUMBER* group_skill_map_gk = group_skill_map[ctx->g][k];
        for(i=0; i<nS; i++) {
            if(/*group_skill_map[ctx->g][k]*/group_skill_map_gk[i]==-1) {
                /*group_skill_map[ctx->g][k]*/group_skill_map_gk[i] = getPI(ctx,i);
            }
            for(m=0; m<nO; m++)
                local_pred_inner[m] += /*group_skill_map[ctx->g][k]*/group_skill_map_gk[i] * getB(ctx,i,m);
        }
        for(m=0; m<nO; m++)
            local_pred[m] += local_pred_inner[m];
    }
    if(n_skills>1) {
        for(m=0; m<nO; m++)
            local_pred[m] /= n_skills;
    }
    projectsimplex(local_pred, nO);
    free(local_pred_inner);
}

void HMMProblemSt::updateValuesLocal(NUMBER*** group_skill_map, NCAT* skills, NPAR n_skills, NUMBER* local_pred, struct context *ctx) {
    NPAR i, j, nS = this->task->nS;//, o = this->task->dat_obs[ctx->t];
    NCAT k, g = ctx->g;
    NUMBER pLe_denom = 0.0;
    NUMBER *pLe = init1D<NUMBER>(nS);
    NUMBER* group_skill_map_gk; // TODO, try to see time change
    
    // update pL
    for(NPAR l=0; l<n_skills; l++) {
        k = skills[l];
        ctx->k = k;
        group_skill_map_gk = group_skill_map[g][k]; // TODO, try to see time change
//        NUMBER* pLbit = gsm(g,k); //BOOST
        
        if(ctx->o>-1) { // known observations //
            // update p(L)
            pLe_denom = 0.0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<nS; i++)
                pLe_denom += /*group_skill_map[g][k]*/group_skill_map_gk[i] * getB(ctx,i,ctx->o);  // TODO: this is local_pred[o]!!!//UNBOOST
//                  pLe_denom += pLbit[i] * getB(ctx,i,o); //BOOST
            for(i=0; i<nS; i++)
                pLe[i] = /*group_skill_map[g][k]*/group_skill_map_gk[i] * getB(ctx,i,ctx->o) / safe0num(pLe_denom); //UNBOOST
//                  pLe[i] = pLbit[i] * getB(ctx,i,o) / safe0num(pLe_denom); //BOOST
            // 2. L = (pLe'*A)';
            for(i=0; i<nS; i++)
            /*group_skill_map[g][k]*/group_skill_map_gk[i]= 0.0; //UNBOOST
//                  pLbit[i]= 0.0; //BOOST
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
            /*group_skill_map[g][k]*/group_skill_map_gk[j] += pLe[i] * getA(ctx,i,j);//A[i][j]; //UNBOOST
//                    pLbit[j] += pLe[i] * getA(ctx,i,j);//A[i][j]; //BOOST
        } else { // unknown observation
            // 2. L = (pL'*A)';
            for(i=0; i<nS; i++)
                pLe[i] = /*group_skill_map[g][k]*/group_skill_map_gk[i]; // copy first; //UNBOOST
//                  pLe[i] = pLbit[i]; // copy first; //BOOST
            for(i=0; i<nS; i++)
            /*group_skill_map[g][k]*/group_skill_map_gk[i] = 0.0; // erase old value //UNBOOST
//                  pLbit[i] = 0.0; // erase old value //BOOST
            for(j=0; j<nS; j++)
                for(i=0; i<nS; i++)
                    /*group_skill_map[g][k]*/group_skill_map_gk[j] += pLe[i] * getA(ctx,i,j);//UNBOOST
//               pLbit[j] += pLe[i] * getA(ctx,i,j);//BOOST
        }// observations
        projectsimplex(/*group_skill_map[g][k]*/group_skill_map_gk, nS); // addition to make sure there's no side effects //UNBOOST
//            projectsimplex(pLbit, nS); // addition to make sure there's no side effects //BOOST
    }
    free(pLe);
}


//void HMMProblemSt::predict(NUMBER* metrics, const char *filename, NPAR* dat_obs, NCAT *dat_group, NCAT *dat_skill, StripedArray<NCAT*> *dat_multiskill) {
void HMMProblemSt::predict(NUMBER* metrics, const char *filename,
                           //NDAT N, NPAR* dat_obs, NCAT *dat_group, NCAT *dat_skill, NCAT *dat_skill_stacked, NCAT *dat_skill_rcount, NDAT *dat_skill_rix,
                           struct task* task,
                           HMMProblemSt **hmms, NPAR nhmms, NPAR *hmm_idx) {
	NDAT t;
	NPAR i, /*o,*/ h, n, m;
	NPAR nS = /*hmms[0]->*/task->nS, nO = /*hmms[0]->*/task->nO;
    NCAT *ar, nK = /*hmms[0]->*/task->nK, nG = /*hmms[0]->*/task->nG;
	char f_multiskill = /*hmms[0]->*/task->multiskill;
	char f_update_known = /*hmms[0]->*/task->update_known;
	char f_update_unknown = /*hmms[0]->*/task->update_unknown;
	int f_predictions = /*hmms[0]->*/task->predictions;
	int f_metrics_target_obs = /*hmms[0]->*/task->metrics_target_obs;
	for(i=1; i<nhmms; i++) {
		if( nS != hmms[i]->task->nS || nO != hmms[i]->task->nO || nK != hmms[i]->task->nK ||
		   nG != hmms[i]->task->nG || // N != hmms[i]->task->N || N_null != hmms[i]->task->N_null || // all internal N's are different
		   f_multiskill != hmms[i]->task->multiskill ||
		   f_update_known != hmms[i]->task->update_known ||
		   f_update_unknown != hmms[i]->task->update_unknown ||
		   f_predictions != hmms[i]->task->predictions ||
		   f_metrics_target_obs != hmms[i]->task->metrics_target_obs) {
			fprintf(stderr,"Error! One of count variables (N, N_null, nS, nO, nK, nG) or flags (multiskill, predictions, metrics_target_obs, update_known, update_unknown) does not have the same value across multiple models\n");
			exit(1);
		}
	}
	
	// initialize and run
//    NUMBER ***predict_hmm = init1D<NUMBER**>(nhmms);
    NUMBER **metrics_hmm = init2D<NUMBER>(nhmms,6);
    NDAT N_all = 0;
    NDAT N_nnul_all = 0;
    NUMBER loglik_all = 0.0;
    NUMBER loglik_nnul_all = 0.0;
    NUMBER sse_all = 0.0;
    NUMBER sse_nnul_all = 0.0;
    NDAT ncorr_all = 0;
    NDAT ncorr_nnul_all = 0;
    for(NPAR h=0;h<nhmms;h++) {
//        predict_hmm[h] = init2D<NUMBER>(hmms[h]->task->N,nO);
        hmms[h]->computeAlphaAndPOParam(metrics_hmm[h]);
        N_all           += hmms[h]->task->N;
        N_nnul_all      += hmms[h]->task->N_null;
        loglik_all      += metrics_hmm[h][0];
        loglik_nnul_all += metrics_hmm[h][1];
        sse_all         += metrics_hmm[h][2];
        sse_nnul_all    += metrics_hmm[h][3];
        ncorr_all       += (NDAT)metrics_hmm[h][4];
        ncorr_nnul_all  += (NDAT)metrics_hmm[h][5];
    }
    
    // print predictions
    FILE *fid = NULL; // file for storing prediction should that be necessary
    HMMProblemSt *hmm = NULL;
    if(f_predictions>0) {
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr,"Can't write output model file %s\n",filename);
            exit(1);
        }
        for(t=0; t<task->N; t++) {
            h = (hmm_idx!=NULL)?hmm_idx[t]:0;
            //o = dat_obs[t];
            hmm = hmms[h];
            getSkillsAtRow(hmm->task, t, &ar, &n);
            if(ar[0]<0) {// if no skills
                for(m=0; m<nO; m++) fprintf(fid,"%12.10f%s",hmm->null_obs_ratios[m],(m<(nO-1))?"\t":"\n");
            } else { // skills present
                for(m=0; m<nO; m++) {
                    fprintf(fid,"%12.10f%s",hmm->task->dat_predict[t][m],(m<(nO-1))?"\t": ((f_predictions==1)?"\n":"\t") );// if we print states of KCs, continut
                }
                if(f_predictions==2) { // if we print out states of KC's as welll
//                    fprintf(stderr,"WARNING! Support for printing pL is not re-implemented \n");
                    for(int l=0; l<n; l++) { // all KC here
                        NDAT tt = l + hmm->task->dat_skill_rix[t];
                        fprintf(fid,"%12.10f%s",hmm->task->dat_predict_k[tt], (l==(n-1) && l==(n-1))?"\n":"\t"); // nnon boost // if end of all states: end line//UNBOOST
                    }
                }
            } // skills presen
        }
        fclose(fid);
    }// print predictions out

    // if necessary guess the obsevaion using Pi and B
    if(f_update_known=='g') {
        fprintf(stderr,"WARNING! Guessing for local update is not re-implemented \n");
//        NUMBER max_local_pred=0;
//        NPAR ix_local_pred=0;
//        for(m=0; m<nO; m++) {
//            if( local_pred[m]>max_local_pred ) {
//                max_local_pred = local_pred[m];
//                ix_local_pred = m;
//            }
//        }
//        o = ix_local_pred;
    }
    if(f_predictions==3) { // if we print out states of KC's as welll
        fprintf(stderr,"WARNING! Writing prediction after the local update is not re-implemented \n");
//        for(int l=0; l<n; l++) { // all KC here
//            fprintf(fid,"%12.10f%s",group_skill_map[g][ ar[l] ][0], (l==(n-1) && l==(n-1))?"\n":"\t"); // nnon boost // if end of all states: end line//UNBOOST
//            //                    fprintf(fid,"%12.10f%s",gsm(g, ar[l] )[0], (l==(n-1) && l==(n-1))?"\n":"\t"); // if end of all states: end line //BOOST
//        }
    }
    
    // tally metrics
	if(metrics != NULL) {
		metrics[0] = loglik_all;
		metrics[1] = loglik_nnul_all;
		metrics[2] = sqrt(sse_all/N_all);
		metrics[3] = sqrt(sse_nnul_all/(N_all-N_nnul_all));
		metrics[4] = (NUMBER)ncorr_all/N_all;
		metrics[5] = (NUMBER)ncorr_nnul_all/(N_all-N_nnul_all);
	}
	
    // recycle
//    for(int h=0;h<nhmms;h++) {
//        free2D(predict_hmm[h],hmms[h]->task->N);
//    }
//    free(predict_hmm);
    free2D(metrics_hmm, nhmms);

}

void HMMProblemSt::postPredictionSave(const char *filename) {
    return;
}

NCAT HMMProblemSt::getModelParamN() {
    return this->model_param_n;
}

void HMMProblemSt::fit() {
    FitNullSkill(/*loglik_rmse, false do RMSE*/);
    switch(this->task->solver)
    {
        case METHOD_BW: // Baum Welch Method
        case METHOD_GD: // Gradient Descent
        case METHOD_CGD: // Conjugate Gradient Descent
        case METHOD_GDL: // Gradient Descent, Lagrange
        case METHOD_GBB: // Brzilai Borwein Gradient Method
            Cycle();
            break;
        default:
            fprintf(stderr,"Solver specified is not supported.\n");
            break;
    }
}

void HMMProblemSt::FitNullSkill() {
    if(this->task->N_null==0) {
        this->null_obs_ratios[0] = 1; // set first obs to 1, simplex preserved
        return; // 0 loglik
    }
    NDAT t;
    NPAR o, m;
    // count distribution of onservations
    for(t=this->task->first_null_skill; t<=this->task->last_null_skill; t++) {
        o = this->task->dat_obs[ t ];
        this->null_obs_ratios[ o ]++;
    }
    // produce per-observation means
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    for(m=0;m<this->task->nO;m++) {
        this->null_obs_ratios[ m ] /= this->task->N_null;
        if( this->null_obs_ratios[m] > this->null_skill_obs_prob ) {
            this->null_skill_obs_prob = this->null_obs_ratios[m];
            this->null_skill_obs = m;
        }
    }
    this->null_skill_obs_prob = safe01num(this->null_skill_obs_prob); // safety for logging
    // metrics don't need to be computed it's handled in computeAlpha... -- main negative log-ikelihood computation
}

void HMMProblemSt::createGradientScaleSimplexAffordances(FitBitSt *fb) {
    NCAT nK = this->task->nK, nO = this->task->nO, nS = this->task->nS, i;
    fb->sclsmplx_n = nK * (1 + 2 * nS);
    fb->sclsmplx_offsets = init1D<NCAT>(fb->sclsmplx_n);
    fb->sclsmplx_sizes = init1D<NPAR>(fb->sclsmplx_n);
    fb->sclsmplx_blocks = init1D<NCAT>(fb->sclsmplx_n);
    fb->param_blocks = init1D<NCAT>(fb->nPARAM);
    fb->sclsmplx_bound_offsets = init1D<NCAT>(fb->sclsmplx_n);
    NCAT c = 0;
    NDAT offset = 0;
    for(NCAT k=0; k<nK; k++) {
        for(i=0; i<(1+nS); i++) {
            fb->sclsmplx_offsets[c] = offset;
            fb->sclsmplx_bound_offsets[c] = i*nS;
            fb->sclsmplx_sizes[c] = (NPAR)nS;
            fb->sclsmplx_blocks[c] = k; // which [skill] block this simplex is part of
            offset+=nS;
            c++;
        }
        for(i=0; i<nS; i++) {
            fb->sclsmplx_offsets[c] = offset;
            fb->sclsmplx_bound_offsets[c] = (1+nS)*nS + i*nO;
            fb->sclsmplx_sizes[c] = (NPAR)nO;
            fb->sclsmplx_blocks[c] = k; // which [skill] block this simplex is part of
            offset+=nO;
            c++;
        }
        for(i=(k*this->skill1_n);i<((1+k)*this->skill1_n);i++) fb->param_blocks[i]=k;
    }
}

void /*NUMBER*/ HMMProblemSt::Cycle() {
	NCAT k, nK = this->task->nK, skill1_n = this->skill1_n;
    NDAT N = this->task->N;
	//
	// fit all as 1 skill first
	//
	if(this->task->single_skill>0) {
        // save originals and force single skill vvvvv
        char saved_multiskill = this->task->multiskill;
        NCAT saved_nK = this->task->nK;
        NCAT* saved_dat_skill = NULL;
        if(this->task->multiskill==0) {
            saved_dat_skill = init1D<NCAT>(N);
            for(NDAT t=0;t<N; t++) saved_dat_skill[t] = this->task->dat_skill[t];
            toZero1D(this->task->dat_skill, N);
        }
        this->task->nK = 1;
        this->task->multiskill = 0;
        // save originals and force single skill ^^^^^
        
        FitBitSt *fb = new FitBitSt(this->task, this->param_skill, this->gradient_skill, (NDAT)(this->skill1_n*nK) );
        createGradientScaleSimplexAffordances(fb);
        fb->fit_results[0].pO0_prefit = -1.0;

        CycleBit(fb); // single computational bit

        if(!this->task->quiet)
            fb->printFitResult(0);
        delete fb;

        // restore originals and copy single skill vvvvv
        this->task->multiskill = saved_multiskill;
        if(this->task->multiskill==0) {
            for(NDAT t=0; t<N; t++) this->task->dat_skill[t] = saved_dat_skill[t];
            delete saved_dat_skill;
        }
        this->task->nK = saved_nK;
        // copy first slot to all other slots
        for(k=1; k<nK; k++) {
            cpy1D(this->param_skill, &this->param_skill[k*skill1_n], (NDAT)skill1_n);
        }
        // restore originals and copy single skill ^^^^^
	}
	//
	// Main fit
	//
    if(this->task->single_skill!=2){
        FitBitSt *fb = new FitBitSt(this->task, this->param_skill, this->gradient_skill, (NDAT)(this->skill1_n*nK) );
        createGradientScaleSimplexAffordances(fb);
        
        CycleBit(fb); // single cycle bit
        
        for(k=0; k<nK && !this->task->quiet; k++) {
//            loglik += fb->fit_results[k].pO*(fb->fit_results[k].pO>0); // reduction'ed
//            if(!this->task->quiet)
            fb->printFitResult(k);
        }
        fb->printFitResult(-1); // overall report
        // loglik = fb->fit_result.pO;
        delete fb;
    }
    
    return;//loglik;
}

void HMMProblemSt::CycleBit(FitBitSt *fb) {
    NCAT k, nK = this->task->nK;
    NCAT n_conv = 0;
    NDAT n_iter = 0;
    
    // inital copy parameter values to the t-1 slice
    NDAT* sizes_k = initToValue1D<NDAT>(nK, skill1_n);   // for full skill paremter vectors
    NDAT* offsets_k = init1D<NDAT>(nK);  // for full skill parameter vectors
    NCAT* blocks_k = init1D<NDAT>(nK);  // index per block
    NUMBER objval = 0.0;
    for(k=0; k<nK; k++) {
        offsets_k[k] = k*this->skill1_n;
        blocks_k[k] = k;
    }
    
    if(this->task->solver==METHOD_BW) {
        objval = computeAlphaAndPOParam(NULL); // just alpha and p(O|param)
        computeBeta();
        computeXiGamma();
    } else {
        objval = computeGradients();
        if( this->task->solver==METHOD_GD || this->task->solver==METHOD_CGD || this->task->solver==METHOD_GBB )
            fb->scaleGradients(false, this->active_set_block);
    }

    // set p(O|param) for pre-step
    fb->fit_result.pO0 = objval; // overall start before this cycle
    fb->fit_result.pO  = objval; // overall finish after this cycle
    if( fb->fit_result.pO0_prefit == -1.0 ) {
        fb->fit_result.pO0_prefit = objval; // overall before fitting started
    }
    for(k=0; k<nK; k++) {
        if(n_iter==0) {
            fb->fit_results[k].pO0 = this->loglik_k[k];
            fb->fit_results[k].pO  = this->loglik_k[k];
            if( fb->fit_results[k].pO0_prefit == -1.0 ) {
                fb->fit_results[k].pO0_prefit = this->loglik_k[k];
                fb->fit_results[k].ndat = this->ndat_k[k];
            }
        }
    }
    while( n_conv!=nK && n_iter<=this->task->maxiter ) {
        // save t-1, t-2 (tm1, tm2) values before making the step or updating the gradient
        cpy1D(fb->PARAMm1, fb->PARAMm2, fb->nPARAM);
        cpy1D(fb->PARAM, fb->PARAMm1, fb->nPARAM);
        if( this->task->solver==METHOD_CGD || this->task->solver==METHOD_GBB)
            cpy1D(fb->GRAD, fb->GRADm1, fb->nPARAM);
        // for congugate gradient descent – work with direction
        if (this->task->solver==METHOD_CGD) {
            if( n_iter==0 ) {
                cpy1D(fb->GRAD, fb->DIR, fb->nPARAM); // copy to DIR
                negate1D(fb->DIR, fb->nPARAM);
            }
            else
                cpy1D(fb->DIR, fb->DIRm1, fb->nPARAM);// fb->copy(FBS_DIR,  FBS_DIRm1);
        }

        
        // Make the step
        if( this->task->solver==METHOD_GD || (n_iter==0 && this->task->solver==METHOD_CGD)  || (n_iter==1 && this->task->solver==METHOD_GBB) )
            doLinearStep(fb, false/*direction*/);
        else if( this->task->solver==METHOD_CGD )
            doConjugateLinearStep(fb);
        else if( this->task->solver==METHOD_GDL )
            doLagrangeStep(fb);
        else if( this->task->solver==METHOD_GBB )
            doBarzilaiBorweinStep(fb);
        else if( this->task->solver==METHOD_BW )
            doBaumWelchStep(fb);
        
        // project to simplex with boundaries if necessary
        fb->projectToSimplex( this->non01constraints, this->active_set_block );
        // check factual convergence (without considering connectivity yet)
        n_conv = fb->checkConvergence(nK, offsets_k, sizes_k, blocks_k, this->active_set_block);

        // decide on active set (propagate dependencies)
        // this->task->connectivities[0] skill-skill
//        NDAT set0to1 = 0; // DEBUG
        for(k=0; k<nK; k++) {
            if( this->active_set_block[k]==1 ) {
                for(NCAT l=0; l<=nK; l++) {
                    if(l!=k && this->task->connectivities[0][k][l]==1) {
//                        if(this->active_set_block[l]==0) set0to1++; // DEBUG
                        this->active_set_block[l] = 1;
                    }
                } // for all potential skills
            } // if the skill is active
        }
//        fprintf(stderr,"[debug] Activated skills %d.\n", set0to1); // DEBUG

        // copy parameter values after we already compared step t-1 with currently computed step t
        // cpy1D(fb->PARAMm1, fb->PARAMm2, fb->nPARAM);
        // cpy1D(fb->PARAM, fb->PARAMm1, fb->nPARAM);
        // for Gradient Descent METHOD_GD, no extra copying
        // if Conjugate Gradient and Barzilain Borwein – copy previous gradient
//        if( this->task->solver==METHOD_CGD || this->task->solver==METHOD_GBB)
//            cpy1D(fb->GRAD, fb->GRADm1, fb->nPARAM);
//        // for congugate gradient descent – work with direction
//        if (this->task->solver==METHOD_CGD) {
//            if( n_iter==0 ) {
//                cpy1D(fb->GRAD, fb->DIRm1, fb->nPARAM);
//                negate1D(fb->DIRm1, fb->nPARAM);
//            }
//            else
//                cpy1D(fb->DIR, fb->DIRm1, fb->nPARAM);// fb->copy(FBS_DIR,  FBS_DIRm1);
//        }

        // move this into the stepper vvvvvv
//        for(k=0; k<=nK; k++) {
//            fb->fit_results[k].iter = n_iter;
//            //frs[k].pOmid = frs[k].pO;
//        }
        // compute p(O|param) post-step after projecting
        if(this->task->solver==METHOD_BW) {
            objval = computeAlphaAndPOParam(NULL); // just alpha and p(O|param)
            computeBeta();
            computeXiGamma();
        } else {
            objval = computeGradients();
            if( this->task->solver==METHOD_GD || this->task->solver==METHOD_CGD || this->task->solver==METHOD_GBB )
                fb->scaleGradients(false, this->active_set_block);
        }
        
        // update objective function & count iterations
        fb->fit_result.pO = objval;
        for(k=0; k<nK; k++) {
            // update objective
            if(this->active_set_block[k]>0) {
                fb->fit_results[k].pO0 = fb->fit_results[k].pO;
                fb->fit_results[k].pO = this->loglik_k[k];
                fb->fit_results[k].iter++;
            }
        }
        n_iter++;
    }// single skill loop
    fb->fit_result.iter = n_iter;
    // this is the case when 2 steps lead to step back to the initial value, i.e. oscillation
    // the code below doesn' make anny sense yet to me.s
//    for(k=0; k<=nK; k++) {
//        if(fb->fit_results[k].iter==2 && fb->fit_results[k].conv && fb->checkConvergenceBitSingle(offsets_k[k], sizes_k[k], &fb->fit_results[0])) {
//            // decrease iteration counter to 1
//            fb->fit_results[k].iter--;
//        }
//        fb->fit_results[k].pO = this->loglik_k[k];
//    }
    free(offsets_k);
    free(sizes_k);
    free(blocks_k);
    return;
}

void HMMProblemSt::doLagrangeStep(FitBitSt *fb) {
    // goal: for PI, and every row of A, and B multiply by negative gradient and average over the sum
        
    NUMBER buf;
    //NCAT active = 0;
    for(NCAT l=0; l<fb->sclsmplx_n; l++) { // for all vectors (simplexes) that are supposed to sum to 1
        if(fb->sclsmplx_blocks[l]==0) continue; // if the simplex belongs to block that is not in active set, skip
        buf = 0;
        for(NCAT i=fb->sclsmplx_offsets[l]; i<(fb->sclsmplx_offsets[l]+fb->sclsmplx_sizes[l]); i++)
            //if( (this->active_set_param[i] + this->unblocked_param[i]) == 2 ) // in active set and unblocked
                buf -= fb->PARAM[i]*fb->GRAD[i];
        for(NCAT i=fb->sclsmplx_offsets[l]; i<(fb->sclsmplx_offsets[l]+fb->sclsmplx_sizes[l]); i++)
            //if( (this->active_set_param[i] + this->unblocked_param[i]) == 2 ) // in active set and unblocked
                fb->PARAM[i] = (buf>0)? (-fb->PARAM[i]*fb->GRAD[i])/buf : fb->PARAM[i];
    }
} // doLagrangeStep

void HMMProblemSt::doBarzilaiBorweinStep(FitBitSt *fb) {
    // Barzilai Borwein value: s' * s / ( s' * (g-g_m1) )
    
    NPAR i;
    NCAT l,k, nK = this->task->nK;

    // first scale down gradients -- MOVED TO Cycle Bit
    // fb->scaleGradients(false, this->active_set_block);
    
    // compute s_k_m1
    NUMBER *s_k_m1 = init1D<NUMBER>(fb->nPARAM);
    NUMBER *alpha_step = init1D<NUMBER>(nK);
    NUMBER *alpha_denom = init1D<NUMBER>(nK);
    for(l=0; l<fb->nPARAM; l++){
        s_k_m1[l] = fb->PARAMm1[l] - fb->PARAMm2[l];
    }

    // compute alpha_step
    for(k=0; k<nK; k++) {
        if (this->active_set_block[k]==0) continue; // if not in the active set, skip
        for(i=0; i<skill1_n; i++) {
            l = skill1_n*k + i;
            alpha_step[k] += s_k_m1[l]*s_k_m1[l];
            alpha_denom[k] += s_k_m1[l]*(fb->GRAD[l] - fb->GRADm1[l]);;
        }
        alpha_step[k] = alpha_step[k] / safe0num(alpha_denom[k]);
    }

    // make the step
    for(k=0; k<nK; k++) {
        if (this->active_set_block[k]==0) continue; // if not in the active set, skip
        for(i=0; i<skill1_n; i++) {
            l = skill1_n*k + i;
            fb->PARAM[l] -= alpha_step[k] * fb->GRAD[l];
        }
    }
    
    // recycle
    free(s_k_m1);
    free(alpha_step);
    free(alpha_denom);
} // doBarzilaiBorweinStep

void HMMProblemSt::doBaumWelchStep(FitBitSt *fb) {
	NCAT *ar, k, l;
    NPAR nS = this->task->nS, nO = this->task->nO;
	NPAR i,j,m, o, n;
	NDAT t, tt, offset;
    
    // use PARAM for new parameter numerator value, use PARAMcopy for new parameter denominator value
    fb->PARAMcopy = init1D<NUMBER>(fb->nPARAM);
    //toZero1D(fb->PARAM, fb->nPARAM); // ONLY ZERO ACTIVE
    for(l=0; l<fb->nPARAM; l++)
        if( this->active_set_block[ fb->param_blocks[l] ]>0 )
            fb->PARAM[l] = 0;
    
    // compute numerators in PARAM and denominators in PARAMcopy
    for(t=0; t<this->task->N; t++) {
        o = this->task->dat_obs[t]; // observation y in the data is 1-right, 0-wrong; HMM-Scalable assumes 1-right, 2-wrong, with -1 taken out, so 0-right, 1-wrong
        // grab skill array (if exists)
        getSkillsAtRow(this->task, t, &ar, &n);

        NCAT n_active = 0;
        for(int l=0; l<n; l++) n_active+=(this->active_set_block[ar[l]]>0);
                
        for(l=0; l<n && n_active==n/**/; l++) {
            k = ar[l];
            tt = l + this->task->dat_skill_rix[t]; // tt – index into stacked skill array
            
            // pi
            offset = k*this->skill1_n;
            for(i=0; i<nS; i++) {
                if(this->backward_ix[tt]==-1) { // first row of skill's sequence.
                    fb->PARAM[offset+i] += this->gamma[tt][i]; // num
                    fb->PARAMcopy[offset+i]++ ;// den -- count first rows
                }
            }
            // A
            offset = k*this->skill1_n + nS;
            for(i=0; i<nS; i++) {
                for(j=0; j<nS; j++) {
                    fb->PARAM[offset+i*nS+j] += this->xi[tt][i][j];
                    fb->PARAMcopy[offset+i*nS+j] += this->gamma[tt][i];
                }
            }
            // B
            offset = k*this->skill1_n + (1+nS)*nS;
            for(i=0; i<nS; i++) {
                for(m=0; m<nO; m++) {
                    fb->PARAM[offset+i*nO+m] += (m==o) * this->gamma[tt][i];
                    fb->PARAMcopy[offset+i*nO+m] += this->gamma[tt][i];
                }
            }
        }// for all skills
        
    }// for all rows
    
    // divide numerators and denominators where needed
    for(l=0; l<fb->nPARAM; l++) { //
        if( this->active_set_block[ fb->param_blocks[l] ]>0 ) {
            fb->PARAM[l] /= safe0num( fb->PARAMcopy[l] );
        }
    }// all parameters

    // recycle
    free(fb->PARAMcopy);
    fb->PARAMcopy = NULL;
}

void HMMProblemSt::doLinearStep(FitBitSt *fb, bool direction) {
    NPAR i;
    NCAT l,k, nK = this->task->nK, skill1_n = this->skill1_n;
   // NDAT c = 0;
    
    // MOVED TO Cycle Bit
    // fb->scaleGradients(direction, this->active_set_block); // it knows whether to scale gradient or direction
    
    fb->PARAMcopy = init1D<NUMBER>(fb->nPARAM);
    fb->GRADcopy = init1D<NUMBER>(fb->nPARAM); // copy of gradient or direction, we will need just the grad's copy saved
    NUMBER * save_GRAD = init1D<NUMBER>((NDAT)fb->nPARAM); // save GRAD since while stepping it would be recomputed
    
    // all of the variables were single value, now we ad `_s` suffix and make them arrays
    NUMBER* e_s = initToValue1D<NUMBER>(nK,this->task->ArmijoSeed); // step seed
    bool* compliesArmijo_s = initToValue1D<bool>(nK, false);
    bool* compliesWolfe2_s = initToValue1D<bool>(nK, false); // second wolfe condition is turned on, if satisfied - honored, if not, just the 1st is used
//    NCAT n_compliesArmijo_or_Wolfe2 = 0;
//    NCAT n_e_lt_ArmijoMinStep = 0;
    NCAT n_compliesArmijo_or_compliesWolfe2_or_n_e_lt_ArmijoMinStep = 0;

    NUMBER* e_Armijo_s = initToValue1D<NUMBER>(nK,-1); // e (step size) at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER* f_xkplus1_Armijo_s = initToValue1D<NUMBER>(nK,0); // f_xkplus1_Armijo at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER* f_xk_s = initToValue1D<NUMBER>(nK,0);
    for(k=0; k<nK; k++) f_xk_s[k] = this->loglik_k[k];
    NUMBER* f_xkplus1_s = initToValue1D<NUMBER>(nK,0);
    NUMBER* p_k_by_neg_p_k_s = initToValue1D<NUMBER>(nK,0);
    NUMBER* p_k_by_neg_p_kp1_s = initToValue1D<NUMBER>(nK,0);

    cpy1D(fb->PARAM,fb->PARAMcopy,fb->nPARAM);
    if(!direction)
        cpy1D(fb->GRAD,fb->GRADcopy,fb->nPARAM);
    else {
        cpy1D(fb->DIR,fb->GRADcopy,fb->nPARAM); // in case we are using DIR, GRADcopy is the DIR
        // negate1D(fb->GRADcopy, fb->nPARAM); // negate, because we treat is as GRADient, direction has inverse direction
        // ^^ was already negated in the Conjugate part
    }
    cpy1D(fb->GRAD, save_GRAD, fb->nPARAM); // copy grads
    // compute p_k * -p_k for all k
    
    NCAT nK_active = 0, nK_active_local = 0 /*set to 1 later to be changed*/;
    
    for(k=0; k<nK; k++) {
        nK_active += this->active_set_block[k]>0;
        if(this->active_set_block[k]>0) {
            for(i=0; i<skill1_n; i++) {
                l = skill1_n*k + i;
                p_k_by_neg_p_k_s[k] -= fb->GRAD[l] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
            }
        }
    } // vvvv REDO
//    for(l=0; l<fb->nPARAM; l++) {
//        k = fb->param_blocks[l];
//        nK_active += this->active_set_block[k]>0;
//        if(this->active_set_block[k]>0) {
//            p_k_by_neg_p_k_s[k] -= fb->GRAD[l] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
//        }
//    }
    nK_active_local = nK_active;
    
    // copy active set for original values would be used locally
    NPAR * active_set_block_local = init1D<NPAR>((NDAT)nK);
    cpy1D(this->active_set_block, active_set_block_local, nK);
    
    int iter = 0; // limit iter steps to 20, via ArmijoMinStep (now 10)
    //while( n_compliesArmijo_or_Wolfe2<nK_active/*nK*/ && n_e_lt_ArmijoMinStep<nK_active/*nK*/ ) {
    while( nK_active_local > 0) { // n_compliesArmijo_or_compliesWolfe2_or_n_e_lt_ArmijoMinStep < nK_active/*nK*/) {
        nK_active_local = 0;
        // n_compliesArmijo_or_Wolfe2 = 0;
        // n_e_lt_ArmijoMinStep = 0;
        //n_compliesArmijo_or_compliesWolfe2_or_n_e_lt_ArmijoMinStep = 0;
        // step
//        c = 0;
//        for(k=0; k<nK; k++) {
//            if(this->active_set_block[k]>0) {
//                for(i=0; i<skill1_n; i++) {
//                    l = skill1_n*k + i;
//                    fb->PARAM[l] = fb->PARAMcopy[l] - e_s[k] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
////                    if(fb->PARAM[l]<0 || fb->PARAM[l]>1) {
////                        c++;
////                    }
//                }
//            }
//        } // vvvv REDO
        for(l=0; l<fb->nPARAM; l++) {
            k = fb->param_blocks[l];
            if(this->active_set_block[k]>0) {
                fb->PARAM[l] = fb->PARAMcopy[l] - e_s[k] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
            }
        }
        
//        if(c>0) {
//            fprintf(stderr, "ERROR! %d values are not within [0, 1] range!\n",c);
//        }

        // project parameters to simplex if needs be --- it's necessary here, since we might make frequent step back-tracks
        fb->projectToSimplex( this->non01constraints, this->active_set_block );
        
        // recompute gradients
        computeGradients();
        fb->scaleGradients(false, this->active_set_block);
        
        NCAT nK_deactiv = 0;
        NCAT nK_reactiv = 0;
        for(k=0; k<nK; k++) {
            if(this->active_set_block[k]>0) {
                // compute f(x_{k+1})
                f_xkplus1_s[k] = this->loglik_k[k];
                // compute Armijo compliance
                compliesArmijo_s[k] = (f_xkplus1_s[k] <= (f_xk_s[k] + (this->task->ArmijoC1 * e_s[k] * p_k_by_neg_p_k_s[k])));
                // compute Wolfe 2
                for(i=0; i<skill1_n; i++) {
                    l = skill1_n*k + i;
                    p_k_by_neg_p_kp1_s[k] -= fb->GRAD[l] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
                }
                compliesWolfe2_s[k] = (p_k_by_neg_p_kp1_s[k] >= this->task->ArmijoC2 * p_k_by_neg_p_k_s[k]);
                
                // decide
                if( compliesArmijo_s[k] && e_Armijo_s[k]==-1 ){
                    e_Armijo_s[k] = e_s[k]; // save the first time Armijo is statisfied, in case we'll roll back to it when Wolfe 2 is finally not satisfied
                    f_xkplus1_Armijo_s[k] = f_xkplus1_s[k];
                }
                // n_compliesArmijo_or_Wolfe2 += (compliesArmijo_s[k] || compliesWolfe2_s[k]);
                // n_e_lt_ArmijoMinStep += (e_s[k] <= this->task->ArmijoMinStep);
                if(compliesArmijo_s[k] || (compliesArmijo_s[k] && compliesWolfe2_s[k]) || (e_s[k] <= this->task->ArmijoMinStep)) {
                    this->active_set_block[k] = 0;
                    //nK_active--;
                    nK_deactiv++;
                    n_compliesArmijo_or_compliesWolfe2_or_n_e_lt_ArmijoMinStep++;// += (compliesArmijo_s[k] || compliesWolfe2_s[k] || (e_s[k] <= this->task->ArmijoMinStep));
                }
                else
                    e_s[k] /= this->task->ArmijoReduceFactor;
            }// if active
        }
        // rehash active set using reachability/connectivity
        // decide on active set (propagate dependencies)
        // this->task->connectivities[0] skill-skill
        for(k=0; k<nK; k++) {
            if(this->active_set_block[k]>0) {
                for(NCAT l=0; l<=nK; l++) {
                    if(l!=k && this->task->connectivities[0][k][l]==1 && this->active_set_block[l]==0) {
                        this->active_set_block[l] = 1;
                        //nK_active++;
                        nK_reactiv++;
                    }
                } // for all potential skills
            }
        }
        // coun active
        for(k=0; k<nK; k++) nK_active_local += this->active_set_block[k]>0;

        iter++;
    } // armijo loop
    
    // reinstate original active set
    cpy1D(active_set_block_local, this->active_set_block, nK);
    // reinstate original gradients
    cpy1D(save_GRAD, fb->GRAD, nK);

//    c = 0;
    for(k=0; k<nK; k++) {
        if(this->active_set_block[k]>0) {
            if(!compliesArmijo_s[k]) { // we couldn't step away from current, copy the inital point back
                e_s[k] = 0;
                f_xkplus1_s[k] = f_xk_s[k];
                cpy1D(&fb->PARAMcopy[k*skill1_n],&fb->PARAM[k*skill1_n],skill1_n);
            } else {
                if(compliesArmijo_s[k] && !compliesWolfe2_s[k]) { // we couldn't step away from current, copy the inital point back
                    e_s[k] = e_Armijo_s[k]; // return the first Armijo-compliant e
                    f_xkplus1_s[k] = f_xkplus1_Armijo_s[k]; // return the first Armijo-compliant f_xkplus1
                    // create new versions of FBS_PAR using e_Armijo as a step
//                    if(!direction) { // if step on GRAD not on DIR -- this is useless, why do it here
//                        cpy1D(&fb->GRADcopy[k*skill1_n],&fb->GRAD[k*skill1_n],skill1_n);
//                    }
                    // step
                    for(i=0; i<skill1_n; i++) {
                        l = skill1_n*k + i;
                        fb->PARAM[l] = fb->PARAMcopy[l] - e_s[k] * fb->GRADcopy[l]; // in case we are using DIR, GRADcopy is the DIR
//                        if(fb->PARAM[l]<0 || fb->PARAM[l]>1) {
//                            c++;
//                        }
                    }
                    // ^^^^^ end of create new versions of FBS_PAR using e_Armijo as a step
                }
            }
        }// if in active set
    }
    // project to simplex
//    if(c>0) {
//        fprintf(stderr, "ERROR! %d values are not within [0, 1] range!\n",c);
//    }
    fb->projectToSimplex( this->non01constraints, this->active_set_block );

//    nK_active = 0;                                                // DEBUG check for reinstatement
//    for(k=0; k<nK; k++) nK_active += this->active_set_block[k]>0; // DEBUG check for reinstatement
    
    // free memory (1)
    free(e_s);
    free(compliesArmijo_s);
    free(compliesWolfe2_s);
    free(e_Armijo_s);
    free(f_xkplus1_Armijo_s);
    free(f_xk_s);
    free(f_xkplus1_s);
    free(p_k_by_neg_p_k_s);
    free(p_k_by_neg_p_kp1_s);
    // free memory (2)
    free(active_set_block_local);
    free(fb->PARAMcopy);
    free(fb->GRADcopy);
    free(save_GRAD);
    fb->PARAMcopy = NULL;
    fb->GRADcopy = NULL;
    return;
} // doLinearStep

void HMMProblemSt::doConjugateLinearStep(FitBitSt *fb) {
    NPAR i=0;
    NCAT nK = this->task->nK, k, l;
   
    // first scale down gradients -- MOVED TO Cycle Bit
    // fb->scaleGradients(false, this->active_set_block);
    
    // compute beta_gradient_direction
    NUMBER *beta_grad = init1D<NUMBER>(nK);
    NUMBER *beta_grad_den = init1D<NUMBER>(nK);

    for(k=0; k<nK; k++) {
        if(this->active_set_block[k]>0) {
            for(i=0; i<skill1_n; i++) {
                l = skill1_n*k + i;
                if(this->task->solver_setting == 1) { // Polak-Ribiere
                    beta_grad[k]     += fb->GRAD[l]*fb->GRAD[l];
                    beta_grad_den[k] += fb->GRADm1[l]*fb->GRADm1[l];
                } else if(this->task->solver_setting == 2) { // Fletcher-Reeves
                    beta_grad[k]     += -fb->GRAD[l]*(-fb->GRAD[l] + fb->GRADm1[l]);
                    beta_grad_den[k] +=  fb->GRADm1[l]*fb->GRADm1[l];
                } else if(this->task->solver_setting == 3) { // Hestenes-Stiefel
                    beta_grad[k]     += fb->GRAD[l]* (-fb->GRAD[l] + fb->GRADm1[l]); // no -, since neg gradient and - is +
                    beta_grad_den[k] += fb->DIRm1[l]*(-fb->GRAD[l] + fb->GRADm1[l]);
                } else if(this->task->solver_setting == 4) { // Dai-Yuan
                    beta_grad[k]     += -fb->GRAD[l]*fb->GRAD[l];
                    beta_grad_den[k] +=  fb->DIRm1[l]*(-fb->GRAD[l] + fb->GRADm1[l]);
                } else {
                    fprintf(stderr,"Wrong stepping algorithm (%d) specified for Conjugate Gradient Descent\n",this->task->solver_setting);
                    exit(1);
                }
            }
            beta_grad[k] = beta_grad[k] / safe0num(beta_grad_den[k]);
            beta_grad[k] = (beta_grad[k]>=0)?beta_grad[k]:0;
        } // if in the active set
    }// for all skills
    
    // compute new negative direction
    toZero1D(fb->DIR, fb->nPARAM);
    for(k=0; k<nK; k++) {
        if(this->active_set_block[k]>0) {
            for(i=0; i<skill1_n; i++) {
                l = skill1_n*k + i;
//              fb->DIR[l] = -fb->GRAD[l] + beta_grad[k] * fb->DIRm1[l];
                fb->DIR[l] =  fb->GRAD[l] - beta_grad[k] * fb->DIRm1[l]; // negative direction
            }
        }
    }
    fb->scaleGradients(true, this->active_set_block);

    // almost identical code from linear step, but using direction
    doLinearStep(fb, true /*direction*/);
     
    // recycle
    free(beta_grad);
    free(beta_grad_den);
} // doConjugateLinearStep

void HMMProblemSt::readNullObsRatio(FILE *fid, struct task* task, NDAT *line_no) {
	NPAR i;
	//
	// read null skill ratios
	//
    fscanf(fid, "Null skill ratios\t");
    this->null_obs_ratios =Calloc(NUMBER, (size_t)this->task->nO);
    this->null_skill_obs      = 0;
    this->null_skill_obs_prob = 0;
	for(i=0; i<task->nO; i++) {
        if( i==(task->nO-1) ) // end
            fscanf(fid,"%lf\n",&this->null_obs_ratios[i] );
        else
            fscanf(fid,"%lf\t",&this->null_obs_ratios[i] );
		
        if( this->null_obs_ratios[i] > this->null_skill_obs_prob ) {
            this->null_skill_obs_prob = this->null_obs_ratios[i];
            this->null_skill_obs = i;
        }
	}
    (*line_no)++;
}

void HMMProblemSt::readModel(const char *filename, bool overwrite) {
	FILE *fid = fopen(filename,"r");
	if(fid == NULL)
	{
		fprintf(stderr,"Can't read model file %s\n",filename);
		exit(1);
	}
	int max_line_length = 1024;
	char *line = Malloc(char,(size_t)max_line_length);
	NDAT line_no = 0;
    struct task inittask;
    set_task_defaults(&inittask);
    
    //
    // read solver info
    //
    if(overwrite)
        readSolverInfo(fid, this->task, &line_no);
    else
        readSolverInfo(fid, &inittask, &line_no);
    //
    // read model
    //
    readModelBody(fid, &inittask, &line_no, overwrite);
		
	fclose(fid);
	free(line);
}

void HMMProblemSt::readModelBody(FILE *fid, struct task* task, NDAT *line_no,  bool overwrite) {
    NPAR i,j,m, nS=this->task->nS, nO=this->task->nO;
    NCAT nK=task->nK;
	NCAT k = 0, idxk = 0, offset;
    std::map<std::string,NCAT>::iterator it;
	string s;
    char col[2048];
    //
    readNullObsRatio(fid, task, line_no);
    //
    // init param
    //
    if(overwrite) {
        this->task->map_group_fwd = new map<string,NCAT>();
        this->task->map_group_bwd = new map<NCAT,string>();
        this->task->map_skill_fwd = new map<string,NCAT>();
        this->task->map_skill_bwd = new map<NCAT,string>();
    }
	//
	// read skills
	//
	for(k=0; k<nK; k++) {
		// read skill label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;
        if(overwrite) {
            this->task->map_skill_fwd->insert(pair<string,NCAT>(s, (NCAT)this->task->map_skill_fwd->size()));
            this->task->map_skill_bwd->insert(pair<NCAT,string>((NCAT)this->task->map_skill_bwd->size(), s));
            idxk = k;
        } else {
            it = this->task->map_skill_fwd->find(s);
            if( it==this->task->map_skill_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid, "%*[^\n]\n");
                fscanf(fid, "%*[^\n]\n");
                fscanf(fid, "%*[^\n]\n");
                (*line_no)+=3;
                continue; // skip this iteration
            }
            else
                idxk = it->second;
        }
        // read PI
        fscanf(fid,"PI\t");
        offset = idxk*this->skill1_n;
        for(i=0; i<(nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            // this->pi[idxk][i] = atof(col);
            this->param_skill[offset+i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        // this->pi[idxk][i] = atof(col);
        this->param_skill[offset+i] = atof(col);
        (*line_no)++;
		// read A
        fscanf(fid,"A\t");
        offset = idxk*this->skill1_n + nS;
		for(i=0; i<nS; i++)
			for(j=0; j<nS; j++) {
                if(i==(nS-1) && j==(nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    // this->A[idxk][i][j] = atof(col);
                    this->param_skill[offset+i*nS + j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    // this->A[idxk][i][j] = atof(col);
                    this->param_skill[offset+i*nS + j] = atof(col);
                }
			}
        (*line_no)++;
		// read B
        fscanf(fid,"B\t");
        offset = idxk*this->skill1_n + (1+nS)*nS;
		for(i=0; i<nS; i++)
			for(m=0; m<nO; m++) {
                if(i==(nS-1) && m==(nO-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    // this->B[idxk][i][m] = atof(col);
                    this->param_skill[offset+i*nO + m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not last one
                    // this->B[idxk][i][m] = atof(col);
                    this->param_skill[offset+i*nO + m] = atof(col);
                }
			}
        (*line_no)++;
	} // for all k
}

// helper function
void HMMProblemSt::toFileAlphaBetaObs(const char *filename, NCAT k) {
    NDAT t, tt=0;
    NPAR nS = this->task->nS;
    FILE *fid = fopen(filename,"w");
    
    if(fid == NULL) {
        fprintf(stderr,"Can't write output model file %s\n",filename);
        exit(1);
    }
    
    // first time skill k starts
    while(this->task->dat_skill_stacked[tt]!=k) {
        tt++;
    }
    
    do {
        t = tt;
        for(NPAR i=0; i<nS; i++) fprintf(fid,"%20.17f\t",this->alpha[tt][i]);
        for(NPAR i=0; i<nS; i++) fprintf(fid,"%20.17f\t",this->beta[tt][i]);
        fprintf(fid,"%i\n",this->task->dat_obs_stacked[tt]);
        tt = this->forward_ix[tt];
    }
    while(this->forward_ix[t]!=-1);
    fclose(fid);
}
