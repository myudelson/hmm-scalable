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

#include "HMMProblemEloK.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemEloK::HMMProblemEloK(struct param *param) {
    //    this->sizes = {param->nK, param->nK, param->nK};
    this->sizes[0] = param->nK;
    this->sizes[1] = param->nK;
    this->sizes[2] = param->nK;
    this->n_params = param->nK * 4 + 1;
    init(param);
}

void HMMProblemEloK::init(struct param *param) {
	this->p = param;
	this->non01constraints = true;
    this->null_obs_ratio = Calloc(NUMBER, (size_t)this->p->nO);
    this->neg_log_lik = 0;
    this->null_skill_obs = 0;
    this->null_skill_obs_prob = 0;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK; NCAT nG = this->p->nG;
    
    NUMBER *a_PI, ** a_A, ** a_B;
    init3Params(a_PI, a_A, a_B, nS, nO);
    
    //
    // setup params
    //
	NPAR i, j, m, idx, offset;
	NUMBER sumPI = 0;
	NUMBER sumA[nS];
	NUMBER sumB[nS];
	for(i=0; i<nS; i++) {
		sumA[i] = 0;
		sumB[i] = 0;
	}
    
    // write default parameters first
	// populate PI
	for(i=0; i<((nS)-1); i++) {
		a_PI[i] = this->p->init_params[i];
		sumPI  += this->p->init_params[i];
	}
	a_PI[nS-1] = 1 - sumPI;
	// populate A
	offset = (NPAR)(nS-1);
	for(i=0; i<nS; i++) {
		for(j=0; j<((nS)-1); j++) {
			idx = (NPAR)(offset + i*((nS)-1) + j);
			a_A[i][j] = this->p->init_params[idx];
			sumA[i]  += this->p->init_params[idx];
		}
		a_A[i][((nS)-1)]  = 1 - sumA[i];
	}
	// polupale B
	offset = (NPAR)((nS-1) + nS*(nS-1));
	for(i=0; i<nS; i++) {
		for(m=0; m<((nO)-1); m++) {
			idx = (NPAR)(offset + i*((nO)-1) + m);
			a_B[i][m] = this->p->init_params[idx];
			sumB[i] += this->p->init_params[idx];
		}
		a_B[i][((nO)-1)]  = 1 - sumB[i];
	}
    // populate K
    offset = (NPAR)( (nS-1) + nS*(nS-1) + nS*(nO-1) );
    this->K = sigmoid(this->p->init_params[offset]); // keep Elo parameters in probabilistic form
    
    // mass produce PI's/PIg's, A's, B's
    if( this->p->do_not_check_constraints==0 && !checkPIABConstraints(a_PI, a_A, a_B)) {
        fprintf(stderr,"params do not meet constraints.\n");
        exit(1);
    }
    this->pi  = init2D<NUMBER>(nK, nS);
    this->A   = init3D<NUMBER>(nK, nS, nS);
    this->B   = init3D<NUMBER>(nK, nS, nO);
    NCAT x;
    for(x=0; x<nK; x++) {
        cpy1D<NUMBER>(a_PI, this->pi[x], nS);
        cpy2D<NUMBER>(a_A,  this->A[x],  nS, nS);
        cpy2D<NUMBER>(a_B,  this->B[x],  nS, nO);
    }
    // destroy setup params
	free(a_PI);
	free2D<NUMBER>(a_A, nS);
	free2D<NUMBER>(a_B, nS);
	
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
			idx = (NPAR)(offset + i*nS + j);
			lbA[i][j] = this->p->param_lo[idx];
			ubA[i][j] = this->p->param_hi[idx];
		}
	// *B
	offset = (NPAR)(nS + nS*nS);
	for(i=0; i<nS; i++)
		for(j=0; j<nO; j++) {
			idx = (NPAR)(offset + i*nS + j);
			lbB[i][j] = this->p->param_lo[idx];
			ubB[i][j] = this->p->param_hi[idx];
		}
    // set Elo values and counts
    this->elo_count_g = init1D<NCAT>((NDAT)nG);
    this->elo_track_g = init1D<NUMBER>((NDAT)nG);
    this->elo_track_g_t = init1D<NUMBER>((NDAT)this->p->N);
}

HMMProblemEloK::~HMMProblemEloK() {
    destroy();
}

void HMMProblemEloK::destroy() {
	// destroy model data - only additional stuff
    this->K = 0.0;
    if(this->elo_track_g != NULL) free(this->elo_track_g);
    if(this->elo_count_g != NULL) free(this->elo_count_g);
    if(this->elo_track_g_t != NULL) free(this->elo_track_g_t);
}// ~HMMProblemEloK

NUMBER** HMMProblemEloK::getPI() {
	return this->pi;
}

NUMBER*** HMMProblemEloK::getA() {
	return this->A;
}

NUMBER*** HMMProblemEloK::getB() {
	return this->B;
}


NUMBER* HMMProblemEloK::getPI(NCAT x) { // same as getPIk(x)
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->pi[x];
}

NUMBER** HMMProblemEloK::getA(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing A_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->A[x];
}

NUMBER** HMMProblemEloK::getB(NCAT x) {
	if( x > (this->p->nK-1) ) {
		fprintf(stderr,"While accessing B_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
		exit(1);
	}
	return this->B[x];
}

NUMBER HMMProblemEloK::getPI(struct data* dt, NPAR i) {
    bool q_neg = (i==0)?false:true;
    NUMBER p = this->pi[dt->k][i], q = this->elo_track_g_t[dt->t];
    return pairing(p,q,q_neg);
    //    return sigmoid( logit( this->pi[dt->k][i] ) + logit( this->PIg[dt->g][i] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemEloK::getA(struct data* dt, NPAR i, NPAR j) {
    bool q_neg = (j==0)?false:true;
    NUMBER p = this->A[dt->k][i][j], q = this->elo_track_g_t[dt->t];
    return pairing(p,q,q_neg);
    //    return sigmoid( logit( this->A[dt->k][i][j] ) + logit( this->Ag[dt->k][i][j] ) );
}

// getters for computing alpha, beta, gamma
NUMBER HMMProblemEloK::getB(struct data* dt, NPAR i, NPAR m) {
    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
    if(m<0)
        return 1;
    bool q_neg = (m==0)?false:true;
    NUMBER p = this->B[dt->k][i][m], q = this->elo_track_g_t[dt->t];
    return pairing(p,q,q_neg);
}

void HMMProblemEloK::setGradPI(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0, ndat = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    //    o = dt->obs[t];
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        o = this->p->dat_obs[ dt->ix[t] ];
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->pi[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = 1 / safe0num( fb->pi[i] * (1-fb->pi[i]) );
            fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:fb->B[i][o]) / safe0num(dt->p_O_param);
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_PI, this->p, (NUMBER)ndat);
}

void HMMProblemEloK::setGradA (FitBit *fb){
    if(this->p->block_fitting[1]>0) return;
    NDAT t, ndat = 0;
    NPAR o, i, j;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        for(t=1; t<dt->n; t++) {
            o = this->p->dat_obs[ dt->ix[t] ];
            for(i=0; i<fb->nS; i++)
                for(j=0; j<fb->nS; j++) {
                    combined = getA(dt,i,j);
                    deriv_logit = 1 / safe0num( fb->A[i][j] * (1-fb->A[i][j]) );
                    fb->gradA[i][j] -= combined * (1-combined) * deriv_logit * dt->beta[t][j] * ((o<0)?1:fb->B[j][o]) * dt->alpha[t-1][i] / safe0num(dt->p_O_param) ;
                }
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_A, this->p, (NUMBER)ndat);
}

void HMMProblemEloK::setGradB (FitBit *fb){
    if(this->p->block_fitting[2]>0) return;
    NDAT t, ndat = 0;
    NPAR o=0, i;
    NUMBER combined, deriv_logit;
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        ndat += dt->n;
        for(t=0; t<dt->n; t++) {
            o = this->p->dat_obs[ dt->ix[t] ];
            if(o<0)
                continue;
            for(i=0; i<fb->nS; i++) {
                combined = getB(dt,i,o);
                deriv_logit = 1 / safe0num( fb->B[i][o] * (1-fb->B[i][o]) );
                fb->gradB[i][o] -= combined * (1-combined) * deriv_logit * dt->alpha[t][i] * dt->beta[t][i] / safe0num(dt->p_O_param * ((o<0)?1:fb->B[i][o]));
            }
        }
    }
    if( this->p->Cslices>0 ) // penalty
        fb->addL2Penalty(FBV_B, this->p, (NUMBER)ndat);
}


void HMMProblemEloK::toFile(const char *filename) {
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
    NPAR i,j,m;
    fprintf(fid,"K\t%12.10f\n",logit(this->K));
	std::map<NCAT,std::string>::iterator it;
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PI\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%12.10f%s",this->pi[k][i],(i==(this->p->nS-1))?"\n":"\t");
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

void HMMProblemEloK::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure==STRUCTURE_ELOK)
        loglik_rmse[0] += GradientDescent();
    else {
        fprintf(stderr,"Solver specified is not supported.\n");
        exit(1);
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

// extracted from prediction with modifications
NDAT HMMProblemEloK::computeGradientsElo(FitBitElo *fbe){
    struct param *param = fbe->param;
    
    NPAR nO = param->nO, nS = param->nS;
    NCAT nG = param->nG, nK = param->nK;
    NPAR *dat_obs = param->dat_obs;
    NCAT *dat_group = param->dat_group;
    NCAT *dat_skill = param->dat_skill;
    NCAT *dat_skill_stacked = param->dat_skill_stacked; // if multiskill==1, stacked array of all multiskills
    NCAT *dat_skill_rcount = param->dat_skill_rcount;  // if multiskill==1, for each multi-skill row count of skills in it
    NDAT *dat_skill_rix = param->dat_skill_rix;     // if multiskill==1, for each data skill row, index into first element in stacked array
    char f_multiskill = param->multiskill;
    int f_metrics_target_obs = param->metrics_target_obs;
    
    fbe->toZero(FBS_GRAD);
    toValue1D<NUMBER>(fbe->elo_track_g, nG, 0.5);
    toValue1D<NUMBER>(fbe->elo_track_g_t, this->p->N, 0.5);
    toZero1D<NUMBER>(fbe->elo_grad_error_g, nG);
    toZero1D<NCAT>(fbe->elo_count_g, nG);

    NDAT t;
    NCAT g, k;

    NPAR i, j, o, isTarget = 0;
    NUMBER *local_pred = init1D<NUMBER>(nO); // local prediction
    NUMBER *pLe = init1D<NUMBER>(nS);// p(L|evidence);
    NUMBER pLe_denom; // p(L|evidence) denominator
    NUMBER ***group_skill_map = init3D<NUMBER>(nG, nK, nS);//UNBOOST
//   ::boost::numeric::ublas::mapped_matrix<NUMBER*> gsm (nG, nK);//BOOST
    
    NUMBER ll = 0.0;
    NUMBER pcor;
    struct data* dt = new data;
    NDAT count = 0;
    
    for(t=0; t<param->N; t++) {
        o = dat_obs[t];
        g = dat_group[t];
        dt->g = g;
        dt->t = t; // statefulness
        
        // ahead-looking value of tracked Elo rating
        fbe->elo_track_g_t[t] = fbe->elo_track_g[g];
        
        isTarget = this->p->metrics_target_obs == o;
        NCAT *ar;
        int n;
        if(f_multiskill==0) {
            k = dat_skill[t];
            ar = &k;
            n = 1;
        } else {
            k = dat_skill_stacked[ dat_skill_rix[t] ];
            ar = &dat_skill_stacked[ dat_skill_rix[t] ];
            n = dat_skill_rcount[t];
        }
       
        // deal with null skill
        if(ar[0]<0) { // if no skill label
            isTarget = this->null_skill_obs==o;
            ll -= isTarget*safelog(this->null_skill_obs_prob) + (1-isTarget)*safelog(1 - this->null_skill_obs_prob);
            continue;
        }
        // check if {g,k}'s were initialized
        for(int l=0; l<n; l++) {
            k = ar[l];
//          NUMBER *z = gsm(g,k); //BOOST
//          if( z==NULL )//BOOST
            if( group_skill_map[g][k][0]==0)//UNBOOST
            {
                dt->k = k;
//                NUMBER * pLbit = Calloc(NUMBER, nS);//BOOST
                
                for(i=0; i<nS; i++) {
                    group_skill_map[g][k][i] = this->getPI(dt,i);//UNBOOST
//                    pLbit[i] = hmm->getPI(dt,i);//BOOST
                    count++;
                }
//              gsm(g,k) = pLbit; //BOOST
            }// pLo/pL not set
        }// for all skills at this transaction
        
        // produce prediction and copy to result
//        hmm->producePCorrect(group_skill_map, local_pred, ar, n, dt); //UNBOOST
        this->producePCorrect(group_skill_map, local_pred, t); //UNBOOST
//        hmm->producePCorrectBoost(&gsm, local_pred, t); //BOOST
        projectsimplex(local_pred, nO); // addition to make sure there's not side effects
        
        
        // update pL
        for(int l=0; l<n; l++) {
            k = ar[l];
            dt->k = k;
//          NUMBER* pLbit = gsm(g,k); //BOOST
            
            if(o>-1) { // known observations //
                // update p(L)
                pLe_denom = 0.0;
                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
                for(i=0; i<nS; i++)
                    pLe_denom += group_skill_map[g][k][i] * this->getB(dt,i,o);  // TODO: this is local_pred[o]!!!//UNBOOST
//                  pLe_denom += pLbit[i] * hmm->getB(dt,i,o); //BOOST
                for(i=0; i<nS; i++)
                    pLe[i] = group_skill_map[g][k][i] * this->getB(dt,i,o) / safe0num(pLe_denom); //UNBOOST
//                  pLe[i] = pLbit[i] * hmm->getB(dt,i,o) / safe0num(pLe_denom); //BOOST
                // 2. L = (pLe'*A)';
                for(i=0; i<nS; i++)
                    group_skill_map[g][k][i]= 0.0; //UNBOOST
//                  pLbit[i]= 0.0; //BOOST
                for(j=0; j<nS; j++)
                    for(j=0; j<nS; j++)
                        for(i=0; i<nS; i++)
                            group_skill_map[g][k][j] += pLe[i] * this->getA(dt,i,j);//A[i][j]; //UNBOOST
//                          pLbit[j] += pLe[i] * hmm->getA(dt,i,j);//A[i][j]; //BOOST
            } else { // unknown observation
                // 2. L = (pL'*A)';
                for(i=0; i<nS; i++)
                    pLe[i] = group_skill_map[g][k][i]; // copy first; //UNBOOST
//                  pLe[i] = pLbit[i]; // copy first; //BOOST
                for(i=0; i<nS; i++)
                    group_skill_map[g][k][i] = 0.0; // erase old value //UNBOOST
//                  pLbit[i] = 0.0; // erase old value //BOOST
                for(j=0; j<nS; j++)
                    for(i=0; i<nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * this->getA(dt,i,j);//UNBOOST
//               pLbit[j] += pLe[i] * hmm->getA(dt,i,j);//BOOST
            }// observations
            projectsimplex(group_skill_map[g][k], nS); // addition to make sure there's not side effects //UNBOOST
//            projectsimplex(pLbit, nS); // addition to make sure there's not side effects //BOOST
        }
        
        pcor = safe01num(local_pred[f_metrics_target_obs]);
        ll -= safelog(pcor) * isTarget  +  safelog(1-pcor) * (1-isTarget);
        
        // update Elo
        // student
        NUMBER error = isTarget-pcor;
        NCAT n_   = fbe->elo_count_g[g];
        NUMBER sig_ = 1;
        // updating elo_track_g[g] in probabilistic space
        NUMBER s = fbe->elo_track_g[g];
        NUMBER k = fbe->ELO[0];
        fbe->elo_track_g[g] = 1 / (1 + ( (1-s)/s )*pow( (1-k)/k, sig_ * error) );
        fbe->elo_count_g[g]++;
//        fbe->elo_track_g_t[t] = fbe->elo_track_g[g]; // updated before the step
        if(n_>0)
            fbe->gradELO[0] -= error * fbe->elo_grad_error_g[g];
        // \frac { sig\cdot error\left( \frac { 1-x }{ x }  \right) ^{ sig\cdot error } }{ \left( \left( \frac { 1-x }{ x }  \right) ^{ sig\cdot error }+s+1 \right) ^{ 2 }\left( 1-x \right) x }
        fbe->elo_grad_error_g[g] += error * sig_ * pow( (1-k)/k, error * sig_) / (  pow(  pow( (1-k)/k, error * sig_) + (1-s)/s + 1, 2  )*(1-k)*k   )  ;
    } // for all data
    delete(dt);
    free(local_pred);
    free(pLe);
    free3D<NUMBER>(group_skill_map, nG, nK);//UNBOOST
//    gsm.clear();//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator1 it1_t;//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator2 it2_t;//BOOST
//    for (it1_t itgsm1 = gsm.begin1(); itgsm1 != gsm.end1(); itgsm1++)//BOOST
//       for (it2_t itgsm2 = itgsm1.begin(); itgsm2 != itgsm1.end(); itgsm2++)//BOOST
//           free( gsm( itgsm2.index1(), itgsm2.index2() ) );//BOOST
    
    fbe->ll = ll;
    return param->N;
} // computeGradients()

FitResult HMMProblemEloK::GradientDescentBitElo(FitBitElo *fbe) {
    FitResult res;
    FitResult *fr = new FitResult;
    fr->iter = 1;
    fr->pO0  = 0.0;
    fr->pO   = 0.0;
    fr->conv = 0; // converged
    fr->ndat = 0;
    // inital copy parameter values to the t-1 slice
    fbe->copy(FBS_PAR, FBS_PARm1);
    while( !fr->conv && fr->iter<=this->p->maxiter ) {
        fr->ndat = computeGradientsElo(fbe);//a_gradPI, a_gradA, a_gradB);
        
        if(fr->iter==1) {
            fr->pO0 = fbe->ll;
            fr->pOmid = fr->pO0;
        }
        // make step
        if( this->p->solver==METHOD_GD || (fr->iter==1 && this->p->solver==METHOD_CGD)  || (fr->iter==1 && this->p->solver==METHOD_GBB) )
            fr->pO = doLinearStepElo(fbe); // step for linked skill 0
        else if( this->p->solver==METHOD_CGD )
            fr->pO = doConjugateLinearStepElo(fbe);
//        else if( this->p->solver==METHOD_GDL )
//            fr->pO = doLagrangeStep(fbe);
//        else if( this->p->solver==METHOD_GBB )
//            fr->pO = doBarzilaiBorweinStep(fbe);
        // converge?
        fr->conv = fbe->checkConvergence(fr);
        // copy parameter values after we already compared step t-1 with currently computed step t
//        if( this->p->solver==METHOD_GBB )
        fbe->copy(FBS_PARm1, FBS_PARm2); // do this for all in order to capture oscillation, e.g. if new param at t is close to param at t-2 (tolerance)
        fbe->copy(FBS_PAR, FBS_PARm1);
        // report if converged
        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr->conv || fr->iter==this->p->maxiter) )) {
            ;//fr->pO = HMMProblem::getSumLogPOPara(xndat, x_data);
        } else {
            // if Conjugate Gradient
            if (this->p->solver==METHOD_CGD) {
                if( fr->iter==1 ) {
                    fbe->copy(FBS_GRAD, FBS_DIRm1); // gradient is not direction, it's negative direction, hence, need to negate it
                    fbe->negate(FBS_DIRm1);
                }
                else fbe->copy(FBS_DIR,  FBS_DIRm1);
                fbe->copy(FBS_GRAD, FBS_GRADm1);
            }
            // if Barzilai Borwein Gradient Method
            if (this->p->solver==METHOD_GBB) {
                fbe->copy(FBS_GRAD, FBS_GRADm1);
            }
        }
        fr->iter ++;
        fr->pOmid = fr->pO;
    }// single skill loop
    // cleanup
    fr->iter--;
    // this is the case when 2 steps lead to step back to the initial value, i.e. oscillation
    if(fr->iter==2 && fr->conv && fbe->checkConvergenceSingle(fr) ) {
        // decrease iteration counter to 1
        fr->iter--;
    }
    
    res.iter = fr->iter;
    res.pO0  = fr->pO0;
    res.pO   = fr->pO;
    res.conv = fr->conv;
    res.ndat = fr->ndat;
    delete fr;
    return res;
}

NUMBER HMMProblemEloK::doLinearStepElo(FitBitElo *fbe) {
    NPAR i, nELO = fbe->nELO; //nS = fbe->nS, ;
    fbe->doLog10ScaleGentle(FBS_GRAD);
    
//    NCAT xndat = fb->xndat;
//    struct data **x_data = fb->x_data;
    fbe->init(FBS_PARcopy);
    fbe->init(FBS_GRADcopy);
    
    NUMBER e = this->p->ArmijoSeed/2; // step seed as twice as short for the Elo
    bool compliesArmijo = false;
    bool compliesWolfe2 = false; // second wolfe condition is turned on, if satisfied - honored, if not, just the 1st is used
    NUMBER e_Armijo = -1; // e (step size) at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER f_xkplus1_Armijo = 0; // f_xkplus1_Armijo at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER f_xk = fbe->ll;
    NUMBER f_xkplus1 = 0;
    
    fbe->copy(FBS_PAR, FBS_PARcopy); // save original parameters
    fbe->copy(FBS_GRAD, FBS_GRADcopy); // save original gradient
    // compute p_k * -p_k
    NUMBER p_k_by_neg_p_k = 0;
    for(i=0; i<nELO; i++)
        p_k_by_neg_p_k -= fbe->gradELO[i]*fbe->gradELO[i];
    
    int iter = 0; // limit iter steps to 20, via ArmijoMinStep (now 10)
    while( !(compliesArmijo && compliesWolfe2) && e > this->p->ArmijoMinStep) {
        // update
        for(i=0; i<nELO; i++) {
            fbe->ELO[i] = fbe->ELOcopy[i] - e * fbe->gradELOcopy[i];
        }
        
        // recompute ll and gradients
        computeGradientsElo(fbe);
        fbe->doLog10ScaleGentle(FBS_GRAD);
        // compute f(x_{k+1})
        f_xkplus1 = fbe->ll;
        // compute Armijo compliance
        compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
        
        // compute Wolfe 2
        NUMBER p_k_by_neg_p_kp1 = 0;
        // NDAT ndat = computeGradients(fb);
        for(i=0; i<nELO; i++)
            p_k_by_neg_p_kp1 -= fbe->gradELOcopy[i]*fbe->gradELO[i];
        
        compliesWolfe2 = (p_k_by_neg_p_kp1 >= this->p->ArmijoC2 * p_k_by_neg_p_k);
        
        if( compliesArmijo && e_Armijo==-1 ){
            e_Armijo = e; // save the first time Armijo is statisfied, in case we'll roll back to it when Wolfe 2 is finnaly not satisfied
            f_xkplus1_Armijo = f_xkplus1;
        }
        
        e /= (compliesArmijo && compliesWolfe2)?1:this->p->ArmijoReduceFactor;
        iter++;
    } // armijo loop
    
//    fb->copy(FBS_GRADcopy, FBS_GRAD); // return the original gradient in its place // unnecessary
    
    if(!compliesArmijo) { // we couldn't step away from current, copy the inital point back
        e = 0;
        fbe->copy(FBS_PARcopy, FBS_PAR);
        f_xkplus1 = f_xk;
    } else if(compliesArmijo && !compliesWolfe2) { // we couldn't step away from current, copy the inital point back
        e = e_Armijo; // return the first Armijo-compliant e
        f_xkplus1 = f_xkplus1_Armijo; // return the first Armijo-compliant f_xkplus1
        // create new versions of FBS_PAR using e_Armijo as a step
        fbe->copy(FBS_GRADcopy, FBS_GRAD); // return the old gradient
        // update
        for(i=0; i<nELO; i++)
            fbe->ELO[i] = fbe->ELOcopy[i] - e * fbe->gradELOcopy[i];
        // ^^^^^ end of create new versions of FBS_PAR using e_Armijo as a step
    }
    fbe->destroy(FBS_PARcopy);
    fbe->destroy(FBS_GRADcopy);
    return f_xkplus1;
} // doLinearStepElo

NUMBER HMMProblemEloK::doConjugateLinearStepElo(FitBitElo *fbe) {
    NPAR i=0, nELO = fbe->nELO;
    // first scale down gradients
    fbe->doLog10ScaleGentle(FBS_GRAD);
    
    // compute beta_gradient_direction
    NUMBER beta_grad = 0, beta_grad_num = 0, beta_grad_den = 0;
    
    switch (this->p->solver_setting) {
        case 1: // Fletcher-Reeves
            for(i=0; i<nELO; i++) {
                beta_grad_num += fbe->gradELO  [i]*fbe->gradELO  [i];
                beta_grad_den += fbe->gradELOm1[i]*fbe->gradELOm1[i];
            }
            break;
        case 2: // Polak–Ribiere
            for(i=0; i<nELO; i++) {
                beta_grad_num += -fbe->gradELO[i]*(-fbe->gradELO[i] + fbe->gradELOm1[i]);
                beta_grad_den +=  fbe->gradELOm1[i]*fbe->gradELOm1[i];
            }
            break;
        case 3: // Hestenes-Stiefel
            for(i=0; i<nELO; i++) {
                beta_grad_num +=  fbe->gradELO[i]* (-fbe->gradELO[i] + fbe->gradELOm1[i]); // no -, since neg gradient and - is +
                beta_grad_den +=  fbe->dirELOm1[i]*(-fbe->gradELO[i] + fbe->gradELOm1[i]);
            }
            break;
        case 4: // Dai-Yuan
            for(i=0; i<nELO; i++) {
                beta_grad_num += -fbe->gradELO [i]*fbe->gradELO  [i];
                beta_grad_den +=  fbe->dirELOm1[i]*(-fbe->gradELO[i] + fbe->gradELOm1[i]);
            }
            break;
        default:
            fprintf(stderr,"Wrong stepping algorithm (%d) specified for Conjugate Gradient Descent\n",this->p->solver_setting);
            break;
    }
    beta_grad = beta_grad_num / safe0num(beta_grad_den);
    beta_grad = (beta_grad>=0)?beta_grad:0;
    // compute new direction
    fbe->toZero(FBS_DIR);
    for(i=0; i<nELO; i++)
        fbe->dirELO[i] = -fbe->gradELO[i] + beta_grad * fbe->dirELOm1[i];
    // scale direction
    fbe->doLog10ScaleGentle(FBS_DIR);
    
    NUMBER e = this->p->ArmijoSeed/2; // step seed as twice as short for the Elo
    bool compliesArmijo = false;
    bool compliesWolfe2 = false; // second wolfe condition is turned on, if satisfied - honored, if not, just the 1st is used
    NUMBER e_Armijo = -1; // e (step size) at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER f_xkplus1_Armijo = 0; // f_xkplus1_Armijo at which Armijo (Wolfe 1) criterion is satisfied, 'cos both criterions are sometimes not met.
    NUMBER f_xk = fbe->ll;
    NUMBER f_xkplus1 = 0;
    
    fbe->init(FBS_PARcopy);
    fbe->init(FBS_GRADcopy);
    
    fbe->copy(FBS_PAR, FBS_PARcopy); // copy parameter
    fbe->copy(FBS_GRAD, FBS_GRADcopy); // copy initial gradient
    // compute p_k * -p_k >>>> now current gradient by current direction
    NUMBER p_k_by_neg_p_k = 0;
    for(i=0; i<nELO; i++)
        p_k_by_neg_p_k += fbe->gradELO[i]*fbe->dirELO[i];
    int iter = 0; // limit iter steps to 20, actually now 10, via ArmijoMinStep
    while( !compliesArmijo && e > this->p->ArmijoMinStep) {
        // update
        for(i=0; i<nELO; i++)
            fbe->ELO[i] = fbe->ELOcopy[i] + e * fbe->dirELO[i];
        // recompute alpha and p(O|param)
        computeGradientsElo(fbe);
        fbe->doLog10ScaleGentle(FBS_GRAD);
        // compute f(x_{k+1})
        f_xkplus1 = fbe->ll;
        // compute Armijo compliance
        compliesArmijo = (f_xkplus1 <= (f_xk + (this->p->ArmijoC1 * e * p_k_by_neg_p_k)));
        // compute Wolfe 2
        NUMBER p_k_by_neg_p_kp1 = 0;
        // gradients already computed
        //computeGradients(fb);
        for(i=0; i<nELO; i++)
            p_k_by_neg_p_kp1 += fbe->dirELO[i]*fbe->gradELO[i];
        compliesWolfe2 = (p_k_by_neg_p_kp1 >= this->p->ArmijoC2 * p_k_by_neg_p_k);
        
        if( compliesArmijo && e_Armijo==-1 ){
            e_Armijo = e; // save the first time Armijo is statisfied, in case we'll roll back to it when Wolfe 2 is finnaly not satisfied
            f_xkplus1_Armijo = f_xkplus1;
        }
        
        e /= (compliesArmijo && compliesWolfe2)?1:this->p->ArmijoReduceFactor;
        iter++;
    } // armijo loop
    
    fbe->copy(FBS_GRADcopy, FBS_GRAD); // return the original gradient in its place
    
    if(!compliesArmijo) { // we couldn't step away from current, copy the inital point back
        e = 0;
        fbe->copy(FBS_PARcopy, FBS_PAR);
        f_xkplus1 = f_xk;
    } else if(compliesArmijo && !compliesWolfe2) { // we couldn't step away from current, copy the inital point back
        e = e_Armijo;
        f_xkplus1 = f_xkplus1_Armijo;
        // create new versions of FBS_PAR using e_Armijo as a step
        // update
        for(i=0; i<nELO; i++)
            fbe->ELO[i] = fbe->ELOcopy[i] + e * fbe->dirELO[i];
    }
    fbe->destroy(FBS_PARcopy);
    fbe->destroy(FBS_GRADcopy);
    return f_xkplus1;
} // doConjugateLinearStepElo

NUMBER HMMProblemEloK::GradientDescent() {
	NCAT k, /*ki, gi, nX, */x;
    NCAT nK = this->p->nK, nG = this->p->nG;
    
    // reset values
    if(this->elo_track_g == NULL)
        this->elo_track_g = init1D<NUMBER>(nG);
    toValue1D<NUMBER>(this->elo_track_g, nG, 0.5);
    if(this->elo_count_g == NULL)
        this->elo_count_g = init1D<NCAT>(nG); // init1DNumber(this->nS);
    toZero1D<NCAT>(this->elo_count_g, nG);
    
    FitResult fr;

    //
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill>0) {
        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
        fb->link( HMMProblem::getPI(0), HMMProblem::getA(0), HMMProblem::getB(0), this->p->nSeq, this->p->k_data);// link skill 0 (we'll copy fit parameters to others
        if(this->p->block_fitting[0]!=0) fb->pi = NULL;
        if(this->p->block_fitting[1]!=0) fb->A  = NULL;
        if(this->p->block_fitting[2]!=0) fb->B  = NULL;

        fb->init(FBS_PARm1);
        fb->init(FBS_PARm2);
        fb->init(FBS_GRAD);
        if(this->p->solver==METHOD_CGD) {
            fb->init(FBS_GRADm1);
            fb->init(FBS_DIRm1);
        }

        NCAT* original_ks = Calloc(NCAT, (size_t)this->p->nSeq);
        for(x=0; x<this->p->nSeq; x++) { original_ks[x] = this->p->all_data[x].k; this->p->all_data[x].k = 0; } // save progonal k's
        FitResult fr = GradientDescentBit(fb);
        for(x=0; x<this->p->nSeq; x++) { this->p->all_data[x].k = original_ks[x]; } // restore original k's
        free(original_ks);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
        delete fb;
    }
	//
	// Main fit
	//
    if( this->p->single_skill!=2 ) { // if not "force single skill"
        NCAT first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
        NCAT iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
        NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK); // counting concecutive number of convergences for skills
        NCAT iter_qual_k = 0; // counting concecutive number of convergences for Elo K
        NCAT* iter_fail_skill = Calloc(NCAT, (size_t)nK); // counting concecutive number of failures to converge for skills
        NCAT iter_fail_k = 0; // counting concecutive number of failures to converge for Elo K
        int skip_k = 0;
        // utilize fitting larger data first
        
        int i = 0; // count runs
//        int parallel_now = this->p->parallel==1; //PAR
//        #pragma omp parallel if(parallel_now) shared(iter_qual_group,iter_qual_skill,iter_fail_skill,iter_fail_group)//PAR
//        {//PAR
            while(skip_k<nK) {
                //
                // Skills first
                //
                if(skip_k<nK) {
//                    #pragma omp for schedule(dynamic) //PAR
                    for(k=0; k<nK; k++) { // for all A,B-by-skill
                        if(iter_qual_skill[k]==iterations_to_qualify || iter_fail_skill[k]==iterations_to_qualify)
                            continue;
                        FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol, this->p->tol_mode);
                        // link and fit
                        fb->link( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), this->p->k_numg[k], this->p->k_g_data[k]);// link skill 0 (we'll copy fit parameters to others
                        if(this->p->block_fitting[0]!=0) fb->pi = NULL;
                        if(this->p->block_fitting[1]!=0) fb->A  = NULL;
                        if(this->p->block_fitting[2]!=0) fb->B  = NULL;

                        fb->Cslice = 0;
                        fb->init(FBS_PARm1);
                        fb->init(FBS_PARm2);
                        fb->init(FBS_GRAD);
                        if(this->p->solver==METHOD_CGD) {
                            fb->init(FBS_GRADm1);
                            fb->init(FBS_DIRm1);
                        }
                        fr = GradientDescentBit(fb);
                        
                        // count number of concecutive failures
                        if(fr.iter==this->p->maxiter) {
                            iter_fail_skill[k]++;
                        } else {
                            iter_fail_skill[k]=0;
                        }

                        // decide on convergence
                        if(i>=first_iteration_qualify || fb->xndat==0) {
                            if(fr.iter==1 /*e<=this->p->tol*/ || fb->xndat==0) { // converged quick, or don't care (others all converged
                                iter_qual_skill[k]++;
                                if(fb->xndat==0) {
                                    iter_qual_skill[k]=iterations_to_qualify;
                                    fr.conv = 1;
                                }
                                if(iter_qual_skill[k]==iterations_to_qualify) {// criterion met, or don't care (others all converged)
//                                   #pragma omp critical(update_skip_k)//PAR
//                                   {//PAR
                                    skip_k++;
//                                   }//PAR
                                    if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                        printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                    }
                                }
                            } else
                                iter_qual_skill[k]=0;
                            // failed too many times in a row
                            if( iter_fail_skill[k]==iterations_to_qualify ) {
//                                #pragma omp critical(update_skip_k)//PAR
//                                {//PAR
                                skip_k++;
//                                }//PAR
                                if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                    printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                                }
                            }
                        } // decide on convergence
                        delete fb;
                    } // for all skills
                }
                //
                // Elo K second
                //
                if(iter_qual_k==iterations_to_qualify || iter_fail_k==iterations_to_qualify)
                    continue;
                FitBitElo *fbe = new FitBitElo(this->p->nG, this->p->nELO, this->p->tol, this->p->tol_mode);
                // link
                fbe->link(&this->K, this->elo_track_g, this->elo_count_g, this->elo_track_g_t, this->p);

                fbe->init(FBS_PARm1);
                fbe->init(FBS_PARm2);
                fbe->init(FBS_GRAD);
                if(this->p->solver==METHOD_CGD) {
                    fbe->init(FBS_GRADm1);
                    fbe->init(FBS_DIRm1);
                }
                fr = GradientDescentBitElo(fbe);
                
                // count number of concecutive failures
                if(fr.iter==this->p->maxiter) {
                    iter_fail_k++;
                } else {
                    iter_fail_k=0;
                }

                // decide on convergence
                if(i>=first_iteration_qualify ) { //can qualify or  student had no skill labelled rows
                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==nK ) { // converged quick, or don't care (others all converged), or
                        iter_qual_k++;
                        if(iter_qual_k==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                            if(skip_k==nK) iter_qual_k=iterations_to_qualify; // K not changing anymore
                            if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                printf("run %2d doneELO               iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,fr.iter,fr.pO0,fr.pO,fr.conv);
                            }
                        }
                    }
                    else
                        iter_qual_k=0;
                    // failed too many times in a row
                    if( iter_fail_k==iterations_to_qualify ) {
                        if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                            printf("run %2d doneELO iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n",i,fr.iter,fr.pO0,fr.pO,fr.conv);
                        }
                    }
                } // decide on convergence
                delete fbe;
                //
                // done with fitting K
                //
//                #pragma omp single//PAR
//                {//PAR
                    i++;
//                }//PAR
            }
//        }//PAR
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_fail_skill != NULL ) free(iter_fail_skill);
        
    } // if not "force single skill
    
    // cleanup
    toValue1D<NUMBER>(this->elo_track_g, nG, 0.5);
    toValue1D<NUMBER>(this->elo_track_g_t, this->p->N, 0.5);
    toZero1D<NCAT>(this->elo_count_g, nG);
    // compute loglik
    return fr.pO;//getSumLogPOPara(this->p->nSeq, this->p->k_data);
}

void HMMProblemEloK::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
	NPAR i,j,m;
	NCAT k = 0, idxk = 0;
	string s;
    char col[2048];
    std::map<std::string,NCAT>::iterator it;
    //
    readNullObsRatio(fid, param, line_no);
    //
    // init param
    //
    if(overwrite) {
        this->p->map_group_fwd = new map<string,NCAT>();
        this->p->map_group_bwd = new map<NCAT,string>();
        this->p->map_skill_fwd = new map<string,NCAT>();
        this->p->map_skill_bwd = new map<NCAT,string>();
    }
	//
	// read Elo K
    fscanf(fid,"K\t");
    fscanf(fid,"%[^\n]\n", col); // last one;
    this->K = sigmoid(atof(col)); // keep all ELO parameters in probabilistic range
    (*line_no)++;
	//
    //
    // read skills
    //
	for(k=0; k<this->p->nK; k++) {
		// read skill label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;
        if(overwrite) {
            this->p->map_skill_fwd->insert(pair<string,NCAT>(s, (NCAT)this->p->map_skill_fwd->size()));
            this->p->map_skill_bwd->insert(pair<NCAT,string>((NCAT)this->p->map_skill_bwd->size(), s));
            idxk = k;
        } else {
            it = this->p->map_skill_fwd->find(s);
            if( it==this->p->map_skill_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                continue; // skip this iteration
            }
            else
                idxk =it->second;
        }
        
        // read PI
        fscanf(fid,"PIk\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->pi[idxk][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->pi[idxk][i] = atof(col);
        (*line_no)++;
		// read A
        fscanf(fid,"Ak\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++) {
                if(i==(this->p->nS-1) && j==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->A[idxk][i][j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->A[idxk][i][j] = atof(col);
                }
			}
        (*line_no)++;
		// read B
        fscanf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nO-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->B[idxk][i][m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->B[idxk][i][m] = atof(col);
                }
			}
        (*line_no)++;
	} // for all k
}

// reserved for any per-row post-processing after producePCorrect, Elo uses it to update its tracked values
void HMMProblemEloK::postprocesslocal(NUMBER target, NUMBER estimate, NDAT t) {
    // update Elo ratings
    NCAT g = this->p->dat_group[t];
    // student
    NUMBER error = target-estimate;
    NUMBER sig_ = 1;
    NUMBER s = this->elo_track_g[g];
    NUMBER k = this->K;
    this->elo_track_g[g] = 1 / (1 + ( (1-s)/s )*pow( (1-k)/k, sig_ * error) );
    return;
}

// reserved for any pre-processing before prediction, Elo uses it to update its tracked values
void HMMProblemEloK::preprocessglobal() {
//    // resetting tracked variables etc.
//    toValue1D<NUMBER>(this->elo_track_g, this->p->nG, 0.5);
//    toValue1D<NUMBER>(this->elo_track_g_t, this->p->N, 0.5);
//    toZero1D<NCAT>(this->elo_count_g, this->p->nG);
//    return;

    // create
    FitBitElo *fbe = new FitBitElo(this->p->nG, this->p->nELO, this->p->tol, this->p->tol_mode);
    // link
    fbe->link(&this->K, this->elo_track_g, this->elo_count_g, this->elo_track_g_t, this->p);
    // init
    fbe->init(FBS_GRAD);
    // run
    computeGradientsElo(fbe);
    // recycle
    delete fbe;
    return;
}

//// nothing for BKT, computes Elo values for Elo-based
//void HMMProblemEloK::preprocess_computeAlphaAndPOParam() {
//    // create
//    FitBitElo *fbe = new FitBitElo(this->p->nG, this->p->nELO, this->p->tol, this->p->tol_mode);
//    // link
//    fbe->link(&this->K, this->elo_track_g, this->elo_count_g, this->elo_track_g_t, this->p);
//    // init
//    fbe->init(FBS_GRAD);
//    // run
//    computeGradientsElo(fbe);
//    // recycle
//    delete fbe;
//    return;
//}
