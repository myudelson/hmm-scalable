/*
 
 Copyright (c) 2012-2015, Michael (Mikhail) Yudelson
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

#include "utils.h"
#include "FitBit.h"
#include "StripedArray.h"

//#include <boost/numeric/ublas/matrix_sparse.hpp>//BOOST
//#include <boost/numeric/ublas/io.hpp>//BOOST

#ifndef _HMMPROBLEM_H
#define _HMMPROBLEM_H

class HMMProblem {
public:
	HMMProblem();
	HMMProblem(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblem();
    NUMBER** getPI();
	NUMBER*** getA();
	NUMBER*** getB();
	NUMBER* getPI(NCAT k);
	NUMBER** getA(NCAT k);
	NUMBER** getB(NCAT k);
	NUMBER* getLbPI();
	NUMBER** getLbA();
	NUMBER** getLbB();
	NUMBER* getUbPI();
	NUMBER** getUbA();
	NUMBER** getUbB();
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j);
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m);
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(FitBit *fb);
    virtual void setGradA (FitBit *fb);
    virtual void setGradB (FitBit *fb);
	virtual void toFile(const char *filename);
	NUMBER getSumLogPOPara(NCAT xndat, struct data **x_data); // generic per k/g-slice
	bool hasNon01Constraints();
    NUMBER getLogLik(); // get log likelihood of the fitted model
    NCAT getNparams(); // get log likelihood of the fitted model
    NUMBER getNullSkillObs(NPAR m); // get log likelihood of the fitted model
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // predicting
//    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NDAT t);
//    virtual void producePCorrectBoost(boost::numeric::ublas::mapped_matrix<NUMBER*> *group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);//BOOST
//	virtual void producePDObs(NUMBER*** group_skill_map, NUMBER* local_pred, struct data* dt); // probability distribution of observations, single skill label
//	virtual void producePDObsBoost(boost::numeric::ublas::mapped_matrix<NUMBER*> *group_skill_map, NUMBER* local_pred, struct data* dt);//BOOST
    static void predict(NUMBER* metrics, const char *filename, NPAR* dat_obs, NCAT *dat_group, NCAT *dat_skill, NCAT *dat_skill_stacked, NCAT *dat_skill_rcount, NDAT *dat_skill_rix, HMMProblem **hmms, NPAR nhmms, NPAR *hmm_idx);
	
	void predictNEW(NUMBER* metrics, const char *filename, NPAR* dat_obs, NCAT *dat_group, NCAT *dat_skill, NCAT *dat_skill_stacked, NCAT *dat_skill_rcount, NDAT *dat_skill_rix);
    void readModel(const char *filename, bool overwrite);
    virtual void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
//    virtual void reorderSequences(NDAT *newnK, NDAT *newnG, bool sort); /*place larger skill and group sequences closer to the beginning*/
protected:
	//
	// Givens
	//
    NCAT n_params; // number of model params
    NCAT sizes[3]; // sizes of arrays of PI,A,B params
    NUMBER *null_obs_ratio;
    NUMBER neg_log_lik; // negative log-likelihood
    NUMBER null_skill_obs; // if null skills are present, what's the default obs to predict
    NUMBER null_skill_obs_prob; // if null skills are present, what's the default obs probability to predict
	NUMBER** pi; // initial state probabilities
	NUMBER*** A; // transition matrix
	NUMBER*** B; // observation matrix
	NUMBER* lbPI; // lower boundary initial state probabilities
	NUMBER** lbA; // lower boundary transition matrix
	NUMBER** lbB; // lower boundary observation matrix
	NUMBER* ubPI; // upper boundary initial state probabilities
	NUMBER** ubA; // upper boundary transition matrix
	NUMBER** ubB; // upper boundary observation matrix
	bool non01constraints; // whether there are lower or upper boundaries different from 0,1 respectively
	struct param *p; // data and params
	//
	// Derived
	//
	virtual void init(struct param *param); // non-fit specific initialization
    virtual void destroy(); // non-fit specific descruction
	void initAlpha(NCAT xndat, struct data** x_data); // generic
	void initXiGamma(NCAT xndat, struct data** x_data); // generic
	void initBeta(NCAT xndat, struct data** x_data); // generic
	NDAT computeAlphaAndPOParam(NCAT xndat, struct data** x_data);
	void computeBeta(NCAT xndat, struct data** x_data);
	void computeXiGamma(NCAT xndat, struct data** x_data);
    void FitNullSkill(NUMBER* loglik_rmse, bool keep_SE); // get loglik and RMSE
    // helpers
    void init3Params(NUMBER* &pi, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
    void toZero3Params(NUMBER* &pi, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
    void free3Params(NUMBER* &pi, NUMBER** &A, NUMBER** &B, NPAR nS);
    void cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO);

    // predicting
	virtual NDAT computeGradients(FitBit *fb);
    virtual NUMBER doLinearStep(FitBit *fb);
    virtual NUMBER doLagrangeStep(FitBit *fb);
    NUMBER doConjugateLinearStep(FitBit *fb);
    NUMBER doBaumWelchStep(FitBit *fb);
    FitResult GradientDescentBit(FitBit *fb); // for 1 skill or 1 group, all 1 skill for all data
    FitResult BaumWelchBit(FitBit *fb);
	// predicting big
    virtual NDAT computeGradientsBig(FitBit **fbs, NCAT nfbs);// global, for all
    NUMBER doLinearStepBig(FitBit **fbs, NCAT nfbs); //
    NUMBER doConjugateLinearStepBig(FitBit **fbs, NCAT nfbs); //
    FitResult GradientDescentBitBig(FitBit **fbs, NCAT nfbs); //
    bool checkConvergenceBig0(FitBit** fbs, NCAT nfbs, NUMBER tol, NUMBER *criterion);
    bool checkConvergenceBig(FitBit** fbs, NCAT nfbs, NUMBER tol, NUMBER *criterion);
    
    NUMBER doBarzilaiBorweinStep(FitBit *fb);
//    NUMBER doBarzilaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
    virtual NUMBER GradientDescent(); // return -LL for the model
    NUMBER BaumWelch(); // return -LL for the model
    void readNullObsRatio(FILE *fid, struct param* param, NDAT *line_no);
	bool checkPIABConstraints(NUMBER* a_PI, NUMBER** a_A, NUMBER** a_B); // all constraints, inc row sums
private:
    // write model
	void toFileSkill(const char *filename);
	void toFileGroup(const char *filename);
};

#endif
