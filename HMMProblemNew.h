/*
 *  HMMProblem.h
 *  HMM
 *
 *  Created by Mikhail Yudelson on 5/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "StripedArray.h"
#include "utils.h"
#include "FitBit.h"
#include <stdio.h>
#include <string>

//#include "liblinear/linear.h"

//#include <Accelerate/Accelerate.h>

#ifndef _HMMPROBLEM_H
#define _HMMPROBLEM_H

class HMMProblem {
public:
	HMMProblem();
	HMMProblem(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblem();
	virtual NUMBER** getPI();
	virtual NUMBER*** getA();
	virtual NUMBER*** getB();
	virtual NUMBER* getPI(NCAT k);
	virtual NUMBER** getA(NCAT k);
	virtual NUMBER** getB(NCAT k);
	NUMBER* getLbPI();
	NUMBER** getLbA();
	NUMBER** getLbB();
	NUMBER* getUbPI();
	NUMBER** getUbA();
	NUMBER** getUbB();
    // getters for computing alpha, beta, gamma, grads, they use t - position in dataset to combine multiple KCs of Group
    virtual NUMBER getPI(NDAT t, NPAR i);  // also works as get pL :)
    virtual NUMBER getA (NDAT t, NPAR i, NPAR j);
    virtual NUMBER getB (NDAT t, NPAR i, NPAR m);
    // getters for computing gradients of alpha, beta, gamma - for all skills in one big sweep
    virtual void setGradPI(struct data* dt, FitBit *fb);
    virtual void setGradA (struct data* dt, FitBit *fb);
    virtual void setGradB (struct data* dt, FitBit *fb);
    void setGradPI(NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb);
    void setGradA(NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb);
    void setGradB(NPAR o, NDAT t, NDAT ix, NUMBER p_O_param, bool is_multiskill, FitBit *fb);
	virtual void toFile(const char *filename);
	NUMBER getSumLogPOPara(NCAT2 ix1, NCAT2 n, struct data *data); // generic per k/g-slice
	NUMBER getSumLogPOPara(NCAT k); // for particular skill
	bool hasNon01Constraints();
    NUMBER getLogLik(); // get log likelihood of the fitted model
    NCAT getNparams(); // get log likelihood of the fitted model
    NUMBER getNullSkillObs(NPAR m); // get log likelihood of the fitted model
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // compute metrics
    void computeMetrics(NUMBER* metrics);
    // predicting
//    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
    void predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill, bool only_unlabeled);
    virtual void readModelBody(FILE *fid, NDAT *line_no);
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
	NUMBER** PI; // initial state probabilities
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
    
    // used for fitting
//    NUMBER *** pL;   // pLearned, first by size of data (N), then by num states (nS) and then by skill (for convenience)
    NUMBER ** pL;   // pLearned, N x nS combined for all involved KCs
    NUMBER ** alpha; // forward variable(s) NN x nS
    NUMBER ** beta;  // backward variables NN x nS
    NUMBER * c;      // alpha scaling factor

	//
	// Derived
	//
//	NUMBER* gradPI; // gradient of initial state probabilities
//	NUMBER** gradA; // gradient of transition matrix
//	NUMBER** gradB; // gradient of observation matrix
	
	virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
    
	void initAlpha();//NCAT xndat, struct data** x_data); // generic
	void initBeta();//(NCAT xndat, struct data** x_data); // generic
	void initXi();//(NCAT xndat, struct data** x_data); // generic
	void initGamma();//(NCAT xndat, struct data** x_data); // generic
    
	void computeAlphaEtAl0(); // compute Alpha, p(O|param), pL too, it is global fastest version w/o subfunction
	void computeAlphaEtAl0(NCAT target_k, bool doZeroCount); // compute Alpha, p(O|param), pL too, for a particular k skill (nK+1st is null skill)
	void computeAlphaEtAl(); // compute Alpha, p(O|param), pL too, it is global fastest version w/o subfunction
	void computeAlphaEtAl(NCAT target_k, bool doZeroCount); // compute Alpha, p(O|param), pL too, for a particular k skill (nK+1st is null skill)
    
	void computeBeta(); //NCAT xndat, struct data** x_data); // it is global
	void computeBeta(NCAT target_k); // compute Beta, p(O|param), pL too, for a particular k skill (nK+1st is null skill)
    
	void computeGamma();//(NCAT xndat, struct data** x_data);
	void computeXi();//(NCAT xndat, struct data** x_data);
    
    void FitNullSkill(NUMBER* loglik_rmse, bool keep_SE); // get loglik and RMSE
//    // helpers
//    void init3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
//    void toZero3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
//    void free3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS);
    void cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO);
    // predicting
	virtual void computeLogLikRMSENullSkill(NUMBER* loglik_rmse, bool keep_SE);
	virtual void computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE);//, NCAT xndat, struct data** x_data);
//	virtual void computeGradients(NCAT2 ix1, NCAT2 n, struct data* data, FitBit *fb, NCAT target_k);
	virtual void computeGradients(FitBit **fbs); // global gradient compute
	virtual void computeGradients(FitBit **fbs, NCAT target_k); // just for one skill gradient compute
    virtual NUMBER doLinearStep(NCAT2 ix1, NCAT2 n, struct data* data, NCAT k_target, FitBit *fb, NCAT copy);//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
    virtual NUMBER doLinearStep(NCAT k, FitBit *fb);
    NUMBER doConjugateLinearStep(NCAT2 ix1, NCAT2 n, struct data* data, NCAT k_target, FitBit *fb, NCAT copy);//NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
    NUMBER doConjugateLinearStep(NCAT k, FitBit *fb);
    NUMBER doBarzalaiBorweinStep(NCAT ndat, struct data* data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
//    bool checkConvergence(NUMBER* PI, NUMBER** A, NUMBER** B, NUMBER* PI_m1, NUMBER** A_m1, NUMBER** B_m1, bool flags[3]);
    // bridge to Liblinear - all objects are liblinear objects
//    void createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space);
//    void createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space, NCAT k); // multiple KC's separately
//    void recycleLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space);
//    void GradientDescent1Skill(FitBit *fb); // return -LL for the model
    FitResult GradientDescentBit(NCAT2 ix1, NCAT2 n, struct data* data, NCAT k_target, NPAR kg_flag, FitBit *fb, bool is1SkillForAll); // for 1 skill or 1 group, all 1 skill for all data
    virtual NCAT GradientDescentBit(FitResult *frs, FitBit **fbs); // for 1 skill or 1 group, all 1 skill for all data
//    void GradientConjugateDescent1Skill(FitBit *fb); // return -LL for the model
    virtual NUMBER GradientDescent0(); // return -LL for the model - old
    virtual NUMBER GradientDescent(); // return -LL for the model - new version
    void readNullObsRatio(FILE *fid, NDAT *line_no);
private:
	bool checkPIABConstraints(NUMBER* a_PI, NUMBER** a_A, NUMBER** a_B); // all constraints, inc row sums
    // fitting methods - helpers (hidden)
    // fitting methods (hidden)
//    NUMBER ConjugateGradientDescent(NPAR kg_flag); // return -LL for the model
    NUMBER BaumWelchSkill();
    void doBaumWelchStep(NCAT ndat, struct data* data, FitBit *fb);//, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B);
//    void predictMetricsGroup(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill);
    // write model
	void toFileSkill(const char *filename);
	void toFileGroup(const char *filename);
};

#endif