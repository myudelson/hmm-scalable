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
#include <stdio.h>
#include <string>

#include "liblinear/linear.h"

//#include <Accelerate/Accelerate.h>

#ifndef _HMMPROBLEM_H
#define _HMMPROBLEM_H

class HMMProblem {
public:
	HMMProblem();
	HMMProblem(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblem();
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
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j);
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m);
	virtual void toFile(const char *filename);
//	static NUMBER getSumNegPOPara(NCAT xndat, struct data **x_data, NUMBER (*f)(NUMBER)); // generic per k/g-slice
	static NUMBER getSumLogPOPara(NCAT xndat, struct data **x_data); // generic per k/g-slice
	bool hasNon01Constraints();
    NUMBER getLogLik(); // get log likelihood of the fitted model
    NCAT getNparams(); // get log likelihood of the fitted model
    NUMBER getNullSkillObs(NPAR m); // get log likelihood of the fitted model
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // compute metrics
    void computeMetrics(NUMBER* metrics);
    // predicting
    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
    void predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill);
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

	//
	// Derived
	//
	NUMBER* gradPI; // gradient of initial state probabilities
	NUMBER** gradA; // gradient of transition matrix
	NUMBER** gradB; // gradient of observation matrix
	
	virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific destruction
	static void initAlpha(NCAT xndat, struct data** x_data, NPAR a_nS); // generic
	static void initXi(NCAT xndat, struct data** x_data, NPAR a_nS); // generic
	static void initGamma(NCAT xndat, struct data** x_data, NPAR a_nS); // generic
	static void initBeta(NCAT xndat, struct data** x_data, NPAR a_nS); // generic
//	static void computeAlphaAndPOParam(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
	void computeAlphaAndPOParam(NCAT xndat, struct data** x_data);
//	void computeAlphaAndPOParamG(NCAT g);
//	virtual void computeBeta(NCAT xndat, struct data** x_data, NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
	void computeBeta(NCAT xndat, struct data** x_data);
	void computeGamma(NCAT xndat, struct data** x_data);
	void computeXi(NCAT xndat, struct data** x_data);
    void FitNullSkill(NUMBER* loglik_rmse, bool keep_SE); // get loglik and RMSE
    // predicting
	virtual void computeLogLikRMSENullSkill(NUMBER* loglik_rmse, bool keep_SE);
	virtual void computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data);
	void computeGradients(NCAT xndat, struct data** x_data, NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB, struct param* param);
    NUMBER doLinearStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
                        NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
    NUMBER doConjugateLinearStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
    NUMBER doBarzalaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
    bool checkConvergence(NUMBER* PI, NUMBER** A, NUMBER** B, NUMBER* PI_m1, NUMBER** A_m1, NUMBER** B_m1, bool flags[3]);
    // bridge to Liblinear - all objects are liblinear objects
    void createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space);
    void createLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space, NCAT k); // multiple KC's separately
    void recycleLiblinearProblem(struct problem &prob, struct parameter &param, struct feature_node *x_space);
private:
	bool checkPIABConstraints(NUMBER* a_PI, NUMBER** a_A, NUMBER** a_B); // all constraints, inc row sums
    
    // fitting methods - helpers (hidden)
    // fitting methods (hidden)
    NUMBER GradientDescentSkill(); // return -LL for the model
    NUMBER ConjugateGradientDescentSkill(); // return -LL for the model
    NUMBER GradientDescentGroup();
    NUMBER BaumWelchSkill();
    void doBaumWelchStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B);
//    void predictMetricsGroup(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill);
    // write model
	void toFileSkill(const char *filename);
	void toFileGroup(const char *filename);
};

#endif