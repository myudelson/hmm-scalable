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

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>


//#include "liblinear/linear.h"

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
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j);
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m);
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
    virtual void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
    virtual void setGradB (struct data* dt, FitBit *fb, NPAR kg_flag);
	virtual void toFile(const char *filename);
	NUMBER getSumLogPOPara(NCAT xndat, struct data **x_data); // generic per k/g-slice
	bool hasNon01Constraints();
    NUMBER getLogLik(); // get log likelihood of the fitted model
    NCAT getNparams(); // get log likelihood of the fitted model
    NUMBER getNullSkillObs(NPAR m); // get log likelihood of the fitted model
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // predicting
    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
    virtual void producePCorrect(boost::numeric::ublas::mapped_matrix<NUMBER*> *group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
    void predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill, bool only_unlabeled);
    void readModel(const char *filename, bool overwrite);
    virtual void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
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
	virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
	void initAlpha(NCAT xndat, struct data** x_data); // generic
	void initXi(NCAT xndat, struct data** x_data); // generic
	void initGamma(NCAT xndat, struct data** x_data); // generic
	void initBeta(NCAT xndat, struct data** x_data); // generic
	void computeAlphaAndPOParam(NCAT xndat, struct data** x_data);
	void computeBeta(NCAT xndat, struct data** x_data);
	void computeGamma(NCAT xndat, struct data** x_data);
	void computeXi(NCAT xndat, struct data** x_data);
    void FitNullSkill(NUMBER* loglik_rmse, bool keep_SE); // get loglik and RMSE
    // helpers
    void init3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
    void toZero3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS, NPAR nO);
    void free3Params(NUMBER* &PI, NUMBER** &A, NUMBER** &B, NPAR nS);
    void cpy3Params(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB, NPAR nS, NPAR nO);
    // predicting
	virtual void computeGradients(NCAT xndat, struct data** x_data, FitBit *fb, NPAR kg_flag);// NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB);
    virtual NUMBER doLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy);
    NUMBER doConjugateLinearStep(NCAT xndat, struct data** x_data, FitBit *fb, NCAT copy);
    NUMBER doBarzalaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
    FitResult GradientDescentBit(NCAT x, NCAT xndat, struct data** x_data, NPAR kg_flag, FitBit *fb, bool is1SkillForAll); // for 1 skill or 1 group, all 1 skill for all data
    virtual NUMBER GradientDescent(); // return -LL for the model
    void readNullObsRatio(FILE *fid, struct param* param, NDAT *line_no);
	bool checkPIABConstraints(NUMBER* a_PI, NUMBER** a_A, NUMBER** a_B); // all constraints, inc row sums
private:
    // fitting methods (hidden)
    NUMBER BaumWelchSkill();
    void doBaumWelchStep(NCAT xndat, struct data** x_data, FitBit *fb);//, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B);
    // write model
	void toFileSkill(const char *filename);
	void toFileGroup(const char *filename);
};

#endif