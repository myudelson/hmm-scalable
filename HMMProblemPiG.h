//
//  HMMProblemPiG.h
//  HMM
//
//  Created by Mikhail Yudelson on 8/24/12.
//
//

#ifndef __HMM__HMMProblemPiG__
#define __HMM__HMMProblemPiG__

#include "StripedArray.h"
#include "HMMProblem.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

//class HMMProblemPiG;
//
//typedef NUMBER (HMMProblemPiG::*getNum1D)(struct data* dt, NPAR i);
//typedef NUMBER (HMMProblemPiG::*getNum2D)(struct data* dt, NPAR i, NPAR j);

class HMMProblemPiG : public HMMProblem {
public:
//	NUMBER** getPI();
//	NUMBER*** getA();
//	NUMBER*** getB();
//	NUMBER* getPI(NCAT k);
//	NUMBER** getA(NCAT k);
//	NUMBER** getB(NCAT k);
	HMMProblemPiG(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblemPiG();
    // writing
	void toFile(const char *filename);
    // predicting
    void predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill);
    // fitting
    void fit(); // return -LL for the model
    //
    NUMBER getPI(struct data* dt, NPAR i);
    NUMBER getA (struct data* dt, NPAR i, NPAR j);
    NUMBER getB (struct data* dt, NPAR i, NPAR m);
private:
	//
	// Derived
	//
    // fitting - supportive
//	static void initGradPI(NUMBER* &a_gradPI, struct param* param); // generic
//	static void initGradAB( NUMBER** &a_gradA, NUMBER** &a_gradB, struct param* param); // generic
//	static void computeAlphaAndPOParamPI(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER*** a_A, NUMBER ***a_B, NPAR a_nS);
//	static void computeAlphaAndPOParamAB(NCAT xndat, struct data** x_data, NUMBER **a_PI, NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
//	static void computeBetaPI(NCAT xndat, struct data** x_data, NUMBER*** a_A, NUMBER ***a_B, NPAR a_nS);
//	static void computeBetaAB(NCAT xndat, struct data** x_data, NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
//	static void computeXiPI(NCAT xndat, struct data** x_data,	NUMBER*** a_A, NUMBER ***a_B, NPAR a_nS);
//	static void computeXiAB(NCAT xndat, struct data** x_data,	NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
	void computeGradientsPI(NCAT xndat, struct data** x_data, NUMBER *a_gradPI);
	void computeGradientsAB(NCAT xndat, struct data** x_data, NUMBER** a_gradA, NUMBER **a_gradB);
    // fitting
    NUMBER GradientDescentPLoGroupOtherSkill();
    NUMBER doLinearStepPLoGroupOtherPI(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER *a_gradPI);
    NUMBER doLinearStepPLoGroupOtherAB(NCAT xndat, struct data** x_data, NUMBER **a_A, NUMBER **a_B,
                                     NUMBER **a_gradA, NUMBER **a_gradB);
    // predicting
//	static void computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER*** a_A, NUMBER ***a_B, struct param* param);
};


#endif /* defined(__HMM__HMMProblemPiG__) */
