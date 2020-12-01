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
//    void predict(NUMBER* metrics, const char *filename, StripedArray<NPAR> *dat_obs, StripedArray<NCAT> *dat_group, StripedArray<NCAT> *dat_skill, StripedArray<NCAT*> *dat_multiskill);
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
