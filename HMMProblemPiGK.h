/*
 *  HMMProblemPiGK.h
 *  HMM
 *
 *  Created by Mikhail Yudelson on 5/3/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */
#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

#ifndef __HMM__HMMProblemPiGK__
#define __HMM__HMMProblemPiGK__

class HMMProblemPiGK : public HMMProblem {
public:
	HMMProblemPiGK(struct param *param); // sizes=={nK, nK, nK} by default
	NUMBER** getPI();
	NUMBER** getPIk();
	NUMBER** getPIg();
	NUMBER* getPI(NCAT k);
	NUMBER* getPIk(NCAT k);
	NUMBER* getPIg(NCAT g);
	NUMBER** getA(NCAT k);
	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
	void toFile(const char *filename);
    // fitting (the only public method)
    void fit(); // return -LL for the model
protected:
	//
	// Givens
	//
	NUMBER** PIg; // initial state probabilities
	//
	// Derived
	//
	NUMBER** gradPIk; // MULTIPLE gradients of initial state probabilities
	NUMBER** gradPIg; // MULTIPLE gradients of initial state probabilities
	NUMBER*** gradA; // MULTIPLE gradients of transition matrix
	NUMBER*** gradB; // MULTIPLE gradients of observation matrix
    bool* fitK; // flags for enabling the fittig of a skill
    bool* fitG; // flags for enabling the fittig of a group
    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	
	void init(struct param *param); // non-fit specific initialization
	void destroy(); // non-fit specific destroy
	void initGrad();
private:
    // fitting methods - helpers (hidden)
//	static void computeAlphaAndPOParam(NCAT xndat, struct data** x_data, NUMBER *a_PIk, NUMBER **a_PIg, NUMBER** a_A, NUMBER **a_B, NPAR a_nS);
//	void computeAlphaAndPOParamPIg(NCAT xndat, struct data** x_data, NUMBER *a_PI);
//	void computeBetaPIg(NCAT xndat, struct data** x_data);
	void computeGradients();
	void computeGradientsK(NCAT k, NUMBER* a_gradPI, NUMBER** a_gradA, NUMBER** a_gradB);
	void computeGradientsG(NCAT g, NUMBER* a_gradPI);
    // fitting methods (hidden)
//    NUMBER GradientDescentPLoSKillGroupOtherSkill0(); // fit in parallel
    NUMBER GradientDescentPLoSKillGroupOtherSkill1(); // fit alternating
    NUMBER doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
                           NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
    NUMBER doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER *a_gradPI);
};



#endif /* defined(__HMM__HMMProblemPiGK__) */
