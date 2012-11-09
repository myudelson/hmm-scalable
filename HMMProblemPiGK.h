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
    virtual ~HMMProblemPiGK();
	NUMBER** getPI();
	NUMBER** getPIk();
	NUMBER** getPIg();
	NUMBER* getPI(NCAT k);
	NUMBER* getPIk(NCAT k);
	NUMBER* getPIg(NCAT g);
	NUMBER** getA(NCAT k);
	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
//    virtual void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
//    virtual void setGradB (struct data* dt, FitBit *fb, NPAR kg_flag);
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
//	NUMBER** gradPI; // MULTIPLE gradients of initial state probabilities
//	NUMBER** gradPIg; // MULTIPLE gradients of initial state probabilities
//	NUMBER*** gradA; // MULTIPLE gradients of transition matrix
//	NUMBER*** gradB; // MULTIPLE gradients of observation matrix
    bool* fitK; // flags for enabling the fittig of a skill
    bool* fitG; // flags for enabling the fittig of a group
    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	
	virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
	void initGrad();
private:
    // fitting methods - helpers (hidden)
//	virtual void computeGradients (NCAT xndat, struct data** x_data, FitBit *fb);
//	virtual void computeGradientsG(NCAT xndat, struct data** x_data, FitBit *fb);
    // fitting methods (hidden)
    NUMBER GradientDescentPLoSKillGroupOtherSkill(NPAR kg_flag); // fit alternating
//    NUMBER doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
//                           NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
//    NUMBER doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER *a_gradPI);
};



#endif /* defined(__HMM__HMMProblemPiGK__) */
