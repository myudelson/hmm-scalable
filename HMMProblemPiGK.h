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
    ~HMMProblemPiGK();
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
    // getters for computing gradients of alpha, beta, gamma
    void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
	void toFile(const char *filename);
    // fitting (the only public method)
    void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER** PIg; // initial state probabilities
	//
	// Derived
	//
    bool* fitK; // flags for enabling the fittig of a skill
    bool* fitG; // flags for enabling the fittig of a group
    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	
	void init(struct param *param); // non-fit specific initialization
	void destroy(); // non-fit specific descruction
    NUMBER GradientDescent(); // fit alternating
private:
};



#endif /* defined(__HMM__HMMProblemPiGK__) */
