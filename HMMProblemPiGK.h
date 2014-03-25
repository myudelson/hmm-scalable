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
//	NUMBER** getPI();
	virtual NUMBER** getPIk();
	virtual NUMBER** getPIg();
//	virtual NUMBER* getPI(NCAT k);
	virtual NUMBER* getPIk(NCAT k);
	virtual NUMBER* getPIg(NCAT g);
//	virtual NUMBER** getA(NCAT k);
//	virtual NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
	void toFile(const char *filename);
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER** PIg; // initial state probabilities
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
	void destroy(); // non-fit specific descruction
    NUMBER GradientDescent(); // fit alternating
private:
};



#endif /* defined(__HMM__HMMProblemPiGK__) */
