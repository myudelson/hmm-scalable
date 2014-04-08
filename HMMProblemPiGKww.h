/*
 *  HMMProblemPiGKww.h
 *  HMM
 *  ww  - version with weighling of student and skill components (both global)
 *
 *  Created by Mikhail Yudelson on 11/3/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "HMMProblemPiGK.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

#ifndef __HMM__HMMProblemPiGKww__
#define __HMM__HMMProblemPiGKww__

class HMMProblemPiGKww : public HMMProblemPiGK {
public:
	HMMProblemPiGKww(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblemPiGKww();
//	virtual NUMBER** getPI();
//	NUMBER** getPIk();
//	NUMBER** getPIg();
//	virtual NUMBER* getPI(NCAT k);
//	NUMBER* getPIk(NCAT k);
//	NUMBER* getPIg(NCAT g);
//	virtual NUMBER** getA(NCAT k);
//	virtual NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
//    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
//    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(FitBit *fb);
    virtual void setGradWW(FitBit *fb);
	void toFile(const char *filename);
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER *ww;
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
//	void destroy(); // non-fit specific descruction
	virtual NDAT computeGradients(FitBit *fb);// NUMBER *a_gradPI, NUMBER** a_gradA, NUMBER **a_gradB);
    NUMBER GradientDescent(); // fit alternating
private:
};



#endif /* defined(__HMM__HMMProblemPiGKww__) */
