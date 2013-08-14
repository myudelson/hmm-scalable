//
//  HMMProblemAGK.h
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
//

#ifndef __HMM__HMMProblemAGK__
#define __HMM__HMMProblemAGK__

#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

class HMMProblemAGK : public HMMProblem {
public:
	HMMProblemAGK(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblemAGK();
	NUMBER** getPI();
	NUMBER*** getA();
	NUMBER*** getAk();
	NUMBER*** getAg();
	NUMBER* getPI(NCAT k);
	NUMBER** getA(NCAT k);
	NUMBER** getAk(NCAT k);
	NUMBER** getAg(NCAT g);
	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    virtual void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
	void toFile(const char *filename);
    // fitting (the only public method)
    void fit(); // return -LL for the model
protected:
	//
	// Givens
	//
	NUMBER*** Ag; // A by group
	//
	// Derived
	//
    virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
//	void initGrad();
    virtual NUMBER GradientDescent(); // fit alternating
    virtual void readModelBody(FILE *fid, struct param* param, NDAT *line_no);
private:
};

#endif /* defined(__HMM__HMMProblemAGK__) */
