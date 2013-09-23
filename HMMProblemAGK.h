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
    ~HMMProblemAGK();
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
    NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
	void toFile(const char *filename);
    // fitting (the only public method)
    void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER*** Ag; // A by group
	//
	// Derived
	//
    void init(struct param *param); // non-fit specific initialization
	void destroy(); // non-fit specific descruction
//	void initGrad();
    NUMBER GradientDescent(); // fit alternating
private:
};

#endif /* defined(__HMM__HMMProblemAGK__) */
