//
//  HMMProblemPiAGK.h
//  HMM
//
//  Created by Mikhail Yudelson on 9/6/12.
//
//

#ifndef __HMM__HMMProblemPiAGK__
#define __HMM__HMMProblemPiAGK__

#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

class HMMProblemPiAGK : public HMMProblem {
public:
	HMMProblemPiAGK(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblemPiAGK();
	NUMBER** getPI();
	NUMBER** getPIk();
	NUMBER** getPIg();
	NUMBER*** getA();
	NUMBER*** getAk();
	NUMBER*** getAg();
	NUMBER* getPI(NCAT k);
	NUMBER* getPIk(NCAT k);
	NUMBER* getPIg(NCAT g);
	NUMBER** getA(NCAT k);
	NUMBER** getAk(NCAT k);
	NUMBER** getAg(NCAT g);
	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
	void toFile(const char *filename);
    // fitting (the only public method)
    void setGradPI(FitBit *fb);
    void setGradA (FitBit *fb);
    void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER** PIg; // PI by group
	NUMBER*** Ag; // A by group
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
    void destroy();
    NUMBER GradientDescent(); // fit alternating
private:
};


#endif /* defined(__HMM__HMMProblemPiAGK__) */
