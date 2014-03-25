//
//  HMMProblemPiABGK.h
//  HMM
//
//  Created by Yudelson, Michael on 5/10/13.
//
//

#ifndef __HMM__HMMProblemPiABGK__
#define __HMM__HMMProblemPiABGK__

#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

class HMMProblemPiABGK : public HMMProblem {
public:
	HMMProblemPiABGK(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblemPiABGK();
	NUMBER** getPI();
	NUMBER** getPIk();
	NUMBER** getPIg();
	NUMBER*** getA();
	NUMBER*** getAk();
	NUMBER*** getAg();
	NUMBER*** getBk();
	NUMBER*** getBg();
	NUMBER* getPI(NCAT k);
	NUMBER* getPIk(NCAT k);
	NUMBER* getPIg(NCAT g);
	NUMBER** getA(NCAT k);
	NUMBER** getAk(NCAT k);
	NUMBER** getAg(NCAT g);
	NUMBER** getB(NCAT k);
	NUMBER** getBk(NCAT k);
	NUMBER** getBg(NCAT g);
    // getters for computing alpha, beta, gamma
    NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
	void toFile(const char *filename);
    // fitting (the only public method)
    void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
    void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
    void setGradB (struct data* dt, FitBit *fb, NPAR kg_flag);
    void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER** PIg; // PI by group
	NUMBER*** Ag; // A by group
	NUMBER*** Bg; // B by group
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
    void destroy();
    //	void initGrad();
    NUMBER GradientDescent(); // fit alternating
private:
};


#endif /* defined(__HMM__HMMProblemPiABGK__) */
