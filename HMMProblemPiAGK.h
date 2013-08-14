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
    virtual ~HMMProblemPiAGK();
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
    void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
    void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
    void fit(); // return -LL for the model
protected:
	//
	// Givens
	//
	NUMBER** PIg; // PI by group
	NUMBER*** Ag; // A by group
	//
	// Derived
	//
    bool* fitK; // flags for enabling the fittig of a skill
    bool* fitG; // flags for enabling the fittig of a group
    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	virtual void init(struct param *param); // non-fit specific initialization
    virtual void destroy();
    NUMBER GradientDescent(); // fit alternating
    virtual void readModelBody(FILE *fid, struct param* param, NDAT *line_no);
private:
};


#endif /* defined(__HMM__HMMProblemPiAGK__) */
