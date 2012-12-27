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
    // getters for computing gradients of alpha, beta, gamma
//    virtual void setGradPI(struct data* dt, FitBit *fb, NPAR kg_flag);
    virtual void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
//    virtual void setGradB (struct data* dt, FitBit *fb, NPAR kg_flag);
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
//	NUMBER** gradPI; // MULTIPLE gradients of initial state probabilities
//	NUMBER*** gradA; // MULTIPLE gradients of transition matrix
//	NUMBER*** gradAg; // MULTIPLE gradients of transition matrix
//	NUMBER*** gradB; // MULTIPLE gradients of observation matrix
//    bool* fitK; // flags for enabling the fittig of a skill
//    bool* fitG; // flags for enabling the fittig of a group
//    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	
    virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
//	void initGrad();
    virtual NUMBER GradientDescent(); // fit alternating
    virtual void readModel(FILE *fid, NDAT *line_no);
private:
//	void computeGradients();
//	virtual void computeGradients (NCAT xndat, struct data** x_data, FitBit *fb);
//	virtual void computeGradientsG(NCAT xndat, struct data** x_data, FitBit *fb);
    // fitting methods (hidden)
//    NUMBER doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
//                             NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
//    NUMBER doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER **a_A, NUMBER **a_gradA);
    // predicting
//	static void computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data, NUMBER*a_PIk, NUMBER** a_Ak, NUMBER*** a_Ag, NUMBER **a_B, struct param* param);
};

#endif /* defined(__HMM__HMMProblemAGK__) */
