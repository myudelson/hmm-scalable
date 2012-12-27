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
//	NUMBER** gradPI; // MULTIPLE gradients of initial state probabilities
//	NUMBER** gradPIg; // MULTIPLE gradients of initial state probabilities
//	NUMBER*** gradA; // MULTIPLE gradients of transition matrix
//	NUMBER*** gradAg; // MULTIPLE gradients of transition matrix
//	NUMBER*** gradB; // MULTIPLE gradients of observation matrix
    bool* fitK; // flags for enabling the fittig of a skill
    bool* fitG; // flags for enabling the fittig of a group
    NCAT* fitK_countG; // number of groups per skill that have a raised fitG flag
	
	virtual void init(struct param *param); // non-fit specific initialization
    virtual void destroy();
//	void initGrad();
    NUMBER GradientDescent(); // fit alternating
    virtual void readModel(FILE *fid, NDAT *line_no);
private:
//	void computeGradients();
//	void computeGradientsK(NCAT k, NUMBER* a_gradPI, NUMBER** a_gradA, NUMBER** a_gradB);
//	void computeGradientsG(NCAT g, NUMBER* a_gradPI, NUMBER** a_gradA);
    // fitting methods (hidden)
//    NUMBER GradientDescentPLoSKillGroupOtherSkill(); // fit alternating
//    NUMBER doLinearStepSkill(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B,
//                           NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB);
//    NUMBER doLinearStepPLoGroup(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER *a_gradPI, NUMBER **a_gradA);
    // predicting
//	static void computeLogLikRMSE(NUMBER* loglik_rmse, bool keep_SE, NCAT xndat, struct data** x_data, NUMBER *a_PIk, NUMBER **a_PIg, NUMBER** a_Ak, NUMBER*** a_Ag, NUMBER **a_B, struct param* param);
};


#endif /* defined(__HMM__HMMProblemPiAGK__) */
