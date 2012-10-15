//
//  HMMProblemKT.h
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
//

#ifndef __HMM__HMMProblemKT__
#define __HMM__HMMProblemKT__

#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

#include "liblinear/linear.h"

class HMMProblemKT : public HMMProblem {
public:
	HMMProblemKT(struct param *param); // sizes=={nK, nK, nK} by default
	virtual ~HMMProblemKT();
//	NUMBER** getPI();
//	NUMBER*** getA();
//	NUMBER*** getB();
//	NUMBER* getPI(NCAT k);
//	NUMBER** getA(NCAT k);
//	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j);
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m);
	virtual void toFile(const char *filename);
    //	static NUMBER getSumNegPOPara(NCAT xndat, struct data **x_data, NUMBER (*f)(NUMBER)); // generic per k/g-slice
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // compute metrics
    void computeMetrics(NUMBER* metrics);
    // predicting
    virtual void producePCorrect(NUMBER*** group_skill_map, NUMBER* local_pred, NCAT* ks, NCAT nks, struct data* dt);
protected:
	//
	// Givens
	//
    NUMBER *liblinear_weights;
    
	//
	// Derived
	//
	virtual void init(struct param *param); // non-fit specific initialization
//	virtual void destroy(); // non-fit specific descruction
private:
    // fitting methods (hidden)
    NUMBER GradientDescentKT1(); // return -LL for the model
    NUMBER GradientDescentKT2(); // return -LL for the model
};

#endif /* defined(__HMM__HMMProblemKT__) */
