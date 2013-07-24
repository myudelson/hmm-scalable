//
//  HMMProblemAGKi.h
//  HMM
//
//  Created by Mikhail Yudelson on 9/13/12.
//
//

#ifndef __HMM__HMMProblemAGKi__
#define __HMM__HMMProblemAGKi__

#include "HMMProblemAGK.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

class HMMProblemAGKi : public HMMProblemAGK {
public:
	HMMProblemAGKi(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblemAGKi();
    virtual void setGradA (struct data* dt, FitBit *fb, NPAR kg_flag);
    void fit(); // return -LL for the model
protected:
    virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
    virtual NUMBER GradientDescent(); // fit alternating
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

#endif /* defined(__HMM__HMMProblemAGKi__) */
