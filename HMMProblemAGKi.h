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
};

#endif /* defined(__HMM__HMMProblemAGKi__) */
