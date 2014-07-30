/*
 
 Copyright (c) 2012-2014, Michael (Mikhail) Yudelson
 All rights reserved.
 
 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of the Michael (Mikhail) Yudelson nor the
 names of other contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL COPYRIGHT HOLDERS AND CONTRIBUTORS BE LIABLE FOR
 ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 */

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
