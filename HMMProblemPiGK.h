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

#include "HMMProblem.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

#ifndef __HMM__HMMProblemPiGK__
#define __HMM__HMMProblemPiGK__

class HMMProblemPiGK : public HMMProblem {
public:
	HMMProblemPiGK(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblemPiGK();
//	NUMBER** getPI();
	virtual NUMBER** getPIk();
	virtual NUMBER** getPIg();
//	virtual NUMBER* getPI(NCAT k);
	virtual NUMBER* getPIk(NCAT k);
	virtual NUMBER* getPIg(NCAT g);
//	virtual NUMBER** getA(NCAT k);
//	virtual NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(FitBit *fb);
	void toFile(const char *filename);
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER** PIg; // initial state probabilities
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
	void destroy(); // non-fit specific descruction
    virtual NUMBER GradientDescent(); // fit alternating
    virtual struct data*** getExdendedData(NCAT xndat, struct data** x_data, NPAR kg_flag, NCAT* xxndat);
//    virtual NUMBER GradientDescentX(); // fit alternating

private:
};



#endif /* defined(__HMM__HMMProblemPiGK__) */
