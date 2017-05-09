/*
 
 Copyright (c) 2012-2017, Michael (Mikhail) Yudelson
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

#include "utils.h"
#include "HMMProblem.h"
#include "FitBitSlicedA.h"
#include "StripedArray.h"

//#include <boost/numeric/ublas/matrix_sparse.hpp>//BOOST
//#include <boost/numeric/ublas/io.hpp>//BOOST

#ifndef _HMMPROBLEM_SLICED_A_H
#define _HMMPROBLEM_SLICED_A_H

class HMMProblemSlicedA : public HMMProblem {
public:
	HMMProblemSlicedA();
	HMMProblemSlicedA(struct param *param); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblemSlicedA();
    NUMBER** getPI();
	NUMBER**** getA();
	NUMBER*** getB();
	NUMBER* getPI(NCAT k);
    NUMBER*** getA(NCAT k);
    NUMBER** getA(NCAT k, NPAR z);
	NUMBER** getB(NCAT k);
	NUMBER*** getLbA();
	NUMBER*** getUbA();
    // getters for computing alpha, beta, gamma
    virtual NUMBER getPI(struct data* dt, NPAR i);
    virtual NUMBER getA (struct data* dt, NPAR i, NPAR j);
    virtual NUMBER getB (struct data* dt, NPAR i, NPAR m);
    // getters for computing gradients of alpha, beta, gamma
    virtual void setGradPI(FitBitSlicedA *fb);
    virtual void setGradA (FitBitSlicedA *fb);
    virtual void setGradB (FitBitSlicedA *fb);
	virtual void toFile(const char *filename);
    // fitting (the only public method)
    virtual void fit(); // return -LL for the model
    // predicting
    void readModel(const char *filename, bool overwrite);
    virtual void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
protected:
	//
	// Givens
	//
	NUMBER**** A; // transition matrix
	NUMBER*** lbA; // lower boundary transition matrix
	NUMBER*** ubA; // upper boundary transition matrix
	//
	// Derived
	//
	virtual void init(struct param *param); // non-fit specific initialization
	virtual void destroy(); // non-fit specific descruction
	NDAT computeAlphaAndPOParam(NCAT xndat, struct data** x_data);
	void computeBeta(NCAT xndat, struct data** x_data);
	void computeXiGamma(NCAT xndat, struct data** x_data);
    // helpers
    void init3Params(NUMBER* &pi, NUMBER*** &A, NUMBER** &B, NPAR nZ, NPAR nS, NPAR nO); // sliced
    void toZero3Params(NUMBER* &pi, NUMBER*** &A, NUMBER** &B, NPAR nZ, NPAR nS, NPAR nO);
    void free3Params(NUMBER* &pi, NUMBER*** &A, NUMBER** &B, NPAR nZ, NPAR nS); // sliced
    void cpy3Params(NUMBER* &soursePI, NUMBER*** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER*** &targetA, NUMBER** &targetB, NPAR nZ, NPAR nS, NPAR nO);// sliced

    // predicting
	virtual NDAT computeGradients(FitBitSlicedA *fb);
    virtual NUMBER doLinearStep(FitBitSlicedA *fb);
    virtual NUMBER doLagrangeStep(FitBitSlicedA *fb);
    NUMBER doConjugateLinearStep(FitBitSlicedA *fb);
    NUMBER doBaumWelchStep(FitBitSlicedA *fb);
    FitResult GradientDescentBit(FitBitSlicedA *fb);  // for 1 skill or 1 group, all 1 skill for all data
    FitResult BaumWelchBit(FitBitSlicedA *fb);
    
    NUMBER doBarzilaiBorweinStep(FitBitSlicedA *fb);
    virtual NUMBER GradientDescent(); // return -LL for the model
    NUMBER BaumWelch(); // return -LL for the model
	bool checkPIABConstraints(NUMBER* a_PI, NUMBER*** a_A, NUMBER** a_B); // all constraints, inc row sums
private:
    // write model
	void toFileSkill(const char *filename);
	void toFileGroup(const char *filename);
};

#endif
