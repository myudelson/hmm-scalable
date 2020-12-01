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

#ifndef __HMM__HMMProblemEloK__
#define __HMM__HMMProblemEloK__

#include "HMMProblem.h"
#include "FitBitElo.h"
#include "StripedArray.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <string>

class HMMProblemEloK : public HMMProblem {
public:
	HMMProblemEloK(struct param *param); // sizes=={nK, nK, nK} by default
    ~HMMProblemEloK();
	NUMBER** getPI();
	NUMBER*** getA();
	NUMBER*** getB();
	NUMBER* getPI(NCAT k);
	NUMBER** getA(NCAT k);
	NUMBER** getB(NCAT k);
    // getters for computing alpha, beta, gamma
    NUMBER getPI(struct data* dt, NPAR i);         // to be redefined
    NUMBER getA (struct data* dt, NPAR i, NPAR j); // same
    NUMBER getB (struct data* dt, NPAR i, NPAR m); // same
	void toFile(const char *filename);
    // fitting (the only public method)
    void setGradPI(FitBit *fb);
    void setGradA (FitBit *fb);
    void setGradB (FitBit *fb);
    void fit(); // return -LL for the model
    void readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite);
    FitResult GradientDescentBitElo(FitBitElo *fb);
    virtual void postprocesslocal(NUMBER target, NUMBER estimate, NDAT t); // reserved for any per-row post-processing after producePCorrect, Elo uses it to update its tracked values
    virtual void preprocessglobal(); // reserved for any pre-processing before prediction, Elo uses it to update its tracked values
//    virtual void preprocess_computeAlphaAndPOParam(); // nothing for BKT, computes Elo values for Elo-based
protected:
	//
	// Givens
	//
	NUMBER K; // Elo sensitivity K
    NUMBER *elo_track_g; // elo g (student) rating
    NCAT *elo_count_g;   // per g (student) count
    NUMBER *elo_track_g_t; // elo g (student) rating for row, not student
    
	//
	// Derived
	//
	void init(struct param *param); // non-fit specific initialization
    void destroy();
    //	void initGrad();
    virtual NUMBER doLinearStepElo(FitBitElo *fbe);
    virtual NUMBER doConjugateLinearStepElo(FitBitElo *fbe);
    virtual NDAT computeGradientsElo(FitBitElo *fbe);
    NUMBER GradientDescent(); // fit alternating
private:
};


#endif /* defined(__HMM__HMMProblemEloK__) */
