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

//
// encapsulator of results of a fitting [sub-]job
//

#ifndef __HMM__FitBitSt__
#define __HMM__FitBitSt__

#include "utilsSt.h"

// return type cummrizing a result of a fir of the [subset of the] data
struct FitResult {
    int iter;   // interations
    NUMBER pO0_prefit; // starting log-likelihood before all iterations
    NUMBER pO0; // starting log-likelihood on current iteration
    NUMBER pO;  // final log-likelihood
    NPAR conv;   // converged? (maybe just went to max-iter)
    NDAT ndat;
};

class FitBitSt {
public:
    NPAR nO, nS; // copies
    NCAT nG, nK; // copies
    
    NCAT   nPARAM; // copy
    NUMBER *PARAM; // usually pointer
    NUMBER *PARAMm1; // previous value
    NUMBER *PARAMm2; // previous previous value
    NUMBER *PARAMcopy; // copy
    NUMBER *GRAD; // gradient
    NUMBER *GRADm1; // previous gradient
    NUMBER *GRADcopy; // gradient copy
    NUMBER *DIR; // step direction
    NUMBER *DIRm1; // previous step direction
    
    FitResult *fit_results;
    FitResult fit_result; // overall

    struct task* task; // sequences of data
    NPAR projecttosimplex; // whether projection to simplex should be done
    NPAR Cslice; // current slice during L2 norm penalty fitting
    NPAR tag; // multippurpose
    
    NCAT* sclsmplx_offsets;// start positions of the scale / simplex vectors in the gradient and parameter arrays
    NCAT* sclsmplx_bound_offsets;// start positions of the simplex boundary vectors in the gradient and parameter arrays
    NPAR* sclsmplx_sizes;// lengths of the scale / simplex vectors in the gradient and parameter arrays
    NCAT* sclsmplx_blocks;// id of a block that the simplex belongs to
    NCAT* param_blocks;// id of a block that the parameter belongs to (parameter to skill)
    NCAT sclsmplx_n; // number of scale / simplex vectors
    
    FitBitSt(struct task *task, NUMBER* a_PARAM, NUMBER* a_GRAD, NCAT a_nPARAM);
    ~FitBitSt();
//    void init(enum FIT_BIT_SLOT fbs);
//    void negate(enum FIT_BIT_SLOT fbs);
//    void link(NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NCAT a_xndat, struct data** a_x_data);
//    void toZero(enum FIT_BIT_SLOT fbs);
//    void destroy(enum FIT_BIT_SLOT fbs);
//    void copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
//    void add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    NCAT checkConvergence(NCAT nbits, NDAT *offsets, NDAT* sizes, NCAT* blocks, NPAR *condition);
    NCAT checkConvergenceSingle(NCAT nbits, NDAT *offsets, NDAT* sizes, NCAT* blocks, NPAR *condition); // without checking for oscillation, actually, afer copying t-1 to t-2     bool checkConvergenceBit(NDAT offset, NDAT size, FitResult *fr);
    NPAR checkConvergenceBit(NDAT offset, NDAT size, FitResult *fr);
    NPAR checkConvergenceBitSingle(NDAT offset, NDAT size, FitResult *fr); // without checking for oscillation, actually, afer copying t-1 to t-2 and t to t-1, it is used to check for oscillation
    void projectToSimplex(bool hasNon01Constraints, NPAR* condition);
    void scaleGradients(bool direction, NPAR* condition);
    void printFitResult(NCAT k);
//    void doLog10ScaleGentle(NDAT offset, NDAT size); // REDUNTANT! utils.h, utilsSt.h have log10-scalers, no need to copy shit
//	void doLog10ScaleGentleByRow(enum FIT_BIT_SLOT fbs);  // REDUNTANT!, utilsSt.h have log10-scalers, no need to copy shit

    // adding penalties
    // void addL2Penalty(enum FIT_BIT_VAR fbv, struct task* task, NUMBER factor);
private:
    NUMBER tol;
    NPAR tol_mode;
};

#endif /* defined(__HMM__FitBitSt__) */
