/*
 
 Copyright (c) 2012-2015, Michael (Mikhail) Yudelson
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

#ifndef __HMM__FitBit__
#define __HMM__FitBit__

#include "utils.h"

#ifndef __HMM__FitBitEnums__
#define __HMM__FitBitEnums__
enum FIT_BIT_SLOT {
    FBS_PAR       = 1, // e.g. PI
    FBS_PARm1     = 2, // e.g. PIm1
    FBS_PARm2     = 7, // e.g. PIm2
    FBS_GRAD      = 3, // e.g. gradPI
    FBS_GRADm1    = 4, // e.g. gradPIm1
    FBS_GRADcopy  = 9, // e.g. gradPIcopy
    FBS_PARcopy   = 5, // e.g. PIcopy
    FBS_DIR       = 8, // e.g. gradPIdir
    FBS_DIRm1     = 6  // e.g. gradPIdirm1
};

enum FIT_BIT_VAR {
    FBV_PI = 1, // PI
    FBV_A  = 2, // A
    FBV_B  = 3  // B
};
#endif /* fit bit enums*/

class FitBit {
public:
    NPAR nO, nS; // copies
    NCAT nG, nK; // copies
    NUMBER* pi; // usually pointer
    NUMBER **A; // usually pointer
    NUMBER **B; // usually pointer
    NUMBER *PIm1; // previous value
    NUMBER **Am1; // previous value
    NUMBER **Bm1; // previous value
    NUMBER *PIm2; // previous previous value
    NUMBER **Am2; // previous previous value
    NUMBER **Bm2; // previous previous value
    NUMBER *gradPI; // gradient
    NUMBER **gradA; // gradient
    NUMBER **gradB; // gradient
    NUMBER *gradPIm1; // previous gradient
    NUMBER **gradAm1; // previous gradient
    NUMBER **gradBm1; // previous gradient
    NUMBER *gradPIcopy; // gradient copy
    NUMBER **gradAcopy; // gradient copy
    NUMBER **gradBcopy; // gradient copy
    NUMBER *PIcopy; // copy
    NUMBER **Acopy; // copy
    NUMBER **Bcopy; // copy
    NUMBER *dirPI; // step direction
    NUMBER **dirA; // step direction
    NUMBER **dirB; // step direction
    NUMBER *dirPIm1; // previous step direction
    NUMBER **dirAm1; // previous step direction
    NUMBER **dirBm1; // previous step direction
    NCAT xndat; // number of sequences of data
    struct data** x_data; // sequences of data
    NPAR projecttosimplex; // whether projection to simplex should be done
    NPAR Cslice; // current slice during L2 norm penalty fitting
    NPAR tag; // multippurpose
    
    FitBit(NPAR a_nS, NPAR a_nO, NCAT a_nK, NCAT a_nG, NUMBER a_tol, NPAR a_tol_mode);
    FitBit(NPAR a_nS, NPAR a_nO, NCAT a_nK, NCAT a_nG, NUMBER a_tol, NPAR a_tol_mode, NPAR a_projecttosimplex);
    ~FitBit();
    void init(enum FIT_BIT_SLOT fbs);
    void negate(enum FIT_BIT_SLOT fbs);
    void link(NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NCAT a_xndat, struct data** a_x_data);
    void toZero(enum FIT_BIT_SLOT fbs);
    void destroy(enum FIT_BIT_SLOT fbs);
    void copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    void add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    bool checkConvergence(FitResult *fr);
    bool checkConvergenceSingle(FitResult *fr); // without checking for oscillation, actually, afer copying t-1 to t-2 and t to t-1, it is used to check for oscillation
    void doLog10ScaleGentle(enum FIT_BIT_SLOT fbs);
	void doLog10ScaleGentleByRow(enum FIT_BIT_SLOT fbs);

    // adding penalties
    void addL2Penalty(enum FIT_BIT_VAR fbv, param* param, NUMBER factor);
private:
    NUMBER tol;
    NPAR tol_mode;

    void init(NUMBER* &pi, NUMBER** &A, NUMBER** &B);
    void negate(NUMBER* &pi, NUMBER** &A, NUMBER** &B);
    void toZero(NUMBER *pi, NUMBER **A, NUMBER **B);
    void destroy(NUMBER* &pi, NUMBER** &A, NUMBER** &B);
    void get(enum FIT_BIT_SLOT fbs, NUMBER* &a_PI, NUMBER** &a_A, NUMBER** &a_B);
    void add(NUMBER *soursePI, NUMBER **sourseA, NUMBER **sourseB, NUMBER *targetPI, NUMBER **targetA, NUMBER **targetB);
    void copy(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB);
};

#endif /* defined(__HMM__FitBit__) */
