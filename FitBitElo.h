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

#ifndef __HMM__FitBitElo__
#define __HMM__FitBitElo__

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
//
//enum FIT_BIT_VAR {
//    FBV_PI = 1, // PI
//    FBV_A  = 2, // A
//    FBV_B  = 3  // B
//};
#endif /* fit bit enums*/

class FitBitElo {
public:
//    NPAR nO, nS; // copies
    NCAT nG;//, nK; // copies
    NPAR   nELO; // copy
    NUMBER *ELO; // usually pointer
    NUMBER *ELOm1; // previous value
    NUMBER *ELOm2; // previous previous value
    NUMBER *gradELO; // gradient
    NUMBER *gradELOm1; // previous gradient
    NUMBER *gradELOcopy; // gradient copy
    NUMBER *ELOcopy; // copy
    NUMBER *dirELO; // step direction
    NUMBER *dirELOm1; // previous step direction
    struct param *param; // problem parameters
    NPAR tag; // multippurpose
    
    NUMBER *elo_track_g; // elo g (student) rating
    NCAT *elo_count_g;   // per g (student) count
    NUMBER *elo_track_g_t; // elo g (student) rating for row, not student
    NUMBER *elo_grad_error_g; // per g (student) error
    NUMBER ll; // loglik

    FitBitElo(/*NPAR a_nS, NPAR a_nO, NCAT a_nK,*/ NCAT a_nG, NPAR a_nELO, NUMBER a_tol, NPAR a_tol_mode);
    ~FitBitElo();
    void init(enum FIT_BIT_SLOT fbs);
    void negate(enum FIT_BIT_SLOT fbs);
    void link(NUMBER *a_ELO, NUMBER *a_elo_track_g, NCAT *a_elo_count_g, NUMBER *a_elo_track_g_t, struct param *a_param);
    void toZero(enum FIT_BIT_SLOT fbs); // since we reset to 0 logits == 0.5 probability, rather use to05 to avoid confusion
    void destroy(enum FIT_BIT_SLOT fbs);
    void copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    void add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    bool checkConvergence(FitResult *fr);
    bool checkConvergenceSingle(FitResult *fr); // without checking for oscillation, actually, afer copying t-1 to t-2 and t to t-1, it is used to check for oscillation
    void doLog10ScaleGentle(enum FIT_BIT_SLOT fbs);
//	void doLog10ScaleGentleByRow(enum FIT_BIT_SLOT fbs);

    // adding penalties
//    void addL2Penalty(enum FIT_BIT_VAR fbv, param* param, NUMBER factor);
private:
    NUMBER tol;
    NPAR tol_mode;

    void init(NUMBER* &a_ELO);
    void negate(NUMBER* &a_ELO);
    void toZero(NUMBER *a_ELO);
    void destroy(NUMBER* &a_ELO);
    void get(enum FIT_BIT_SLOT fbs, NUMBER* &a_ELO);
    void add(NUMBER *sourseELO, NUMBER *targetELO);
    void copy(NUMBER* &sourseELO, NUMBER* &targetELO);
};

#endif /* defined(__HMM__FitBitElo__) */
