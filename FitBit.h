//
//  FitBit.h
//  HMM
//
//  Created by Mikhail Yudelson on 10/12/12.
//
//

#ifndef __HMM__FitBit__
#define __HMM__FitBit__

#include <iostream>
#include "utils.h"

enum FIT_BIT_SLOT {
    FBS_PAR       = 1, // e.g. PI
    FBS_PARm1     = 2, // e.g. PIm1
    FBS_GRAD      = 3, // e.g. gradPI
    FBS_GRADm1    = 4, // e.g. gradPIm1
    FBS_PARcopy   = 5, // e.g. PIcopy
    FBS_DIRm1     = 6  // e.g. gradPIdirm1
};

enum FIT_BIT_VAR {
    FBV_PI = 1, // PI
    FBV_A  = 2, // A
    FBV_B  = 3  // B
};

class FitBit {
public:
    NPAR nO, nS; // copies
    NCAT nG, nK; // copies
    NUMBER *PI; // usually pointer #1
    NUMBER **A; // usually pointer
    NUMBER **B; // usually pointer
    NUMBER *PIm1; // previous value #2
    NUMBER **Am1; // previous value
    NUMBER **Bm1; // previous value
    NUMBER *gradPI; // gradient #3
    NUMBER **gradA; // gradient
    NUMBER **gradB; // gradient
    NUMBER *gradPIm1; // previous gradient #4
    NUMBER **gradAm1; // previous gradient
    NUMBER **gradBm1; // previous gradient
    NUMBER *PIcopy; // previous value #5
    NUMBER **Acopy; // previous value
    NUMBER **Bcopy; // previous value
    NUMBER *dirPIm1; // previous step direction #6
    NUMBER **dirAm1; // previous step direction
    NUMBER **dirBm1; // previous step direction
    NPAR projecttosimplex; // whether projection to simplex should be done
    
    FitBit(NPAR a_nS, NPAR a_nO, NCAT a_nK, NCAT a_nG, NUMBER a_tol);
    FitBit(NPAR a_nS, NPAR a_nO, NCAT a_nK, NCAT a_nG, NUMBER a_tol, NPAR a_projecttosimplex);
    ~FitBit();
    void init(enum FIT_BIT_SLOT fbs);
    void linkPar(NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B);
    void toZero(enum FIT_BIT_SLOT fbs);
    void destroy(enum FIT_BIT_SLOT fbs);
    void copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    void add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs);
    bool checkConvergence();
    void doLog10ScaleGentle(enum FIT_BIT_SLOT fbs);
private:
    NUMBER tol;

    void init(NUMBER* &PI, NUMBER** &A, NUMBER** &B);
    void toZero(NUMBER *PI, NUMBER **A, NUMBER **B);
    void destroy(NUMBER* &PI, NUMBER** &A, NUMBER** &B);
    void get(enum FIT_BIT_SLOT fbs, NUMBER* &a_PI, NUMBER** &a_A, NUMBER** &a_B);
    void add(NUMBER *soursePI, NUMBER **sourseA, NUMBER **sourseB, NUMBER *targetPI, NUMBER **targetA, NUMBER **targetB);
    void copy(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB);
};

#endif /* defined(__HMM__FitBit__) */
