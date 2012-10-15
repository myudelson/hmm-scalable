//
//  FitBit.cpp
//  HMM
//
//  Created by Mikhail Yudelson on 10/12/12.
//
//

#include "FitBit.h"

FitBit::FitBit(struct param *p) {
    this->nS = p->nO;
    this->nO = p->nS;
    this->nG = p->nG;
    this->nK = p->nK;
    this->PI = NULL;
    this->A = NULL;
    this->B = NULL;
    this->PIm1 = NULL;
    this->Am1 = NULL;
    this->Bm1 = NULL;
    this->gradPI = NULL;
    this->gradA = NULL;
    this->gradB = NULL;
    this->gradPIm1 = NULL;
    this->gradAm1 = NULL;
    this->gradBm1 = NULL;
    this->gradPIsum = NULL;
    this->gradAsum = NULL;
    this->gradBsum = NULL;
    this->gradPIsumm1 = NULL;
    this->gradAsumm1 = NULL;
    this->gradBsumm1 = NULL;
    this->dirPIm1 = NULL;
    this->dirAm1 = NULL;
    this->dirBm1 = NULL;
}

FitBit::~FitBit() {
//    if(this->PI != NULL) free(this->PI); // these are usually linked
//    if(this->A != NULL) free2DNumber(this->A, this->nS); // these are usually linked
//    if(this->B != NULL) free2DNumber(this->B, this->nS); // these are usually linked
    if(this->PIm1 != NULL) free(this->PIm1);
    if(this->Am1 != NULL) free2DNumber(this->Am1, this->nS);
    if(this->Bm1 != NULL) free2DNumber(this->Bm1, this->nS);
    if(this->gradPI != NULL) free(this->gradPI);
    if(this->gradA != NULL) free2DNumber(this->gradA, this->nS);
    if(this->gradB != NULL) free2DNumber(this->gradB, this->nS);
    if(this->gradPIm1 != NULL) free(this->gradPIm1);
    if(this->gradAm1 != NULL) free2DNumber(this->gradAm1, this->nS);
    if(this->gradBm1 != NULL) free2DNumber(this->gradBm1, this->nS);
    if(this->gradPIsum != NULL) free(this->gradPIsum);
    if(this->gradAsum != NULL) free2DNumber(this->gradAsum, this->nS);
    if(this->gradBsum != NULL) free2DNumber(this->gradBsum, this->nS);
    if(this->gradPIsumm1 != NULL) free(this->gradPIsumm1);
    if(this->gradAsumm1 != NULL) free2DNumber(this->gradAsumm1, this->nS);
    if(this->gradBsumm1 != NULL) free2DNumber(this->gradBsumm1, this->nS);
    if(this->dirPIm1 != NULL) free(this->dirPIm1);
    if(this->dirAm1 != NULL) free2DNumber(this->dirAm1, this->nS);
    if(this->dirBm1 != NULL) free2DNumber(this->dirBm1, this->nS);
}

void FitBit::init(NUMBER* &a_PI, NUMBER** &a_A, NUMBER** &a_B) {
    a_PI = init1DNumber(this->nS);
    a_A  = init2DNumber(this->nS, this->nS);
    a_B  = init2DNumber(this->nS, this->nO);
}

void FitBit::linkPar(NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B) {
    this->PI = a_PI;
    this->A  = a_A;
    this->B  = a_B;
}

void FitBit::toZero(NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B) {
    toZero1DNumber(a_PI, this->nS);
    toZero2DNumber(a_A,  this->nS, this->nS);
    toZero2DNumber(a_B,  this->nS, this->nO);
}

void FitBit::copy(NUMBER* &soursePI, NUMBER** &sourseA, NUMBER** &sourseB, NUMBER* &targetPI, NUMBER** &targetA, NUMBER** &targetB){
    cpy1DNumber(soursePI, targetPI, this->nS);
    cpy2DNumber(sourseA,  targetA,  this->nS, this->nS);
    cpy2DNumber(sourseB,  targetB,  this->nS, this->nO);
}

void FitBit::add(NUMBER *soursePI, NUMBER **sourseA, NUMBER **sourseB, NUMBER *targetPI, NUMBER **targetA, NUMBER **targetB){
    add1DNumbersWeighted(soursePI, targetPI, this->nS, 1.0);
    add2DNumbersWeighted(sourseA,  targetA,  this->nS, this->nS, 1.0);
    add2DNumbersWeighted(sourseB,  targetB,  this->nS, this->nO, 1.0);
}

void FitBit::destroy(NUMBER* &a_PI, NUMBER** &a_A, NUMBER** &a_B) {
    free(PI);
    free2DNumber(A, this->nS);
    free2DNumber(B, this->nS);
    PI = NULL;
    A  = NULL;
    B  = NULL;
}

void FitBit::init(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            init(this->PI, this->A, this->B);
            break;
        case FBS_PARm1:
            init(this->PIm1, this->Am1, this->Bm1);
            break;
        case FBS_GRAD:
            init(this->gradPI, this->gradA, this->gradB);
            break;
        case FBS_GRADm1:
            init(this->gradPIm1, this->gradAm1, this->gradBm1);
            break;
        case FBS_GRADsum:
            init(this->gradPIsum, this->gradAsum, this->gradBsum);
            break;
        case FBS_GRADsumm1:
            init(this->gradPIsumm1, this->gradAsumm1, this->gradBsumm1);
            break;
        case FBS_DIRm1:
            init(this->dirPIm1, this->dirAm1, this->dirBm1);
            break;
        default:
            break;
    }
}

void FitBit::toZero(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            toZero(this->PI, this->A, this->B);
            break;
        case FBS_PARm1:
            toZero(this->PIm1, this->Am1, this->Bm1);
            break;
        case FBS_GRAD:
            toZero(this->gradPI, this->gradA, this->gradB);
            break;
        case FBS_GRADm1:
            toZero(this->gradPIm1, this->gradAm1, this->gradBm1);
            break;
        case FBS_GRADsum:
            toZero(this->gradPIsum, this->gradAsum, this->gradBsum);
            break;
        case FBS_GRADsumm1:
            toZero(this->gradPIsumm1, this->gradAsumm1, this->gradBsumm1);
            break;
        case FBS_DIRm1:
            toZero(this->dirPIm1, this->dirAm1, this->dirBm1);
            break;
        default:
            break;
    }
}

void FitBit::destroy(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            destroy(this->PI, this->A, this->B);
            break;
        case FBS_PARm1:
            destroy(this->PIm1, this->Am1, this->Bm1);
            break;
        case FBS_GRAD:
            destroy(this->gradPI, this->gradA, this->gradB);
            break;
        case FBS_GRADm1:
            destroy(this->gradPIm1, this->gradAm1, this->gradBm1);
            break;
        case FBS_GRADsum:
            destroy(this->gradPIsum, this->gradAsum, this->gradBsum);
            break;
        case FBS_GRADsumm1:
            destroy(this->gradPIsumm1, this->gradAsumm1, this->gradBsumm1);
            break;
        case FBS_DIRm1:
            destroy(this->dirPIm1, this->dirAm1, this->dirBm1);
            break;
        default:
            break;
    }
}

void FitBit::get(enum FIT_BIT_SLOT fbs, NUMBER* &a_PI, NUMBER** &a_A, NUMBER** &a_B) {
    switch (fbs) {
        case FBS_PAR:
            a_PI = this->PI;
            a_A  = this->A;
            a_B  = this->B;
            break;
        case FBS_PARm1:
            a_PI = this->PIm1;
            a_A  = this->Am1;
            a_B  = this->Bm1;
            break;
        case FBS_GRAD:
            a_PI = this->gradPI;
            a_A  = this->gradA;
            a_B  = this->gradB;
            break;
        case FBS_GRADm1:
            a_PI = this->gradPIm1;
            a_A  = this->gradAm1;
            a_B  = this->gradBm1;
            break;
        case FBS_GRADsum:
            a_PI = this->gradPIsum;
            a_A  = this->gradAsum;
            a_B  = this->gradBsum;
            break;
        case FBS_GRADsumm1:
            a_PI = this->gradPIsumm1;
            a_A  = this->gradAsumm1;
            a_B  = this->gradBsumm1;
            break;
        case FBS_DIRm1:
            a_PI = this->dirPIm1;
            a_A  = this->dirAm1;
            a_B  = this->dirBm1;
            break;
        default:
            break;
    }
}

void FitBit::copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs) {
    NUMBER *soursePI;
    NUMBER **sourseA;
    NUMBER **sourseB;
    get(sourse_fbs, soursePI, sourseA, sourseB);
    NUMBER *targetPI;
    NUMBER **targetA;
    NUMBER **targetB;
    get(target_fbs, targetPI, targetA, targetB);
    
    copy(soursePI, sourseA, sourseB, targetPI, targetA, targetB);
}

void FitBit::add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs) {
    NUMBER *soursePI;
    NUMBER **sourseA;
    NUMBER **sourseB;
    get(sourse_fbs, soursePI, sourseA, sourseB);
    NUMBER *targetPI;
    NUMBER **targetA;
    NUMBER **targetB;
    get(target_fbs, targetPI, targetA, targetB);
    
    add(soursePI, sourseA, sourseB, targetPI, targetA, targetB);
}

bool FitBit::checkConvergence(struct param *p, bool flags[3]) {
	NUMBER critetion = 0;
	for(NPAR i=0; i<this->nS; i++)
	{
		if(flags[0]) critetion += pow( this->PI[i]-this->PIm1[i], 2 )/*:0*/;
		for(NPAR j=0; flags[1] && j<this->nS; j++) {
			critetion += pow(this->A[i][j] - this->Am1[i][j],2);
		}
		for(NPAR k=0; flags[2] && k<this->nO; k++) {
			critetion += pow(this->B[i][k] - this->Bm1[i][k],2);
		}
	}
	return sqrt(critetion) < p->tol;
}


