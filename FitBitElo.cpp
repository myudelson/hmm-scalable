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

#include "FitBitElo.h"

FitBitElo::FitBitElo(/*NPAR a_nS, NPAR a_nO, NCAT a_nK,*/ NCAT a_nG, NPAR a_nELO, NUMBER a_tol, NPAR a_tol_mode) {
//    this->nS = a_nS;
//    this->nO = a_nO;
    this->nG = a_nG;
//    this->nK = a_nK;
    this->nELO = a_nELO;
    this->tol = a_tol;
    this->tol_mode = a_tol_mode;
    this->ELO = NULL;
    this->ELOm1 = NULL;
    this->ELOm2 = NULL;
    this->gradELO = NULL;
    this->gradELOm1 = NULL;
    this->gradELOcopy = NULL;
    this->ELOcopy = NULL;
    this->dirELO = NULL;
    this->dirELOm1 = NULL;
    this->param = NULL;
    this->tag = 0;
    this->elo_track_g = NULL; //init1D<NUMBER>((NDAT)this->nG);
    this->elo_count_g = NULL; //init1D<NCAT>((NDAT)this->nG);
    this->elo_grad_error_g = init1D<NUMBER>((NDAT)this->nG);
    this->ll = 0.0;
}

FitBitElo::~FitBitElo() {
    if(this->ELOm1 != NULL) free(this->ELOm1);
    if(this->ELOm2 != NULL) free(this->ELOm2);
    if(this->gradELO != NULL) free(this->gradELO);
    if(this->gradELOm1 != NULL) free(this->gradELOm1);
    if(this->gradELOcopy != NULL) free(this->gradELOcopy);
    if(this->ELOcopy != NULL) free(this->ELOcopy);
    if(this->dirELO != NULL) free(this->dirELO);
    if(this->dirELOm1 != NULL) free(this->dirELOm1);
//    if(this->elo_track_g != NULL) free(this->elo_track_g);
//    if(this->elo_count_g != NULL) free(this->elo_count_g);
    if(this->elo_grad_error_g != NULL) free(this->elo_grad_error_g);
}

void FitBitElo::init(NUMBER* &a_ELO) {
    if(a_ELO == NULL)
        a_ELO = init1D<NUMBER>((NDAT)this->nELO); // init1DNumber(this->nS);
    else
        toZero1D<NUMBER>(a_ELO, (NDAT)this->nELO);
}

void FitBitElo::negate(NUMBER* &a_ELO) {
    for(NPAR i=0; i<this->nELO; i++) a_ELO[i] = -a_ELO[i];
}

void FitBitElo::link(NUMBER *a_ELO, NUMBER *a_elo_track_g, NCAT *a_elo_count_g, NUMBER *a_elo_track_g_t, struct param *a_param) {
    this->ELO = a_ELO;
    this->param = a_param;
    this->elo_track_g = a_elo_track_g;
    this->elo_count_g = a_elo_count_g;
    this->elo_track_g_t = a_elo_track_g_t;
}

void FitBitElo::toZero(NUMBER *a_ELO) {
    if(a_ELO != NULL) toZero1D<NUMBER>(a_ELO, (NDAT)this->nELO);
}

void FitBitElo::copy(NUMBER* &sourseELO, NUMBER* &targetELO){
    cpy1D<NUMBER>(sourseELO, targetELO, (NDAT)this->nELO);
}

void FitBitElo::add(NUMBER *sourseELO, NUMBER *targetELO){
    add1DNumbersWeighted(sourseELO, targetELO, this->nELO, 1.0);
}

void FitBitElo::destroy(NUMBER* &a_ELO) {
    if(a_ELO != NULL) free(a_ELO);
    a_ELO = NULL;
}

void FitBitElo::init(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            init(this->ELO);
            break;
        case FBS_PARm1:
            init(this->ELOm1);
            break;
        case FBS_PARm2:
            init(this->ELOm2);
            break;
        case FBS_GRAD:
            init(this->gradELO);
            break;
        case FBS_GRADm1:
            init(this->gradELOm1);
            break;
        case FBS_GRADcopy:
            init(this->gradELOcopy);
            break;
        case FBS_PARcopy:
            init(this->ELOcopy);
            break;
        case FBS_DIR:
            init(this->dirELO);
            break;
        case FBS_DIRm1:
            init(this->dirELOm1);
            break;
        default:
            break;
    }
}

void FitBitElo::negate(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            negate(this->ELO);
            break;
        case FBS_PARm1:
            negate(this->ELOm1);
            break;
        case FBS_PARm2:
            negate(this->ELOm2);
            break;
        case FBS_GRAD:
            negate(this->gradELO);
            break;
        case FBS_GRADm1:
            negate(this->gradELOm1);
            break;
        case FBS_GRADcopy:
            negate(this->gradELOcopy);
            break;
        case FBS_PARcopy:
            negate(this->ELOcopy);
            break;
        case FBS_DIR:
            negate(this->dirELO);
            break;
        case FBS_DIRm1:
            negate(this->dirELOm1);
            break;
        default:
            break;
    }
}

void FitBitElo::toZero(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
        case FBS_PAR:
            toZero(this->ELO);
            break;
        case FBS_PARm1:
            toZero(this->ELOm1);
            break;
        case FBS_PARm2:
            toZero(this->ELOm2);
            break;
        case FBS_GRAD:
            toZero(this->gradELO);
            break;
        case FBS_GRADm1:
            toZero(this->gradELOm1);
            break;
        case FBS_GRADcopy:
            toZero(this->gradELOcopy);
            break;
        case FBS_PARcopy:
            toZero(this->ELOcopy);
            break;
        case FBS_DIR:
            toZero(this->dirELO);
            break;
        case FBS_DIRm1:
            toZero(this->dirELOm1);
            break;
        default:
            break;
    }
}

void FitBitElo::destroy(enum FIT_BIT_SLOT fbs){
    switch (fbs) {
//         case FBS_PAR:
//             destroy(this->ELO);
//             break;
         case FBS_PARm1:
             destroy(this->ELOm1);
             break;
         case FBS_PARm2:
             destroy(this->ELOm2);
             break;
         case FBS_GRAD:
             destroy(this->gradELO);
             break;
         case FBS_GRADm1:
             destroy(this->gradELOm1);
             break;
         case FBS_GRADcopy:
             destroy(this->gradELOcopy);
             break;
         case FBS_PARcopy:
             destroy(this->ELOcopy);
             break;
         case FBS_DIR:
             destroy(this->dirELO);
             break;
         case FBS_DIRm1:
             destroy(this->dirELOm1);
             break;
         default:
             break;
     }
}

void FitBitElo::get(enum FIT_BIT_SLOT fbs, NUMBER* &a_ELO) {
    switch (fbs) {
        case FBS_PAR:
            a_ELO = this->ELO;
            break;
        case FBS_PARm1:
            a_ELO = this->ELOm1;
            break;
        case FBS_PARm2:
            a_ELO = this->ELOm2;
            break;
        case FBS_GRAD:
            a_ELO = this->gradELO;
            break;
        case FBS_GRADm1:
            a_ELO = this->gradELOm1;
            break;
        case FBS_GRADcopy:
            a_ELO = this->gradELOcopy;
            break;
        case FBS_PARcopy:
            a_ELO = this->ELOcopy;
            break;
        case FBS_DIR:
            a_ELO = this->dirELO;
            break;
        case FBS_DIRm1:
            a_ELO = this->dirELOm1;
            break;
        default:
            break;
    }
}

void FitBitElo::copy(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs) {
    NUMBER *sourseELO = NULL;
    get(sourse_fbs, sourseELO);
    NUMBER *targetELO = NULL;
    get(target_fbs, targetELO);
    
    copy(sourseELO, targetELO);
}

void FitBitElo::add(enum FIT_BIT_SLOT sourse_fbs, enum FIT_BIT_SLOT target_fbs) {
    NUMBER *sourseELO = NULL;
    get(sourse_fbs, sourseELO);
    NUMBER *targetELO = NULL;
    get(target_fbs, targetELO);

    add(sourseELO, targetELO);
}

bool FitBitElo::checkConvergence(FitResult *fr) {
    
    NUMBER criterion = 0;
    NUMBER criterion_oscil = 0; // oscillation between PAR and PARm1, i.e. PAR is close to PARm2
    bool res = false;
    switch (this->tol_mode) {
        case 'p':
            for(NPAR i=0; i<this->nELO; i++)
            {
                criterion += pow( this->ELO[i]-this->ELOm1[i], 2 )/*:0*/;
            }
            
            if( (!(sqrt(criterion) < this->tol)) && ( this->ELOm2 != NULL) ) {
                for(NPAR i=0; i<this->nELO; i++)
                {
                    criterion_oscil += pow( this->ELO[i]-this->ELOm2[i], 2 )/*:0*/;
                }
                res = sqrt(criterion_oscil) < this->tol; // double the truth or false
            }
            else
                res = sqrt(criterion) < this->tol; // double the truth or false
            break;
        case 'l':
            criterion = (fr->pOmid-fr->pO)/fr->ndat;
            res = criterion < this->tol;
            break;
    }
    return res;
}
// without checking for oscillation, actually, afer copying t-1 to t-2 and t to t-1, it is used to check for oscillation
bool FitBitElo::checkConvergenceSingle(FitResult *fr) {
    
    NUMBER criterion = 0;
    bool res = false;
    switch (this->tol_mode) {
        case 'p':
            for(NPAR i=0; i<this->nELO; i++)
            {
                criterion += pow( this->ELO[i]-this->ELOm1[i], 2 )/*:0*/;
            }
            break;
        case 'l':
            criterion = (fr->pOmid-fr->pO)/fr->ndat;
            res = criterion < this->tol;
            break;
    }
    return res;
}

void FitBitElo::doLog10ScaleGentle(enum FIT_BIT_SLOT fbs) {
    // fbs - gradient or direction
    NUMBER *a_slot = NULL;
    get(fbs, a_slot);
    
    NUMBER* slot = Calloc(NUMBER, this->nELO);
    NUMBER* par = Calloc(NUMBER, this->nELO);
    for(NPAR i=0; i<this->nELO; i++) {
        par[i] = this->ELO[i];
        slot[i] = a_slot[i];
    }
    
    NUMBER scale = doLog10Scale1DGentle(slot, par, this->nELO);
    
    for(NPAR i=0; i<this->nELO; i++) {
        a_slot[i] *= scale;
    }
    free(slot);
    free(par);
}
