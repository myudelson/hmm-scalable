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

#include "utilsSt.h"
#include "FitBitSt.h"

FitBitSt::FitBitSt(struct task *task, NUMBER* a_PARAM, NUMBER* a_GRAD, NCAT a_nPARAM) {
    this->nS = task->nS;
    this->nO = task->nO;
    this->nG = task->nG;
    this->nK = task->nK;
    this->tol = task->tol;
    this->tol_mode = task->tol_mode;
    this->nPARAM = a_nPARAM;
    this->PARAM = a_PARAM;
    this->PARAMm1 = NULL;
    this->PARAMm2 = NULL;
    this->GRAD = a_GRAD;
    this->GRADm1 = NULL;
    this->GRADcopy = NULL;
    this->PARAMcopy = NULL;
    this->DIR = NULL;
    this->DIRm1 = NULL;
    this->task = task;
    this->projecttosimplex = 1;
    this->Cslice = 0;
    this->tag = 0;
    
    this->fit_results = new FitResult[this->nK];
    for(NCAT k=0; k<this->nK; k++) {
        this->fit_results[k].iter = 0;
        this->fit_results[k].pO0_prefit = -1.0; // means never assigned
        this->fit_results[k].pO0  = 0.0;
        this->fit_results[k].pO   = 0.0;
        this->fit_results[k].conv = 0;
        this->fit_results[k].ndat = 0;
    }
    this->fit_result.iter = 0;
    this->fit_result.pO0_prefit = -1.0; // means never assigned
    this->fit_result.pO0  = 0.0;
    this->fit_result.pO   = 0.0;
    this->fit_result.conv = 0;
    this->fit_result.ndat = 0;

    this->sclsmplx_offsets = NULL;// start positions of the scaled vectors in the gradient and parameter arrays
    this->sclsmplx_sizes = NULL;// lengths of the scaled vectorsin the gradient and parameter arrays
    this->sclsmplx_bound_offsets = NULL;
    this->sclsmplx_blocks = NULL;
    this->param_blocks = NULL;
    this->sclsmplx_n = 0;

    this->PARAMm1 = init1D<NUMBER>(this->nPARAM);
    this->PARAMm2 = init1D<NUMBER>(this->nPARAM);
    
    if(this->task->solver==METHOD_GD || this->task->solver==METHOD_CGD || this->task->solver==METHOD_GBB || this->task->solver==METHOD_GDL) {
        toValue1D<NUMBER>(this->GRAD, this->nPARAM, 0); // init1D<NUMBER>(this->nPARAM);
    }
    if(this->task->solver==METHOD_CGD) {
        this->DIR = init1D<NUMBER>(this->nPARAM);
        this->DIRm1 = init1D<NUMBER>(this->nPARAM);
        this->GRADm1 = init1D<NUMBER>(this->nPARAM);
    }
    if(this->task->solver==METHOD_GBB) {
        this->GRADm1 = init1D<NUMBER>(this->nPARAM);
    }
}

FitBitSt::~FitBitSt() {
    // this->PARAM is external reference, do not destroy it
    if(this->PARAMm1 != NULL) free(this->PARAMm1);
    if(this->PARAMm2 != NULL) free(this->PARAMm2);
//    if(this->GRAD != NULL) free(this->GRAD); is external reference, do not destroy it
    if(this->GRADm1 != NULL) free(this->GRADm1);
    if(this->GRADcopy != NULL) free(this->GRADcopy);
    if(this->PARAMcopy != NULL) free(this->PARAMcopy);
    if(this->DIR != NULL) free(this->DIR);
    if(this->DIRm1 != NULL) free(this->DIRm1);
    
    delete [] this->fit_results;

    if(this->sclsmplx_offsets != NULL) free(this->sclsmplx_offsets);
    if(this->sclsmplx_sizes != NULL) free(this->sclsmplx_sizes);
    if(this->sclsmplx_blocks != NULL) free(this->sclsmplx_blocks);
    if(this->param_blocks != NULL) free(this->param_blocks);
    if(this->sclsmplx_bound_offsets != NULL) free(this->sclsmplx_bound_offsets);
}

void FitBitSt::printFitResult(NCAT k) {
    if(k==-1) {
        printf(" all skills, dat %8d, iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", this->fit_result.ndat, this->fit_result.iter, this->fit_result.pO0_prefit, this->fit_result.pO, this->fit_result.conv);
    } else {
        printf("skill %5d, dat %8d, iter#%3d p(O|param)= %15.7f >> %15.7f, conv=%d\n", k, this->fit_results[k].ndat, this->fit_results[k].iter, this->fit_results[k].pO0_prefit, this->fit_results[k].pO, this->fit_results[k].conv);
    }
}

NPAR FitBitSt::checkConvergenceBit(NDAT offset, NDAT size, FitResult *fr) {
    NUMBER criterion = 0;
    NUMBER criterion_oscil = 0; // oscillation between PAR and PARm1, i.e. PAR is close to PARm2
    NCAT i;
    NPAR res = 0;
    if (fr->iter >= this->task->first_iteration_qualify) {
        if( fr->pO0<=fr->pO) { // rough, won't go up
            res = 1;
        } else {
            switch (this->tol_mode) {
                case 'p':
                    for(i=offset; i<(offset+size); i++) {
                        criterion += pow( this->PARAM[i]-this->PARAMm1[i], 2 );
                    }
                    if( (!(sqrt(criterion) < this->tol)) && ( this->PARAMm2 != NULL ) ) {
                        for(i=offset; i<(offset+size); i++) {
                            criterion_oscil += pow( this->PARAM[i]-this->PARAMm2[i], 2 )/*:0*/;
                        }
                        res = 1*(sqrt(criterion_oscil) < this->tol); // double the truth or false
                    }
                    else
                        res = 1*(sqrt(criterion) < this->tol); // double the truth or false
                    break;
                case 'l':
                    criterion = (fr->pO0-fr->pO)/fr->ndat;
                    res = 1*(criterion < this->tol);
                    break;
            }
        }
    }
    fr->conv = res;
    return res;
}

// without checking for oscillation, actually, afer copying t-1 to t-2 and t to t-1, it is used to check for oscillation
NPAR FitBitSt::checkConvergenceBitSingle(NDAT offset, NDAT size, FitResult *fr) {
    NCAT i;
    NUMBER criterion = 0;
    NPAR res = 0;
    if(fr->pO0==fr->pO) {
        res = 1;
    } else {
        switch (this->tol_mode) {
            case 'p':
                for(i=offset; i<(offset+size); i++) {
                    criterion += pow( this->PARAM[i]-this->PARAMm1[i], 2 );
                }
                 res = 1*(criterion < this->tol);
            case 'l':
                criterion = (fr->pO0-fr->pO)/fr->ndat;
                res = 1*(criterion < this->tol);
                break;
        }
    }
    fr->conv = res;
    return res;
}

NCAT FitBitSt::checkConvergence(NCAT nbits, NDAT *offsets, NDAT* sizes, NCAT* blocks, NPAR *condition){
    NCAT n_conv = 0;
    bool res;
    for(NCAT i=0; i<nbits; i++) {
        if( condition[ blocks[i] ]==1 ) {
            res = checkConvergenceBit(offsets[i], sizes[i], &this->fit_results[i]);
//            if(blocks[i]<0 || blocks[i]>this->nK) {
//                int a = 0;
//            } // DEBUG
            condition[ blocks[i] ] = !res; // set active set value
        } else {
            res = 1;
        }
        n_conv += res;
    }
    return n_conv;
}
// without checking for oscillation, actually, afer copying t-1 to t-2     bool checkConvergenceBit(NDAT offset, NDAT size, FitResult *fr);
NCAT FitBitSt::checkConvergenceSingle(NCAT nbits, NDAT *offsets, NDAT* sizes, NCAT* blocks, NPAR *condition) {
    NCAT n_conv = 0;
    bool res;
    for(NCAT i=0; i<nbits; i++) {
        if( condition[ blocks[i] ]==1 ) {
            res = checkConvergenceBitSingle(offsets[i], sizes[i], &this->fit_results[i]);
//            if(blocks[i]<0 || blocks[i]>this->nK) {
//                int a = 0;
//            } // DEBUG
            condition[ blocks[i] ] = !res; // set active set value
        } else {
            res = 1;
        }
        n_conv += res;
    }
    return n_conv;
}


// mass-scaling via log10 using pre-constructed division of parameter space into vecros of PI, rows of A and B
// if !direction – scale GRAD, if direction – scale DIR
void FitBitSt::scaleGradients(bool direction, NPAR* condition) {
    NUMBER *TOSCALE = NULL;
    if(!direction) TOSCALE = this->GRAD;
    else TOSCALE = this->DIR;
    
    if(this->sclsmplx_n<=0 || this->sclsmplx_sizes==NULL || this->sclsmplx_offsets==NULL || this->sclsmplx_blocks==NULL) {
        fprintf(stderr,"Error, gradient scaling affordances were nnot created.\n");
        exit(1);
    }
    for(NCAT i=0; i<this->sclsmplx_n; i++) {
        if( condition[ this->sclsmplx_blocks[i] ]>0 ) {
            doLog10Scale1DGentle(&TOSCALE[ this->sclsmplx_offsets[i] ],
                                 &this->PARAM[ this->sclsmplx_offsets[i] ],
                                 this->sclsmplx_sizes[i]);
        }
    }
}

// mass-projection to simplex using pre-constructed division of parameter space into vecros of PI, rows of A and B
void FitBitSt::projectToSimplex(bool hasNon01Constraints, NPAR* condition) {
    if(this->sclsmplx_n<=0 || this->sclsmplx_sizes==NULL || this->sclsmplx_offsets==NULL || this->sclsmplx_blocks==NULL) {
        fprintf(stderr,"Error, projection to simplex affordances were nnot created.\n");
        exit(1);
    };
    if(hasNon01Constraints) {
        for(NCAT i=0; i<this->sclsmplx_n; i++) {
//            if(this->sclsmplx_blocks[i]<0 || this->sclsmplx_blocks[i]>this->nK) {
//                int a = 0;
//            } // DEBUG
            if( condition[ this->sclsmplx_blocks[i] ]>0 )
                projectsimplexbounded(&this->PARAM[ this->sclsmplx_offsets[i] ],
                                      &this->task->param_values_lb[ this->sclsmplx_bound_offsets[i] ],
                                      &this->task->param_values_ub[ this->sclsmplx_bound_offsets[i] ],
                                      this->sclsmplx_sizes[i]);
        }
    } else {
        for(NCAT i=0; i<this->sclsmplx_n; i++) {
//            if(this->sclsmplx_blocks[i]<0 || this->sclsmplx_blocks[i]>this->nK) {
//                int a = 0;
//            } // DEBUG
            if( condition[ this->sclsmplx_blocks[i] ] )
                projectsimplex(&this->PARAM[ this->sclsmplx_offsets[i] ], this->sclsmplx_sizes[i]);
        }
    }
}
