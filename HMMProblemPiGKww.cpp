/*
 *  HMMProblemPiGKww.h
 *  HMM
 *  ww  - version with weighling of student and skill components (both global)
 *
 *  Created by Mikhail Yudelson on 11/3/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "HMMProblemPiGKww.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string>
#include "utils.h"
#include <math.h>
#include <map>

HMMProblemPiGKww::HMMProblemPiGKww(struct param *param) : HMMProblemPiGK(param) {
    init(param);
}

void HMMProblemPiGKww::init(struct param *param) {
//    HMMProblemPiGK::init(param);
    this->ww = init1D<NUMBER>(3);
	this->ww[0] = 0.5;//1;
	this->ww[1] = 0.5;//1;
    this->ww[2] = 0;//0;
}

HMMProblemPiGKww::~HMMProblemPiGKww() {
    // no double destroy
    if(this->ww  != NULL) free(this->ww);
}

//void HMMProblemPiGKww::destroy() {
//}// ~HMMProblemPiGKww
//
//NUMBER** HMMProblemPiGKww::getPI() { // same as getPIk
//	return this->PI;
//}

//NUMBER** HMMProblemPiGKww::getPIk() {
//	return this->PI;
//}
//
//NUMBER** HMMProblemPiGKww::getPIg() {
//	return this->PIg;
//}
//
//NUMBER* HMMProblemPiGKww::getPI(NCAT x) { // same as getPIk(x)
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->PI[x];
//}
//
//NUMBER** HMMProblemPiGKww::getA(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing A, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->A[x];
//}

//NUMBER** HMMProblemPiGKww::getB(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing B, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->B[x];
//}
//
//NUMBER* HMMProblemPiGKww::getPIk(NCAT x) {
//	if( x > (this->p->nK-1) ) {
//		fprintf(stderr,"While accessing PI_k, skill index %d exceeded last index of the data %d.\n", x, this->p->nK-1);
//		exit(1);
//	}
//	return this->PI[x];
//}
//
//NUMBER* HMMProblemPiGKww::getPIg(NCAT x) {
//	if( x > (this->p->nG-1) ) {
//		fprintf(stderr,"While accessing PI_g, skill index %d exceeded last index of the data %d.\n", x, this->p->nG-1);
//		exit(1);
//	}
//	return this->PIg[x];
//}

NUMBER HMMProblemPiGKww::getPI(struct data* dt, NPAR i) {
//    NUMBER p = this->PI[dt->k][i] * 2 * this->ww[0][i], q = this->PIg[dt->g][i] * 2 * this->ww[1][i];
//    return 1/( 1 + (1-p)*(1-q)/(p*q) );

//    NUMBER lw0 = logit(this->ww[0][i]), lw1 = logit(this->ww[1][i]);
//    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
//    return 1/( 1 + pow(1-p,lw0)*pow(1-q,lw1)/(pow(p,lw0)*pow(q,lw1)) );
    
//    NUMBER p = this->PI[dt->k][i] * 2 * this->ww[0], q = this->PIg[dt->g][i] * 2 * this->ww[1];
//    return 1/( 1 + (1-p)*(1-q)/(p*q) );

//    // a logit(k) + b logit(u)
//    NUMBER p = pow((1-this->PI[dt->k][i])/safe0num(this->PI[dt->k][i]), this->ww[0]), q = pow((1-this->PIg[dt->g][i])/safe0num(this->PIg[dt->g][i]), this->ww[1]);
//    return 1/( 1 + p*q/exp(this->ww[2]) );
    
    // logit1(a) logit(k) + logit1(b) logit(u)
    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
    return 1/( 1 + pow( (1-p)/p,1+logit(this->ww[0])) * pow( (1-q)/q,1+logit(this->ww[1])) * ( (0.5-this->ww[2])/safe0num(0.5+this->ww[2]) )/* */  );
    
//    // logit1(a) logit(k) + logit1(b) logit(u)
//    NUMBER p = this->PI[dt->k][i], q = this->PIg[dt->g][i];
//    NUMBER a = 1+sqrt(3)*(this->ww[0]-0.5)/sqrt(1-pow(this->ww[0]-0.5,2)),
//           b = 1+sqrt(3)*(this->ww[1]-0.5)/sqrt(1-pow(this->ww[1]-0.5,2)),
//           c = 1+sqrt(3)*(this->ww[2]-0.5)/sqrt(1-pow(this->ww[2]-0.5,2));
//    
//    return 1/( 1 + pow( (1-p)/p,a) * pow( (1-q)/q,b) * exp(-c)  );
    
    
    
}

//// getters for computing alpha, beta, gamma
//NUMBER HMMProblemPiGKww::getA(struct data* dt, NPAR i, NPAR j) {
//    return this->A[dt->k][i][j];
//}
//
//// getters for computing alpha, beta, gamma
//NUMBER HMMProblemPiGKww::getB(struct data* dt, NPAR i, NPAR m) {
//    // special attention for "unknonw" observations, i.e. the observation was there but we do not know what it is
//    // in this case we simply return 1, effectively resulting in no change in \alpha or \beta vatiables
//    if(m<0)
//        return 1;
//    return this->B[dt->k][i][m];
//}

void HMMProblemPiGKww::setGradPI(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit;
    //    o = dt->obs[t];
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        o = this->p->dat_obs->get( dt->ix[t] );
    
//    // a logit(k) + b logit(u)
//    if(kg_flag == 0) { // k
//        for(i=0; i<this->p->nS; i++) {
//            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
//            deriv_logit = this->ww[0] / safe0num( this->PI[ dt->k ][i] * ( 1 - this->PI[ dt->k ][i] ) );
//			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,this->PI[ dt->k ][i], 0.5); // PENALTY;
//        }
//    }
//    else // g
//        for(i=0; i<this->p->nS; i++) {
//            combined = getPI(dt,i);
//            deriv_logit = this->ww[1] / safe0num( this->PIg[ dt->g ][i] * ( 1 - this->PIg[ dt->g ][i] ) );
//			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,this->PIg[ dt->g ][i], 0.5); // PENALTY;
//        }
    
    // logit1(a) logit(k) + logit1(b) logit(u)
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            deriv_logit = (1+logit(this->ww[0])) / safe0num( fb->PI[i] * ( 1 - fb->PI[i] ) );
			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) +
                L2penalty(this->p,fb->PI[i], 0.5); // PENALTY;
        }
    }
    
//    // x / sqrt(1+ x^2)
//    if(kg_flag == 0) { // k
//        for(i=0; i<this->p->nS; i++) {
//            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
//            deriv_logit = (1+sqrt(3)*(this->ww[0]-0.5)/sqrt(1-pow(this->ww[0]-0.5,2))) / safe0num( this->PI[ dt->k ][i] * ( 1 - this->PI[ dt->k ][i] ) );
//			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,this->PI[ dt->k ][i], 0.5); // PENALTY;
//        }
//    }
//    else // g
//        for(i=0; i<this->p->nS; i++) {
//            combined = getPI(dt,i);
//            deriv_logit = (1+sqrt(3)*(this->ww[1]-0.5)/sqrt(1-pow(this->ww[1]-0.5,2))) / safe0num( this->PIg[ dt->g ][i] * ( 1 - this->PIg[ dt->g ][i] ) );
//			fb->gradPI[i] -= combined * (1-combined) * deriv_logit * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param) + L2penalty(this->p,this->PIg[ dt->g ][i], 0.5); // PENALTY;
//        }
}

void HMMProblemPiGKww::setGradWW(FitBit *fb){
    if(this->p->block_fitting[0]>0) return;
    NDAT t = 0;
    NPAR i, o;
    NUMBER combined, deriv_logit0, deriv_logit1, deriv_logit2;
    //    o = dt->obs[t];
    struct data* dt;
    for(NCAT x=0; x<fb->xndat; x++) {
        dt = fb->x_data[x];
        if( dt->cnt!=0 ) continue;
        o = this->p->dat_obs->get( dt->ix[t] );
        for(i=0; i<fb->nS; i++) {
            combined = getPI(dt,i);//sigmoid( logit(this->PI[k][i]) + logit(this->PIg[g][i]) );
            
//            // a logit(k) + b logit(u)
//            fb->gradPI[0] -= combined * (1-combined) * logit( this->PI [ dt->k ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
//            fb->gradPI[1] -= combined * (1-combined) * logit( this->PIg[ dt->g ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
//            fb->gradPI[2] -= combined * (1-combined) * 1                              * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
            
            
            // logit1(a) logit(k) + logit1(b) logit(u)
            deriv_logit0 = 1/ safe0num( this->ww[0] * (1-this->ww[0]) );
            deriv_logit1 = 1/ safe0num( this->ww[1] * (1-this->ww[1]) );
            deriv_logit2 = 1/ safe0num( -0.25 + this->ww[2] * this->ww[2] );
            fb->gradPI[0] -= combined * (1-combined) * deriv_logit0 * logit( this->PI [ dt->k ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
            fb->gradPI[1] -= combined * (1-combined) * deriv_logit1 * logit( this->PIg[ dt->g ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
//            fb->gradPI[2] -= combined * (1-combined) * deriv_logit2 * 1                              * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
            
//            // x / sqrt( 1 + x^2 )
//            deriv_logit0 = sqrt(3) / safe0num(pow( 1 - pow(this->ww[0]-0.5,2), 1.5 ));
//            deriv_logit1 = sqrt(3) / safe0num(pow( 1 - pow(this->ww[1]-0.5,2), 1.5 ));
////            deriv_logit2 = sqrt(3) / safe0num(pow( 1 - pow(this->ww[2]-0.5,2), 1.5 ));
//            fb->gradPI[0] -= combined * (1-combined) * deriv_logit0 * logit( this->PI [ dt->k ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
//            fb->gradPI[1] -= combined * (1-combined) * deriv_logit1 * logit( this->PIg[ dt->g ][i] ) * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
////            fb->gradPI[2] -= combined * (1-combined) * deriv_logit2 * 1                              * dt->beta[t][i] * ((o<0)?1:getB(dt,i,o)) / safe0num(dt->p_O_param);
        }
    }
}

void HMMProblemPiGKww::toFile(const char *filename) {
	FILE *fid = fopen(filename,"w");
	if(fid == NULL) {
		fprintf(stderr,"Can't write output model file %s\n",filename);
		exit(1);
	}
    
    // write solved id
    writeSolverInfo(fid, this->p);
    
	fprintf(fid,"Null skill ratios\t");
	for(NPAR m=0; m<this->p->nO; m++)
		fprintf(fid," %10.7f%s",this->null_obs_ratio[m],(m==(this->p->nO-1))?"\n":"\t");
    /*vvdiffvv*/
	fprintf(fid,"Student and skill component ratios\t");
	fprintf(fid,"%10.7f\t%10.7f\t%10.7f\n",this->ww[0],this->ww[1],this->ww[2]);
//	fprintf(fid,"%10.7f\t%10.7f\n",this->ww[0],this->ww[1]);
    /*^^diff^^*/
	NCAT k, g;
    NPAR i,j,m;
	std::map<NCAT,std::string>::iterator it;
	for(g=0;g<this->p->nG; g++) {
		it = this->p->map_group_bwd->find(g);
		fprintf(fid,"%d\t%s\n",g,it->second.c_str());
		fprintf(fid,"PIg\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PIg[g][i],(i==(this->p->nS-1))?"\n":"\t");
    }
	for(k=0;k<this->p->nK;k++) {
		it = this->p->map_skill_bwd->find(k);
		fprintf(fid,"%d\t%s\n",k,it->second.c_str());
		fprintf(fid,"PIk\t");
		for(i=0; i<this->p->nS; i++)
			fprintf(fid,"%10.8f%s",this->PI[k][i],(i==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++)
				fprintf(fid,"%10.8f%s",this->A[k][i][j],(i==(this->p->nS-1) && j==(this->p->nS-1))?"\n":"\t");
		fprintf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nO; m++)
				fprintf(fid,"%10.8f%s",this->B[k][i][m],(i==(this->p->nS-1) && m==(this->p->nO-1))?"\n":"\t");
	}
	fclose(fid);
}

void HMMProblemPiGKww::fit() {
    NUMBER* loglik_rmse = init1D<NUMBER>(2);
    FitNullSkill(loglik_rmse, false /*do not need RMSE/SE*/);
    if(this->p->structure == STRUCTURE_PIgkww) {
        loglik_rmse[0] += GradientDescent();
    } else {
        fprintf(stderr,"Solver specified is not supported.\n");
    }
    this->neg_log_lik = loglik_rmse[0];
    free(loglik_rmse);
}

NDAT HMMProblemPiGKww::computeGradients(FitBit *fb){
    fb->toZero(FBS_GRAD);
    
    NDAT ndat = computeAlphaAndPOParam(fb->xndat, fb->x_data);
    computeBeta(fb->xndat, fb->x_data);

    /*vv diff vv*/
    if(fb->PI != NULL && this->p->block_fitting[0]==0) setGradPI(fb);
    if(fb->PI != NULL && this->p->block_fitting[0]==0) setGradWW(fb);
    /*^^ diff ^^*/

    if(fb->A  != NULL && this->p->block_fitting[1]==0) setGradA(fb);
    if(fb->B  != NULL && this->p->block_fitting[2]==0) setGradB(fb);
    return ndat;
} // computeGradients()

NUMBER HMMProblemPiGKww::GradientDescent() {
	NCAT k, g, x;
    NPAR nS = this->p->nS, nO = this->p->nO; NCAT nK = this->p->nK, nG = this->p->nG;
    NUMBER loglik = 0;
    FitResult fr;
    FitBit *fb = new FitBit(this->p->nS, this->p->nO, this->p->nK, this->p->nG, this->p->tol);
    fb->init(FBS_PARm1);
    fb->init(FBS_GRAD);
    if(this->p->solver==METHOD_CGD) {
        fb->init(FBS_GRADm1);
        fb->init(FBS_DIRm1);
    }

//    // all together
//    NUMBER ** copyPI = init2D<NUMBER>((NDAT)this->sizes[0], (NDAT)nS);
//    NUMBER ** copyPIg = init2D<NUMBER>(nG, (NDAT)nS);
//    NUMBER *** copyA =  init3D<NUMBER>((NDAT)this->sizes[1], (NDAT)nS, (NDAT)nS);
//    NUMBER *** copyB =  init3D<NUMBER>((NDAT)this->sizes[2], (NDAT)nS, (NDAT)nO);

    // fit bit and fir result for ww
    FitResult frww;
    FitBit *fbww = new FitBit(3, nO, nK, nG, this->p->tol, 1/*do project to simplex*/);
    fbww->init(FBS_PARm1);
    fbww->init(FBS_GRAD);
    if(this->p->solver==METHOD_CGD) {
        fbww->init(FBS_GRADm1);
        fbww->init(FBS_DIRm1);
    }
    
    //
	// fit all as 1 skill first, set group gradients to 0, and do not fit them
	//
	if(this->p->single_skill>0) {
        fb->link( HMMProblem::getPI(0), HMMProblem::getA(0), HMMProblem::getB(0), this->p->nSeq, this->p->k_data);// link skill 0 (we'll copy fit parameters to others
        NCAT* original_ks = Calloc(NCAT, this->p->nSeq);
        for(x=0; x<this->p->nSeq; x++) { original_ks[x] = this->p->all_data[x].k; this->p->all_data[x].k = 0; } // save progonal k's
        fr = GradientDescentBit(fb);
        for(x=0; x<this->p->nSeq; x++) { this->p->all_data[x].k = original_ks[x]; } // restore original k's
        free(original_ks);
        if( !this->p->quiet )
            printf("single skill iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n", fr.iter,fr.pO0,fr.pO,fr.conv);
    }
    
	//
	// Main fit
	//
    if( this->p->single_skill!=2 ) { // if not "force single skill"
        NCAT first_iteration_qualify = this->p->first_iteration_qualify; // at what iteration, qualification for skill/group convergence should start
        NCAT iterations_to_qualify = this->p->iterations_to_qualify; // how many concecutive iterations necessary for skill/group to qualify as converged
        NCAT* iter_qual_skill = Calloc(NCAT, (size_t)nK);
        NCAT* iter_qual_group = Calloc(NCAT, (size_t)nG);
        int skip_k = 0, skip_g = 0;
        
        
        int i = 0; // count runs
        while(skip_k<nK || skip_g<nG || frww.conv==0) {
            //
            // Skills first
            //
            for(k=0; k<nK && skip_k<nK; k++) { // for all A,B-by-skill
                if(iter_qual_skill[k]==iterations_to_qualify)
                    continue;
//                NCAT xndat = this->p->k_numg[k];
//                struct data** x_data = this->p->k_g_data[k];
                // link and fit
                fb->link( HMMProblem::getPI(k), HMMProblem::getA(k), HMMProblem::getB(k), this->p->k_numg[k], this->p->k_g_data[k]);// link skill 0 (we'll copy fit parameters to others
//                cpy1D<NUMBER>(this->PI[k],copyPI[k],nS); /*prep hide*/
                fr = GradientDescentBit(fb);
                // decide on convergence
                if(i>=first_iteration_qualify) {
                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_g==nG) { // converged quick, or don't care (others all converged
                        iter_qual_skill[k]++;
                        if(iter_qual_skill[k]==iterations_to_qualify || skip_g==nG) {// criterion met, or don't care (others all converged)
                            if(skip_g==nG) iter_qual_skill[k]=iterations_to_qualify; // G not changing anymore
                            skip_k++;
                            if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                computeAlphaAndPOParam(fb->xndat, fb->x_data);
                                printf("run %2d skipK %4d skill %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_k,k,fr.iter,fr.pO0,fr.pO,fr.conv);
                            }
                        }
                    }
                    else
                        iter_qual_skill[k]=0;
                } // decide on convergence
//                swap1D<NUMBER>(this->PI[k],copyPI[k],nS); /*hide*/
            } // for all skills
            //
            // PIg second
            //
            for(g=0; g<nG && skip_g<nG; g++) { // for all PI-by-user
                if(iter_qual_group[g]==iterations_to_qualify)
                    continue;
//                NCAT xndat = this->p->g_numk[g];
//                struct data** x_data = this->p->g_k_data[g];
//                cpy1D<NUMBER>(this->PIg[g],copyPIg[g],nS); /*prep hide*/
                // vvvvvvvvvvvvvvvvvvvvv ONLY PART THAT IS DIFFERENT FROM others
                fb->link(this->getPIg(g), NULL, NULL, this->p->g_numk[g], this->p->g_k_data[g]);
                // ^^^^^^^^^^^^^^^^^^^^^
                // decide on convergence
                fr = GradientDescentBit(fb);
                if(i>=first_iteration_qualify) {
                    if(fr.iter==1 /*e<=this->p->tol*/ || skip_k==nK) { // converged quick, or don't care (others all converged
                        iter_qual_group[g]++;
                        if(iter_qual_group[g]==iterations_to_qualify || skip_k==nK) {// criterion met, or don't care (others all converged)
                            if(skip_k==nK) iter_qual_group[g]=iterations_to_qualify; // K not changing anymore
                            skip_g++;
                            if( !this->p->quiet && ( /*(!conv && iter<this->p->maxiter) ||*/ (fr.conv || fr.iter==this->p->maxiter) )) {
                                computeAlphaAndPOParam(fb->xndat, fb->x_data);
                                printf("run %2d skipG %4d group %4d iter#%3d p(O|param)= %15.7f -> %15.7f, conv=%d\n",i,skip_g,g,fr.iter,fr.pO0,fr.pO,fr.conv);
                            }
                        }
                    }
                    else
                        iter_qual_group[g]=0;
                } // decide on convergence
//                swap1D<NUMBER>(this->PIg[g],copyPIg[g],nS); /*hide*/
            } // for all groups
            
//            // vvv temporary explorative
//            NCAT xndat = this->p->nSeq;
//            struct data** x_data = this->p->k_data;
//            
//            this->ww[0]=1; this->ww[1]=1; this->ww[2]=0;
//            computeAlphaAndPOParam(xndat, x_data);
//            NUMBER l1 = loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
//            
//            this->ww[0]=1.01; this->ww[1]=1; this->ww[2]=0;
//            computeAlphaAndPOParam(xndat, x_data);
//            NUMBER l2 = loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
//
//            this->ww[0]=1.05; this->ww[1]=1; this->ww[2]=0;
//            computeAlphaAndPOParam(xndat, x_data);
//            NUMBER l3 = loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
//            // ^^^ temporary explorative
//            
            
            
            //
            // ww third
            //
//            NCAT xndat = this->p->nSeq;
//            struct data** x_data = this->p->k_data;
            fbww->link(this->ww, NULL, NULL, this->p->nSeq, this->p->k_data);
            computeAlphaAndPOParam(fb->xndat, fb->x_data);
//            loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
            bool non01constraints = this->non01constraints;
            this->non01constraints = false;
            frww = GradientDescentBit(fbww);
            this->non01constraints = non01constraints;
//            loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
            
//            // unhide
//            cpy2D<NUMBER>(copyPI,this->PI,nK,nS);
//            cpy2D<NUMBER>(copyPIg,this->PIg,nG,nS);
//            cpy3D<NUMBER>(copyA,this->A,nK,nS,nS);
//            cpy3D<NUMBER>(copyB,this->B,nK,nS,nS);
            i++;
        }
        // recycle qualifications
        if( iter_qual_skill != NULL ) free(iter_qual_skill);
        if( iter_qual_group != NULL) free(iter_qual_group);
    } // if not "force single skill
    
    delete fb;
    delete fbww;
    // compute loglik
    loglik = getSumLogPOPara(this->p->nSeq, this->p->k_data);
    return loglik;

//    // all together
//    if(copyPI!=NULL) free2D(copyPI,(NDAT)this->sizes[0]);
//    if(copyPIg!=NULL) free2D(copyPIg,(NDAT)nG);
//    if(copyA !=NULL) free3D(copyA, (NDAT)this->sizes[1], (NDAT)nS);
//    if(copyB !=NULL) free3D(copyB, (NDAT)this->sizes[2], (NDAT)nS);
}

void HMMProblemPiGKww::readModelBody(FILE *fid, struct param* param, NDAT *line_no, bool overwrite) {
	NPAR i,j,m;
	NCAT k = 0, g = 0, idxk = 0, idxg = 0;
	string s;
    std::map<std::string,NCAT>::iterator it;
    char col[2048];
    // null skill ratios
    readNullObsRatio(fid, param, line_no);
    // ww weights
    fscanf(fid, "Student and skill component ratios\t");
    this->ww[0] = 0.5;//1;
    this->ww[1] = 0.5;//1;
    this->ww[2] = 0;//0;
    fscanf(fid,"%lf\t",&this->ww[0] );
    fscanf(fid,"%lf\t",&this->ww[1] );
    fscanf(fid,"%lf\n",&this->ww[2] );
    (*line_no)++;
    //
    // init param
    //
    if(overwrite) {
        this->p->map_group_fwd = new map<string,NCAT>();
        this->p->map_group_bwd = new map<NCAT,string>();
        this->p->map_skill_fwd = new map<string,NCAT>();
        this->p->map_skill_bwd = new map<NCAT,string>();
    }
	//
	// read grouped PIg
	//
    for(g=0; g<this->p->nG; g++) {
		// read group label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;        
        if(overwrite) {
            this->p->map_group_fwd->insert(pair<string,NCAT>(s, (NCAT)this->p->map_group_fwd->size()));
            this->p->map_group_bwd->insert(pair<NCAT,string>((NCAT)this->p->map_group_bwd->size(), s));
            idxg = g;
        } else {
            it = this->p->map_group_fwd->find(s);
            if( it==this->p->map_group_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                continue; // skip this iteration
            }
            else
                idxg =it->second;
        }

        // read PIg
        fscanf(fid,"PIg\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->PIg[idxg][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->PIg[idxg][i] = atof(col);
        (*line_no)++;
    }
    //
    // read skills
    //
	for(k=0; k<this->p->nK; k++) {
		// read skill label
        fscanf(fid,"%*s\t%[^\n]\n",col);
        s = string( col );
        (*line_no)++;
        if(overwrite) {
            this->p->map_skill_fwd->insert(pair<string,NCAT>(s, (NCAT)this->p->map_skill_fwd->size()));
            this->p->map_skill_bwd->insert(pair<NCAT,string>((NCAT)this->p->map_skill_bwd->size(), s));
            idxk = k;
        } else {
            it = this->p->map_skill_fwd->find(s);
            if( it==this->p->map_skill_fwd->end() ) { // not found, skip 3 lines and continue
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                fscanf(fid,"%*s\n");
                continue; // skip this iteration
            }
            else
                idxk =it->second;
        }
        
        // read PI
        fscanf(fid,"PIk\t");
        for(i=0; i<(this->p->nS-1); i++) { // read 1 less then necessary
            fscanf(fid,"%[^\t]\t",col);
            this->PI[idxk][i] = atof(col);
        }
        fscanf(fid,"%[^\n]\n",col);// read last one
        this->PI[idxk][i] = atof(col);
        (*line_no)++;
		// read A
        fscanf(fid,"A\t");
		for(i=0; i<this->p->nS; i++)
			for(j=0; j<this->p->nS; j++) {
                if(i==(this->p->nS-1) && j==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->A[idxk][i][j] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->A[idxk][i][j] = atof(col);
                }
			}
        (*line_no)++;
		// read B
        fscanf(fid,"B\t");
		for(i=0; i<this->p->nS; i++)
			for(m=0; m<this->p->nS; m++) {
                if(i==(this->p->nS-1) && m==(this->p->nS-1)) {
                    fscanf(fid,"%[^\n]\n", col); // last one;
                    this->B[idxk][i][m] = atof(col);
                }
                else {
                    fscanf(fid,"%[^\t]\t", col); // not las one
                    this->B[idxk][i][m] = atof(col);
                }
			}
        (*line_no)++;
	} // for all k
}
