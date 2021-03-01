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
#include "StripedArray.h"
#include "HMMProblemSt.h"

//#include <boost/numeric/ublas/matrix_sparse.hpp>//BOOST
//#include <boost/numeric/ublas/io.hpp>//BOOST

#ifndef _HMMPROBLEMELOST_H
#define _HMMPROBLEMELOST_H

class HMMProblemEloSt : public HMMProblemSt {
public:
    HMMProblemEloSt();
    HMMProblemEloSt(struct task* task); // sizes=={nK, nK, nK} by default
    virtual ~HMMProblemEloSt();
    // getters for parameters in context
    virtual NUMBER getPI(struct context* ctx, NPAR i);
    virtual NUMBER getA (struct context* ctx, NPAR i, NPAR j);
    virtual NUMBER getB (struct context* ctx, NPAR i, NPAR m);
    // getters for computing gradients of alpha, beta, gamma
//	virtual void toFile(const char *filename);
////	NUMBER getSumLogPOPara(NCAT xndat, struct data **x_data); // generic per k/g-slice
////    NUMBER getLogLik(); // get log likelihood of the fitted model
//    NCAT getModelParamN(); // get log likelihood of the fitted model
////    NUMBER getNullSkillObs(NPAR m); // get log likelihood of the fitted model
//    // fitting (the only public method)
//    virtual void fit(); // return -LL for the model
//    // predicting
//    virtual void producePCorrect(NUMBER*** group_skill_map, NCAT* skills, NPAR n_skills, NUMBER* local_pred, struct context *ctx);
    virtual void updateValuesLocal(NUMBER*** group_skill_map, NCAT* skills, NPAR n_skills, NUMBER* local_pred, struct context *ctx);
    virtual void postPredictionSave(const char *filename);
//    static void predict(NUMBER* metrics, const char *filename,
//                        //NDAT N, NPAR* dat_obs, NCAT *dat_group, NCAT *dat_skill, NCAT *dat_skill_stacked, NCAT *dat_skill_rcount, NDAT *dat_skill_rix,
//                        struct task* task,
//                        HMMProblemEloSt **hmms, NPAR nhmms, NPAR *hmm_idx);
//	
//    void readModel(const char *filename, bool overwrite);
//    virtual void readModelBody(FILE *fid, struct task* task, NDAT *line_no, bool overwrite);
    NUMBER (HMMProblemEloSt::*sensitivity)(struct context*) = NULL; // pointer to the function that constructs sensitiviy Elo multiplier for a data row
    NUMBER sensitivity_K(struct context *ctx); // straight/direct formulation of sensitivity
    NUMBER sensitivity_U1b(struct context *ctx); // Uncertainty \frac{1}{1+b*n} U1b, as a variant of Uab \frac{a}{1+b*n}
    NUMBER sensitivity_U1bexp0_01(struct context *ctx); // Uncertainty \frac{1}{1+b*e^{0.01*n}} U1bexp0_01 as a variant of \frac{a}{1+b*e^{c*n}}
    NUMBER sensitivity_O2Z_IN_B(struct context *ctx); // sensitivity 1 to 0 in B steps in lenear fashion
protected:
//	//
//	// Givens
//	//
//    NCAT model_param_n; // number of model params
//    NUMBER* null_obs_ratios;
//    NUMBER null_skill_obs; // if null skills are present, what's the default obs to predict
//    NUMBER null_skill_obs_prob; // if null skills are present, what's the default obs probability to predict
//    NUMBER* param_skill; // skill parameters as concatenated vectorized forms (PI, A, B){nK}
//    NUMBER* gradient_skill; // skill parameters as concatenated vectorized forms (PI, A, B){nK}
    NUMBER* param_elo;
    NUMBER* gradient_elo;
    NPAR param_elo_n;
    NUMBER *elo_track_g; // elo g (student) rating
    NCAT *elo_count_g;   // per g (student) count
    NUMBER *elo_track_t; // elo g (student) rating for row, not student

//    NCAT skill1_n; // number of params per 1 skill slot
//	NUMBER* lb_param_skill; // lower boundary of parameters in vector form (PI, A, B){1}
//	NUMBER* ub_param_skill; // upper boundary of parameters in vector form (PI, A, B){1}
//	bool non01constraints; // whether there are lower or upper boundaries different from 0,1 respectively
//    NPAR *active_set_block; // which param block (whole skill param set – nK of them, whole student/group param set) is still actively fit
//    NPAR *unblocked_simplex; // which parameter simplex is not blocked from change (priors, A or B row, for students 1+2*nS of them)
//	struct task *task; // data and params
//	//
//	// Derived
//	//
//    NUMBER** alpha; // forward variables, 2D array of Nst rows nS values in each ((1..nS)){Nst}
//    NUMBER** po_param_gk; // p(O|param) likelihood of student/skill sequence, indexed by g and k
//    NUMBER* loglik_k; // p(O|param) log-likelihood of student/skill sequences, indexed by k -- it is logarithm scale, not probability like po_param_gk
//    NDAT* ndat_k; // number of datapoints, indexed by k
//    NUMBER** beta; // backward variables, 2D array of Nst rows nS values in each ((1..nS)){Nst}
//    NUMBER** gamma; // gamma variables, 2D array of Nst rows nS values in each ((1..nS)){Nst}
//    NUMBER*** xi; // xi variables, 3D array of Nst rows nS x nS in each
//    NUMBER* scale; // scale variable, 1D array of Nst rows 1 values in each
//    NDAT* backward_ix; // Nst sized array of pointing to previous student-skill repetition, -1 means none before
//    NDAT* forward_ix; // Nst sized array of pointing to next student-skill repetition, -1 means no more further
//    bool is_fwd_bwd_built; // flag for noting whether backward_ix & forward_ix are already built
//
    virtual void init(struct task* task); // non-fit specific initialization
    virtual void initAlphaEtAl(); // reset variables before making a computeAlphaAndPOParam run

//	NUMBER computeAlphaAndPOParam(NUMBER *metrics_res); // return loglikelihood, alternatively predictions, and metrics too
//	void computeBeta();
//	void computeXiGamma();
//    void FitNullSkill(); // get neg-log-lik, SSE, ncorrect for null skills
//
//    // predicting
//    virtual NUMBER computeGradients(); // returns loglik from the alphas
//    virtual void createGradientScaleSimplexAffordances(FitBitSt *fb); // set scaling vectors of offsets, length, and number of them in fb
//    virtual void /*NUMBER*/ Cycle(); // return -LL for the model
//    virtual void CycleBit(FitBitSt *fb); // for 1 skill or 1 group, all 1 skill for all data
//    virtual void doLinearStep(FitBitSt *fb, bool direction); // direction=true when, after computing direction as per CGD, we make a linear step
//    virtual void doLagrangeStep(FitBitSt *fb);
//    virtual void doConjugateLinearStep(FitBitSt *fb);
//    virtual void doBaumWelchStep(FitBitSt *fb);
//    virtual void doBarzilaiBorweinStep(FitBitSt *fb);
////    NUMBER doBarzsilaiBorweinStep(NCAT xndat, struct data** x_data, NUMBER *a_PI, NUMBER **a_A, NUMBER **a_B, NUMBER *a_PI_m1, NUMBER **a_A_m1, NUMBER **a_B_m1, NUMBER *a_gradPI_m1, NUMBER **a_gradA_m1, NUMBER **a_gradB_m1, NUMBER *a_gradPI, NUMBER **a_gradA, NUMBER **a_gradB, NUMBER *a_dirPI_m1, NUMBER **a_dirA_m1, NUMBER **a_dirB_m1);
////    NUMBER BaumWelch(); // return -LL for the model
//    void readNullObsRatio(FILE *fid, struct task* task, NDAT *line_no);
//	bool checkSkillConstraints(NUMBER* param_skill); // all constraints, inc row sums
private:
    // write model
	void toFileSkill(const char *filename);
//	void toFileGroup(const char *filename);
};

#endif
