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

#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <map>
#include <list>
#include "utilsSt.h"
#include "InputUtilSt.h"
#include "HMMProblemSt.h"
////#include "HMMProblemPiG.h"
//#include "HMMProblemPiGK.h"
//#include "HMMProblemPiGKww.h"
//#include "HMMProblemAGK.h"
////#include "HMMProblemAGKi.h"
//#include "HMMProblemPiAGK.h"
//#include "HMMProblemPiABGK.h"
////#include "HMMProblemKT.h"
//#include "HMMProblemSlicedAB.h"
//#include "HMMProblemSlicedA.h"
//#include "StripedArray.h"
//#include "HMMProblemComp.h"
//#include "HMMProblemEloK.h"
#include "HMMProblemEloSt.h"
//#include "SparseArray2D.h"
//#include <boost/numeric/ublas/matrix_sparse.hpp>//BOOST
//#include <boost/numeric/ublas/io.hpp>//BOOST
using namespace std;

struct task task;
void exit_with_help();

void parse_arguments_step1(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name, char *console_file_name); // things that do not need data file read
void parse_arguments_step2(int argc, char **argv, FILE *fid_console); // things that do need data file read, namely, number of observations

bool read_and_structure_data(const char *filename, FILE *fid_console);
NUMBER cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
NUMBER cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
NUMBER cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
//NUMBER cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR
//NUMBER cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR
//NUMBER cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR

static int max_line_length;
static char * line;
static char* readline(FILE *fid) {
	int length = 0;
	
	if(fgets(line,max_line_length,fid) == NULL)
		return NULL;
	
	while(strrchr(line,'\n') == NULL && strrchr(line,'\r') == NULL) // do take both line endings
	{
		max_line_length *= 2;
		line = (char *) realloc(line, (size_t)max_line_length);
		length = (int) strlen(line);
		if(fgets(line+length,max_line_length-length,fid) == NULL)
			break;
	}
	return line;
}

// temporary experimental: IRT-like for fitting pLo in liblinear
/*
void write_pLo_irt() {
    FILE *fid0 = fopen("uopx12_irt.txt","w");
    NPAR **group_skill_mask = init2D<NPAR>(task.nG, task.nK);
    NCAT g_k, g, k;
    NDAT t;
    data *dat;
    NPAR obs;
    for(g=0; g<task.nG; g++) {
        g_k = task.g_numk[g];
        for(k=0; k<g_k; k++) {
            dat = task.g_k_data[g][ k ];
            t = dat->ix[0];
            NCAT *ar;
            int n = 0;
            if(task.multiskill==0) {
                k = task.dat_skill[t];
                ar = &k;
                n = 1;
            } else {
//                ar = &task.dat_multiskill->get(t)[1];
//                n = task.dat_multiskill->get(t)[0];
                k = task.dat_skill_stacked[ task.dat_skill_rix[t] ];
                ar = &task.dat_skill_stacked[ task.dat_skill_rix[t] ];
                n = task.dat_skill_rcount[t];
                qsortNcat(ar, (NPAR)n);
            }
            obs = task.dat_obs[ dat->ix[0] ]; //->get( dat->ix[0] );
            NPAR count = 0; // 557687 >> 499117
            for(int l=0; l<n; l++) {
                count = (NPAR)(count + (group_skill_mask[g][ ar[l] ] == 1));
                if(count<n) {
                    fprintf(fid0,"%s %u:1", ((1-obs)==0)?"-1":"+1",dat->g+1);
                    
                    for(int l=0; l<n; l++) {
                        fprintf(fid0, " %u:1",ar[l]+task.nG+1);
                        group_skill_mask[g][ ar[l] ] = 1;
                    }
                    fprintf(fid0,"\n");
                }
            }
        }
    }
    fclose(fid0);
    free2D(group_skill_mask, task.nG);
}*/

int main (int argc, char ** argv) {
	clock_t tm_all = clock();//overall time //SEQ
//    double _tm_all = omp_get_wtime(); //PAR
    
	char input_file[1024]; // data
	char output_file[1024]; // model
    char colsole_file[1024]; // console copy
	char predict_file[1024]; // predictions
    
	set_task_defaults(&task);
    
    
    // parse parameters, step 1
	parse_arguments_step1(argc, argv, input_file, output_file, predict_file, colsole_file);

    FILE *fid_console = NULL;
    if(task.duplicate_console==1)
        fid_console = fopen(colsole_file,"w");
    
    if(!task.quiet) {
        printf("trainhmmst starting...\n");
        if(task.duplicate_console==1) fprintf(fid_console, "trainhmmst starting...\n");
    }

    clock_t tm_read = clock();//overall time //SEQ
//    double _tm_read = omp_get_wtime(); //PAR
    int read_ok = read_and_structure_data(input_file, fid_console);
    tm_read = (clock_t)(clock()-tm_read);//SEQ
    
//    _tm_read = omp_get_wtime()-_tm_read;//PAR

    if( ! read_ok )
        return 0;
    
    // once we know nO (number of observations) parse parameters, step 2
    parse_arguments_step2(argc, argv, fid_console);
    
//    write_pLo_irt();
    
    
//    //
//    // read item mean % correct
//    //
//    FILE *fid = fopen("a89_kts_train01_voc_i.txt","r");
//	max_line_length = 1024;
//	char *col;
//    string item;
//    map<string,NCAT>::iterator it;
//    task.item_complexity = Calloc(NUMBER, task.map_step_bwd->size());
//    line = (char *)malloc(max_line_length);// Malloc(char,max_line_length);
//    while( readline(fid)!=NULL) {
//		// Group
//		col = strtok(line,"\t\n\r");
//		item = string( col );
//		it = task.map_step_fwd->find(item);
//		if( it==task.map_step_fwd->end() ) { // not found
//            fprintf(stderr,"DID NOT FIND THE STEP!!\n");
//            return false;
//		}
//		else {
////            if( it->second > task.map_step_bwd->size()) {
////                int z = 0 ;
////            }
//            NUMBER v = atof( strtok(NULL,"\t\n\r") );
//            task.item_complexity[ it->second ] = v;
//		}
//    }
//    free(line);
//    fclose(fid);

    if(!task.quiet) {
        printf("input read, nO=%d, nG=%d, nK=%d, nI=%d, nZ=%d\n",task.nO, task.nG, task.nK, task.nI, task.nZ);
        if(task.duplicate_console==1) fprintf(fid_console, "input read, nO=%d, nG=%d, nK=%d, nI=%d, nZ=%d\n",task.nO, task.nG, task.nK, task.nI, task.nZ);
    }
    
    //    // write time
    //    if(task.time==1) {
    //        const char * fn = "a89_kts_times.txt";
    ////        const char * fn = "a89_uskts_times.txt";
    //        write_time_interval_data(&param, fn);
    //    }
    
    // erase blocking labels
//    zeroLabels(&task);

    clock_t tm_fit = 0; //SEQ
    clock_t tm_predict = 0; //SEQ
//    double _tm_fit;//PAR
//    double _tm_predict;//PAR
    
    if(task.cv_folds==0) { // not cross-validation
        // create problem
        HMMProblemSt *hmm = NULL;
        switch(task.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
//            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmm = new HMMProblemSt(&task);
                break;
                //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
                //                hmm = new HMMProblemPiG(&param);
                //                break;
//            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
//                hmm = new HMMProblemSlicedAB(&param);
//                break;
//            case STRUCTURE_SKAslc: // Conjugate Gradient Descent
//                hmm = new HMMProblemSlicedA(&param);
//                break;
//            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                hmm = new HMMProblemPiGK(&param);
//                break;
//            case STRUCTURE_PIgkww: // Gradient Descent, pLo=f(K,G), other by K
//                hmm = new HMMProblemPiGKww(&param);
//                break;
//            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmm = new HMMProblemPiAGK(&param);
//                break;
//            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
//                hmm = new HMMProblemAGK(&param);
//                break;
////                case STRUCTURE_Agki: // Gradient Descent, pT=f(K,G), other by K
////                    hmm = new HMMProblemAGKi(&param);
////                    break;
//            case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
//                hmm = new HMMProblemPiABGK(&param);
//                break;
//                //            case BKT_GD_T: // Gradient Descent with Transfer
//                //                hmm = new HMMProblemKT(&param);
//                //                break;
//            case STRUCTURE_COMP: // Gradient Descent, pT=f(K,G), other by K
//                hmm = new HMMProblemComp(&param);
//                break;
            case STRUCTURE_ELO: // Gradient Descent, pT=f(K,G), other by K
                hmm = new HMMProblemEloSt(&task);
                break;
            default:
                fprintf(stderr,"Solver specified (%d) is not supported!\n",task.structure);
                exit(1);
                break;
        }
        tm_fit = clock(); //SEQ
//        _tm_fit = omp_get_wtime(); //PAR
        hmm->fit();
        tm_fit = clock()-tm_fit;//SEQ
//        _tm_fit = omp_get_wtime()-_tm_fit;//PAR
        
        // write model
        hmm->toFile(output_file);
        
        if(task.metrics>0 || task.predictions>0) {
            NUMBER* metrics = Calloc(NUMBER, (size_t)8); // LL, LLnonull, RMSE, RMSEnonull, Acc, Acc_nonull, AIC, BIC;
            // takes care of predictions and metrics, writes predictions if task.predictions==1
            
            // temporary
//            if(task.per_kc_rmse_acc) {
//                task.kc_counts = Calloc(NDAT, (size_t)task.nK);
//                task.kc_rmse = Calloc(NUMBER, (size_t)task.nK);
//                task.kc_acc  = Calloc(NUMBER, (size_t)task.nK);
//            }

            // NUMBER l1 = hmm->getSumLogPOPara(task.nSeq, task.k_data);
//            printf("hmm-style ll_no_null %15.7f\n",l1);
            
            tm_predict = clock(); //SEQ
//            _tm_predict = omp_get_wtime(); //PAR
//            if(task.structure==STRUCTURE_SKAslc) {
//                ((HMMProblemSlicedA *)hmm)->predict(metrics, predict_file, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix);
//            } else if(task.structure==STRUCTURE_SKABslc) {
//                ((HMMProblemSlicedAB *)hmm)->predict(metrics, predict_file, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix);
//            } else {
                HMMProblemSt::predict(metrics, predict_file,
                                      //task.N, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix,
                                      &task,
                                      &hmm, 1, NULL);
            
//            }
            tm_predict = clock()-tm_predict;//SEQ
//            _tm_predict = omp_get_wtime()-_tm_predict;//PAR
            
            if( task.metrics>0 /*&& !task.quiet*/) {
                printf("trained model LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
                       metrics[0], metrics[1], // ll's
                       2*hmm->getModelParamN() + 2*metrics[0], // AIC
                       hmm->getModelParamN()*safelog(task.N) + 2*metrics[0], //BIC
                       metrics[2], metrics[3], // rmse's
                       metrics[4], metrics[5]); // acc's
                if(task.duplicate_console==1) fprintf(fid_console, "trained model LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
                       metrics[0], metrics[1], // ll's
                       2*hmm->getModelParamN()+ 2*metrics[0], // AIC
                       hmm->getModelParamN()*safelog(task.N) + 2*metrics[0],
                       metrics[2], metrics[3], // rmse's
                       metrics[4], metrics[5]); // acc's

            }
            free(metrics);

            // temporary
//            if(task.per_kc_rmse_acc) {
//                for(NCAT i=0; i<task.nK; i++) {
//                    printf("KC %4u RMSE=%8.6f Acc=%8.6f\n",i,task.kc_rmse[i],task.kc_acc[i]);
//                    if(task.duplicate_console==1) fprintf(fid_console, "KC %4u RMSE=%8.6f Acc=%8.6f\n",i,task.kc_rmse[i],task.kc_acc[i]);
//                }
//                free(task.kc_counts);
//                free(task.kc_rmse);
//                free(task.kc_acc);
//            }
        
        } // if predict or metrics
        
        delete hmm;
    } else { // cross-validation
        NUMBER* metrics = Calloc(NUMBER, (size_t)7); // AIC, BIC, RMSE, RMSE no null
		NUMBER n_par = 0;
        switch (task.cv_strat) {
            case CV_GROUP:
                n_par = cross_validate(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                n_par = cross_validate(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            case CV_ITEM:
                n_par = cross_validate_item(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                n_par = cross_validate_item(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            case CV_NSTR:
                n_par = cross_validate_nstrat(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                n_par = cross_validate_nstrat(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            default:
                
                break;
        }

		printf("%d-fold cross-validation: LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
			   task.cv_folds,
			   metrics[0], metrics[1], // ll's
			   2*n_par + 2*metrics[0], n_par*safelog(task.N) + 2*metrics[0],
			   metrics[2], metrics[3], // rmse's
			   metrics[4], metrics[5]); // acc's
		if(task.duplicate_console==1) fprintf(fid_console, "%d-fold cross-validation: LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
											   task.cv_folds,
											   metrics[0], metrics[1], // ll's
											   2*n_par + 2*metrics[0], n_par*safelog(task.N) + 2*metrics[0],
											   metrics[2], metrics[3], // rmse's
											   metrics[4], metrics[5]); // acc's
		
		
        free(metrics);
    }
	// free data
	destroy_input_data(&task);
	
//	if(task.quiet == 0) {
        printf("timing: overall %f seconds, read %f, fit %f, predict %f\n",(NUMBER)((clock()-tm_all)/CLOCKS_PER_SEC), (NUMBER)tm_read/CLOCKS_PER_SEC,  (NUMBER)tm_fit/CLOCKS_PER_SEC,  (NUMBER)tm_predict/CLOCKS_PER_SEC);//SEQ
        if(task.duplicate_console==1) fprintf(fid_console, "timing: overall %f seconds, read %f, fit %f, predict %f\n",(NUMBER)((clock()-tm_all)/CLOCKS_PER_SEC), (NUMBER)tm_read/CLOCKS_PER_SEC,  (NUMBER)tm_fit/CLOCKS_PER_SEC,  (NUMBER)tm_predict/CLOCKS_PER_SEC);//SEQ
//        if(task.duplicate_console==1) fprintf(fid_console, "timing: overall %lf sec, read %lf sec, fit %lf sec, predict %lf sec\n",omp_get_wtime()-_tm_all, _tm_read, _tm_fit, _tm_predict);//PAR
//        printf("timing: overall %lf sec, read %lf sec, fit %lf sec, predict %lf sec\n",omp_get_wtime()-_tm_all, _tm_read, _tm_fit, _tm_predict);//PAR
//    }
    
    if(task.duplicate_console==1)
        fclose(fid_console);
    return 0;
}

void exit_with_help() {
	printf(
		   "Usage: trainhmmst [options] input_file [[output_file] predicted_response_file]\n"
           "options:\n"
           "-s : structure.solver[.solver setting], structures: 1-by skill, 2-by user;\n"
           "     solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent;\n"
           "     Conjugate Gradient Descent has 3 settings: 1-Polak-Ribiere,\n"
           "     2-Fletcher–Reeves, 3-Hestenes-Stiefel, 4-Dai-Yuan.\n"
           "     For example '-s 1.3.1' would be by skill structure (classical) with\n"
           "     Conjugate Gradient Descent and Hestenes-Stiefel formula, '-s 2.1' would be\n"
           "     by student structure fit using Baum-Welch method.\n"
           "-e : tolerance of termination criterion (0.01 for parameter change default);\n"
           "     could be compuconvergeted by the change in log-likelihood per datapoint, e.g.\n"
           "     '-e 0.00001,l'.\n"
           "-i : maximum iterations (200 by default)\n"
           "-q : quiet mode, without output, 0-no (default), or 1-yes\n"
           "-n : number of hidden states, should be 2 or more (default 2)\n"
           "-0 : initial parameters comma-separated for priors, transition, and emission\n"
           "     probabilities skipping the last value from each vector (matrix row) since\n"
           "     they sum up to 1; default 0.5,1.0,0.4,0.8,0.2\n"
           "-l : lower boundaries for parameters, comma-separated for priors, transition,\n"
           "     and emission probabilities (without skips); default 0,0,1,0,0,0,0,0,0,0\n"
           "-u : upper boundaries for params, comma-separated for priors, transition,\n"
           "     and emission probabilities (without skips); default 0,0,1,0,0,0,0,0,0,0\n"
           "-c : specification of the C weight and cetroids for L2 penalty, empty (default).\n"
           "     For standard BKT - 4 comma-separated numbers: C weight of the penalty and\n"
           "     centroids, for PI, A, and B matrices respectively. If used for iBKT with\n"
           "     student effects, 8 values will be used with 4 additional values for student\n"
           "     effect matrices. For example, '-c 1.0,0.5,0.5,0.0'.\n"
           "-f : fit as one skill, 0-no (default), 1 - fit as one skill and use params as\n"
           "     starting point for multi-skill, 2 - force one skill\n"
           "-m : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To \n"
           "     specify observation for which metrics to be reported, list it after ','.\n"
           "     For example '-m 0', '-m 1' (by default, observation 1 is assumed), '-m 1,2'\n"
           "     (compute metrics for observation 2). Incompatible with-v option.\n"
           "-v : cross-validation folds, stratification, and target state to validate\n"
           "     against, default 0 (no cross-validation),\n"
           "     examples '-v 5,i,2' - 5 fold, item-stratified c.-v., predict state 2,\n"
           "     '-v 10' - 10-fold subject-stratified c.-v. predict state 1 by default,\n"
           "     alternatively '-v 10,g,1', and finally '-v 5,n,2,' - 5-fold unstratified\n"
           "     c.-v. predicting state 1.\n"
           "-p : report model predictions on the train set 0-no (default), 1-yes; 2-yes,\n"
           "     plus output state probability; works with -v and -m parameters.\n"
           "-U : controls how update to the probability distribution of the states is\n"
           "     updated. Takes the following format '-U r|g[,t|g]', where first\n"
           "     character controls how prediction treats known observations, second -- how\n"
           "     prediction treats unknown observations, and third -- whether to output\n"
           "     probabilities of priors. Dealing with known observations 'r' - reveal\n"
           "     actual observations for the update of state probability distribution (makes\n"
           "     sense for modeling how an actual system would work), 'g' - 'guessing' the\n"
           "     observation based on the predicted outcomes (arg max) -- more appropriate\n"
           "     when comparing models (so that no information about observation is never\n"
           "     revealed). Dealing with unknown observations (marked as '.' -- dot): 't' --\n"
           "     use transition matrix only, 'g' -- 'guess' the observation.\n"
           "     Default (if ommitted) is '-U r,t'.\n"
           "     For examle, '-U g,g would require 'guessing' of what the observation was\n"
           "     using model parameters and the running value of the probabilities of state\n"
           "     distributions.\n"
           "-d : delimiter for multiple skills per observation; 0-single skill per\n"
           "     observation (default), otherwise -- delimiter character, e.g. '-d ~'.\n"
           "-b : treat input file as binary input file (specifications TBA).\n"
           "-B : block re-estimation of prior, transitions, or emissions parameters\n"
           "     respectively (defailt is '-B 0,0,0'), to block re-estimation of transition\n"
           "     probabilities specify '-B 0,1,0'.\n"
           "-P : use parallel processing, defaul - 0 (no parallel processing), 1 - fit\n"
           "     separate skills/students separately, 2 - fit separate sequences within\n"
           "     skill/student separately.\n"
           "-o : in addition to printing to console, print output to the file specified\n"
           "     default is empty.\n"
		   );
	exit(1);
}

void parse_arguments_step1(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name, char *console_file_name) {
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
    
    // at this time we do not know nO (the number of observations) yet
	int i;
    int n;
    char *ch, *ch2;
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break; // end of options stop parsing
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 'e':
				task.tol = atof( strtok(argv[i],",\t\n\r") );
                ch = strtok(NULL,",\t\n\r"); // could be NULL
                if(ch != NULL)
                    task.tol_mode = ch[0];
				if(task.tol<0) {
					fprintf(stderr,"ERROR! Fitting tolerance cannot be negative\n");
					exit_with_help();
				}
                if(task.tol>10) {
                    fprintf(stderr,"ERROR! Fitting tolerance cannot be >10\n");
                    exit_with_help();
                }
                if(task.tol_mode!='p' && task.tol_mode!='l') {
                    fprintf(stderr,"ERROR! Tolerance mode '%c' is not allowed\n",task.tol_mode);
                    exit_with_help();
                }
				break;
			case 't':
				task.sliced = (NPAR)atoi(argv[i]);
				if(task.sliced!=0 && task.sliced!=1) {
					fprintf(stderr,"ERROR! Time parameter should be either 0 (off) or 1(on)\n");
					exit_with_help();
				}
				break;
			case 'i':
				task.maxiter = atoi(argv[i]);
				if(task.maxiter<10) {
					fprintf(stderr,"ERROR! Maximum iterations should be at least 10\n");
					exit_with_help();
				}
				break;
			case 'q':
				task.quiet = (NPAR)atoi(argv[i]);
				if(task.quiet!=0 && task.quiet!=1) {
					fprintf(stderr,"ERROR! Quiet param should be 0 or 1\n");
					exit_with_help();
				}
				break;
			case 'n':
				task.nS = (NPAR)atoi(argv[i]);
				if(task.nS<2) {
					fprintf(stderr,"ERROR! Number of hidden states should be at least 2\n");
					exit_with_help();
				}
                if(task.nS != 2) {
                    task.stat_specd_gt2 = true;
                }
				break;
			case 'S':
				task.is_scaled = (NPAR)atoi(argv[i]);
				if(task.is_scaled < 0 || task.is_scaled > 1) {
					fprintf(stderr,"ERROR! Scaling flag should be either 0 (off) or 1 (in)\n");
					exit_with_help();
				}
				break;
			case 's':
				task.structure = (NPAR)atoi( strtok(argv[i],".\t\n\r") );
                ch = strtok(NULL,".\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL)
                    task.solver = (NPAR)atoi(ch);
                ch = strtok(NULL,"\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL)
                    task.solver_setting = (NPAR)atoi(ch);
                if( task.structure != STRUCTURE_SKILL && // task.structure != STRUCTURE_GROUP &&
                   task.structure != STRUCTURE_PIg   && task.structure != STRUCTURE_PIgk  &&
                   task.structure != STRUCTURE_PIAgk && task.structure != STRUCTURE_Agk &&
                   task.structure != STRUCTURE_PIABgk && task.structure != STRUCTURE_Agki &&
                   task.structure != STRUCTURE_PIgkww && task.structure != STRUCTURE_SKABslc &&
                   task.structure != STRUCTURE_SKAslc && task.structure != STRUCTURE_COMP &&
                   task.structure != STRUCTURE_ELO ) {
                    fprintf(stderr, "Model Structure specified (%d) is out of range of allowed values\n",task.structure);
					exit_with_help();
                }
                if( task.solver != METHOD_BW  && task.solver != METHOD_GD &&
                   task.solver != METHOD_CGD && task.solver != METHOD_GDL &&
                   task.solver != METHOD_GBB) {
                    fprintf(stderr, "Method specified (%d) is out of range of allowed values\n",task.solver);
					exit_with_help();
                }
//                if( (task.structure == STRUCTURE_SKABslc || task.structure == STRUCTURE_SKAslc) && task.solver == METHOD_BW) {
//                    fprintf(stderr, "Method specified (%d) is not defined for this structure (%d) \n",task.solver,task.structure);
//                    exit_with_help();
//                }
                if( task.structure == STRUCTURE_COMP && task.solver == METHOD_BW) {
                    fprintf(stderr, "Method specified (%d) is not defined for this structure (%d) \n",task.solver,task.structure);
                    exit_with_help();
                }
//                if( task.structure == STRUCTURE_ELO && !(task.solver == METHOD_GD || task.solver == METHOD_CGD) ) {
//                    fprintf(stderr, "Method specified (%d) is not defined for this structure (%d) \n",task.solver,task.structure);
//                    exit_with_help();
//                }
                if( task.solver == METHOD_BW && ( task.solver != STRUCTURE_SKILL /*&& task.solver != STRUCTURE_GROUP*/ ) ) {
                    fprintf(stderr, "Baum-Welch solver does not support model structure specified (%d)\n",task.solver);
					exit_with_help();
                }
                if( task.solver == METHOD_CGD  &&
                   ( task.solver_setting != 1 && task.solver_setting != 2 &&
                    task.solver_setting != 3 && task.solver_setting !=4 )
                   ) {
                    fprintf(stderr, "Conjugate Gradient Descent setting specified (%d) is out of range of allowed values\n",task.solver_setting);
					exit_with_help();
                }
				break;
            case 'f':
                task.single_skill = (NPAR)atoi(argv[i]);
                break;
			case 'm':
                task.metrics = atoi( strtok(argv[i],",\t\n\r"));
                ch = strtok(NULL, "\t\n\r");
                if(ch!=NULL) {
                    task.metrics_target_obs = (NPAR)(atoi(ch)-1);
                }
				if(task.metrics<0 || task.metrics>1) {
					fprintf(stderr,"value for -m should be either 0 or 1.\n");
					exit_with_help();
				}
				if(task.metrics_target_obs<0) {// || task.metrics_target_obs>(task.nO-1)) {
					fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",task.metrics_target_obs+1);
					exit_with_help();
				}
                break;
			case 'b':
                task.binaryinput = atoi( strtok(argv[i],"\t\n\r"));
                break;
			case 'v':
				task.cv_folds   = (NPAR)atoi( strtok(argv[i],",\t\n\r"));
                ch2 = strtok(NULL, ",\t\n\r");
                if(ch2!=NULL)
                    task.cv_strat = ch2[0];
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    task.cv_target_obs = (NPAR)(atoi(ch)-1);
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    strcpy(task.cv_folds_file, ch);
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    task.cv_inout_flag = ch[0];
                
				if(task.cv_folds<2) {
					fprintf(stderr,"number of cross-validation folds should be at least 2\n");

					exit_with_help();
				}
				if(task.cv_folds>10) {
					fprintf(stderr,"please keep number of cross-validation folds less than or equal to 10\n");
					exit_with_help();
				}
                if(task.cv_strat != CV_GROUP && task.cv_strat != CV_ITEM && task.cv_strat != CV_NSTR){
					fprintf(stderr,"cross-validation stratification parameter '%c' is illegal\n",task.cv_strat);
					exit_with_help();
                }
				if(task.cv_target_obs<0) {// || task.cv_target_obs>(task.nO-1)) {
					fprintf(stderr,"target observation to be cross-validated against cannot be '%d'\n",task.cv_target_obs+1);
					exit_with_help();
				}
                if( task.cv_inout_flag!='i' && task.cv_inout_flag!='o') {
					fprintf(stderr,"cross-validation folds input/output flag should be ither 'o' (out) or 'i' (in), while it is '%c'\n",task.cv_inout_flag);
					exit_with_help();
                }
                
				break;
            case  'p':
                task.predictions = atoi(argv[i]);
                if(task.predictions<0 || task.predictions>3) {
                    fprintf(stderr,"a flag of whether to report predictions for training data (-p) should be 0, 1, 2, or 3\n");
                    exit_with_help();
                }
                break;
            case  'U':
                task.update_known = *strtok(argv[i],",\t\n\r");
                ch = strtok(NULL, ",\t\n\r");
                task.update_unknown = ch[0];
                
                if( (task.update_known!='r' && task.update_known!='g') ||
                   (task.update_unknown!='t' && task.update_unknown!='g') ) {
                    fprintf(stderr,"specification of how probabilities of states should be updated (-U) is incorrect, it sould be r|g[,t|g] \n");
                    exit_with_help();
                }
                break;
            case  'd':
				task.multiskill = argv[i][0]; // just grab first character (later, maybe several)
                break;
            case  'P':
				n = atoi(argv[i]);
                if(n!=0 && n!=1 && n!=2) {
					fprintf(stderr,"parallel processing flag (-P) should be 0 or 1\n");
					exit_with_help();
                }
                task.parallel = (NPAR)n;
                break;
            case 'r': // coordinate descend parameters
                // if two first_iteration_qualify,iterations_to_qualify
                // if one iterations_to_qualify (first_iteration_qualify==0)
                n = atoi( strtok(argv[i],",\t\n\r") );
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch==NULL) {// one parameter
                    task.first_iteration_qualify = 0;
                    task.iterations_to_qualify   = (NPAR)n;
                } else {
                    task.first_iteration_qualify = (NPAR)n;
                    task.iterations_to_qualify   = (NPAR)atoi(ch);
                }
                break;
            case 'R': // coordinate descend parameter for hard iterations limit of skill/group/other fitting
                n = atoi(argv[i]);
                if(n<=0) {
                    fprintf(stderr,"skill/group gradient descent iteration limit (-R) should be an integer >0\n");
                    exit_with_help();
                }
                if(n<=0) {
                    fprintf(stderr,"skill/group gradient descent iteration limit (-R) should be an integer >0\n");
                    exit_with_help();
                }
                task.iterations_limit = (NPAR)n;
                break;
            /*
            case 'c': {
                    // this version is with single center of gravity per Pi, A, and B
                    StripedArray<NUMBER> * tmp_array = new StripedArray<NUMBER>();
                    ch = strtok(argv[i],",\t\n\r");
                    while( ch != NULL ) {
                        tmp_array->add( atof(ch) );
                        ch = strtok(NULL,",\t\n\r");
                    }
                    if( (tmp_array->getSize() % 4) != 0 ) {
                        fprintf(stderr,"The number of regularization parameters should be a multiple of 4 and it is %d\n",tmp_array->getSize());
                        exit_with_help();
                    }
                    task.Cslices = (NPAR) tmp_array->getSize() / 4;
                    task.Cw = Calloc(NUMBER, (size_t)task.Cslices);
                    task.Ccenters = Calloc(NUMBER, (size_t)(task.Cslices * 3) );
                    int c1 = 0, c2 = 0, i = 0;
                    for(int l=0; l<(int)tmp_array->getSize() / 4; l++) {
                        task.Cw[c1++] = tmp_array->get((NDAT)i++);
                        for(int j=0; j<3; j++)
                            task.Ccenters[c2++] = tmp_array->get((NDAT)i++);
                    }
                    delete tmp_array;
                }
                break;
                */
            case 'c': // just to keep it a valid option
                break;
            case 'o':
                task.duplicate_console = 1;
                strcpy(console_file_name,argv[i]);
                break;
            case '0':
                task.init_reset = true;
                break;
            case 'l': // just to keep it a valid option
                break;
            case 'u': // just to keep it a valid option
                break;
            case 'B': // just to keep it a valid option
                break;
            case 'k': // just to keep it a valid option
                break;
			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
    // post-process checks
    // -v and -m collision
    if(task.cv_folds>0 && task.metrics>0) { // correct for 0-start coding
        fprintf(stderr,"values for -v and -m cannot be both non-zeros\n");
        exit_with_help();
    }
    // scaling
    if(task.is_scaled == 1 && task.solver != METHOD_BW) {
        task.is_scaled = 0;
        printf("Scaling can only be enabled for Baum-Welch method. Setting it to off\n");
    }
    // specifying >2 states via -n and mandatory specification of -0 (initial parameters)
    if(task.nS > 2 && !task.init_reset) {
        fprintf(stderr,"when >2 latent states specified via '-n', initial values of parameters have to be explicitly set via '-0'!\n");
        exit_with_help();
    }
    // STRUCTURE_SKABslc solver and -t 1 should be set together
    if( (task.sliced==1) != (task.structure == STRUCTURE_SKABslc || task.structure == STRUCTURE_SKAslc) ) {
        fprintf(stderr,"Error! sliced parameter ('-t 1') and STRUCTURE_SKABslc or STRUCTURE_SKAslc structure should be either both set on or off.\n");
        exit_with_help();
    }
    
	// next argument should be input file name
	if(i>=argc) // if not
		exit_with_help(); // leave
	
	strcpy(input_file_name, argv[i++]); // copy and advance
	
	if(i>=argc) { // no output file name specified
		strcpy(output_file_name,"output.hmm");
		strcpy(predict_file_name,"predict_hmm.txt"); // the predict file too
	}
	else {
		strcpy(output_file_name,argv[i++]); // copy and advance
		if(i>=argc) // no prediction file name specified
			strcpy(predict_file_name,"predict_hmm.txt"); // the predict file too
		else
			strcpy(predict_file_name,argv[i]);
	}
}

void parse_arguments_step2(int argc, char **argv, FILE *fid_console) {
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
    
    // at this time we do know nO (the number of observations)
	int i;
    int n, expected_n;
    char *ch;
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break; // end of options stop parsing
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case '0': // init_params
				int len;
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++) {
					n += (argv[i][j]==',')?1:0;
                    if( (argv[i][j] >= 'a' && argv[i][j] <= 'z') || (argv[i][j] >= 'A' && argv[i][j] <= 'Z') ) {
                        strcpy(task.initfile, argv[i]);
                        break;
                    }
                }
                // init parameters
                if(task.init_param_values!=NULL) free(task.init_param_values);
                task.init_param_values = Calloc(NUMBER, (size_t)n);
                task.init_param_values_n = (NPAR)n;
                // read params and write to params
                task.init_param_values[0] = atof( strtok(argv[i],",\t\n\r") );
                for(int j=1; j<n; j++) {
                    task.init_param_values[j] = atof( strtok(NULL,",\t\n\r") );
                }
				break;
			case 'l': // lower boundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (NPAR)((argv[i][j]==',')?1:0);
				// init params
				if(task.param_values_lb!=NULL) free(task.param_values_lb);
				task.param_values_lb = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				task.param_values_lb[0] = atof( strtok(argv[i],",\t\n\r") );
                for(int j=1; j<n; j++) {
					task.param_values_lb[j] = atof( strtok(NULL,",\t\n\r") );
//                    if(task.param_lo[j] >0) {
//                        int a = 0;
//                    }
                }
                // check if the number of lower boundaries corresponds
                expected_n = task.nS +
                    ((task.sliced==1 && (task.structure==STRUCTURE_SKAslc || task.structure==STRUCTURE_SKABslc) )?task.nZ:1)*task.nS*task.nS +
                    ((task.sliced==1 && task.structure==STRUCTURE_SKABslc )?task.nZ:1)*task.nS*task.nO;
                if( n != expected_n ) {
//                    fprintf(stderr,"Structure=%d, nS=%d, nO=%d, nZ=%d, The expected n=%d specified=%d\n",task.structure,task.nS,task.nO,task.nZ,expected_n,n);
                    fprintf(stderr,"The expected number of lower parameter boundaries is %d and it is %d\n",expected_n,n);
                    exit_with_help();
                }
                task.lb_specd = true;
				break;
			case 'u': // upper boundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (argv[i][j]==',')?1:0;
				// init params
				if(task.param_values_ub!=NULL) free(task.param_values_ub);
				task.param_values_ub = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				task.param_values_ub[0] = atof( strtok(argv[i],",\t\n\r") );
                for(int j=1; j<n; j++) {
					task.param_values_ub[j] = atof( strtok(NULL,",\t\n\r") );
//                    if(task.param_hi[j] < 1) {
//                        int a = 0;
//                    }
                }
                // check if the number of upper boundaries corresponds
                expected_n = task.nS +
                    ((task.sliced==1 && (task.structure==STRUCTURE_SKAslc || task.structure==STRUCTURE_SKABslc) )?task.nZ:1)*task.nS*task.nS +
                    ((task.sliced==1 && task.structure==STRUCTURE_SKABslc )?task.nZ:1)*task.nS*task.nO;
                if( n != expected_n ) {
                    fprintf(stderr,"The expected number of upper parameter boundaries is %d and it is %d\n",expected_n,n);
                    exit_with_help();
                }
                task.ub_specd = true;
				break;
			case 'B': // block fitting
                // first
				task.block_fitting[0] = (NPAR)atoi( strtok(argv[i],",\t\n\r") );
                if(task.block_fitting[0]!=0 && task.block_fitting[0]!=1) {
                    fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                    if(task.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                    exit_with_help();
                }
                // second
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL) {
                    task.block_fitting[1] = (NPAR)atoi(ch);
                    if(task.block_fitting[1]!=0 && task.block_fitting[1]!=1) {
                        fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        if(task.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        exit_with_help();
                    }
                }
                else {
                    fprintf(stderr,"There should be 3 blockig the fitting flags specified.\n");
                    if(task.duplicate_console==1) fprintf(fid_console,"There should be 3 blockig the fitting flags specified.\n");
                    exit_with_help();
                }
                // third
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL) {
                    task.block_fitting[2] = (NPAR)atoi(ch);
                    if(task.block_fitting[2]!=0 && task.block_fitting[2]!=1) {
                        fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        if(task.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        exit_with_help();
                    }
                }
                else {
                    fprintf(stderr,"There should be 3 blockig the fitting flags specified.\n");
                    if(task.duplicate_console==1) fprintf(fid_console,"There should be 3 blockig the fitting flags specified.\n");
                    exit_with_help();
                }
				break;
			case 'c': {
                    // this version with per-parameter center of gravity in Pi, A, and B
                    NPAR nS = task.nS;
                    NPAR nO = task.nO;
                    int nPerSlice = 1 + nS + nS*nS + nS*nO; // parameters per slice
                    int nCenters = nS + nS*nS + nS*nO; // centers per slice
                    StripedArray<NUMBER> * tmp_array = new StripedArray<NUMBER>();
                    ch = strtok(argv[i],",\t\n\r");
                    while( ch != NULL ) {
                        tmp_array->add( atof(ch) );
                        ch = strtok(NULL,",\t\n\r");
                    }
                    if( (tmp_array->getSize() % nPerSlice) != 0 ) {
                        fprintf(stderr,"The number of regularization parameters should be a multiple of %d and it is %d\n",nPerSlice,tmp_array->getSize());
                        exit_with_help();
                    }
                    task.Cslices = (NPAR) ( tmp_array->getSize() / nPerSlice );
                    task.Cw = Calloc(NUMBER, (size_t)task.Cslices);
                    task.Ccenters = Calloc(NUMBER, (size_t)(task.Cslices * nCenters) );
                    int c1 = 0/*counter for weights*/, c2 = 0/*counter for centers*/, i = 0/*counter for tmp_array*/;
                    for(int l=0; l<(int)tmp_array->getSize() / nPerSlice; l++) {
                        task.Cw[c1++] = tmp_array->get((NDAT)i++);
                        for(int j=0; j<nCenters; j++)
                            task.Ccenters[c2++] = tmp_array->get((NDAT)i++);
                    }
                    delete tmp_array;
                }
                break;
            case 'k':
                // ELO parameters specified as a comma-delimited list of values, first is an integer, the rest are real
                len = (int)strlen( argv[i] );
                // count delimiters
                n = 1; // start with 1
                for(int j=0;j<len;j++)
                    n += (NPAR)((argv[i][j]==',')?1:0);
                if( n < 3 ) {
                    fprintf(stderr,"To use Elo, one needs at lest three parameters; %d were specified\n",n);
                    exit_with_help();
                }
                // init params
                task.elo_param_values_n = (NPAR)(n-2); // first parameter is Elo type
                if(task.elo_param_values!=NULL) free(task.elo_param_values);
                task.elo_param_values = Calloc(NUMBER, (size_t)(n-2));
                // read params and write to params
                task.elo_type = (NPAR)atoi( strtok(argv[i],",\t\n\r") );
                task.elo_scope = (NCAT)atof( strtok(NULL,",\t\n\r") );
                for(int j=2; j<n; j++) {
                    task.elo_param_values[j-2] = atof( strtok(NULL,",\t\n\r") );
                }
                break;
                
        } // end switch
    }// end for
    
    // post parse actions
    
    // post-argument checks
    if( task.cv_target_obs>(task.nO-1)) {
        fprintf(stderr,"target observation to be cross-validated against cannot be '%d'\n",task.cv_target_obs+1);
        if(task.duplicate_console==1) fprintf(fid_console,"target observation to be cross-validated against cannot be '%d'\n",task.cv_target_obs+1);
        exit_with_help();
    }
    if(task.metrics_target_obs>(task.nO-1)) {
        fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",task.metrics_target_obs+1);
        if(task.duplicate_console==1) fprintf(fid_console,"target observation to compute metrics against cannot be '%d'\n",task.metrics_target_obs+1);
        exit_with_help();
    }
    
    if( !( (task.elo_param_values_n==0 && task.elo_type==0 && task.structure != STRUCTURE_ELO) || // not an Elo
           (task.elo_param_values_n==1 && task.elo_type==1 && task.structure == STRUCTURE_ELO) || // Simple Elo
           (task.elo_param_values_n==1 && task.elo_type==2 && task.structure == STRUCTURE_ELO) || // relation Uncertainty \frac{1}{1+b*n} Elo
           (task.elo_param_values_n==1 && task.elo_type==3 && task.structure == STRUCTURE_ELO) || // relation Uncertainty \frac{1}{1+b*e^{0.01*n}} Elo
           (task.elo_param_values_n==1 && task.elo_type==4 && task.structure == STRUCTURE_ELO)    // linear trend from 1 to 0 in b steps
          ) ) {
        fprintf(stderr,"Elo type specified (%d) does not match the number of Elo parameters (%d)\n",task.elo_type,task.elo_param_values_n);
        exit_with_help();
    }
    
    if( task.structure == STRUCTURE_ELO && task.nS != 2) {
        fprintf(stderr,"Elo solvers only work with number of states equal to 2, while %d were specified!",task.nS);
        exit_with_help();
    }
    
    // upper and lower bounds
    if( (task.sliced==1) && (task.structure == STRUCTURE_SKABslc || task.structure == STRUCTURE_SKAslc) &&
       !(task.lb_specd && task.ub_specd)) {
        fprintf(stderr,"Error! when parameter slicing is enabled, lower and upper boundaries for all parameters have to be set explicitly.\n");
        exit_with_help();
    }
    if(!task.lb_specd && (task.nS!=2 || task.nO!=2) ) { // if not specified, and it's not 2-state 2-obs case, set to 0
        if(task.param_values_lb!=NULL) free(task.param_values_lb);
        task.param_values_lb = Calloc(NUMBER, (size_t)( task.nS*(1+task.nS+task.nO) ) );
    }
 
    if(!task.ub_specd && (task.nS!=2 || task.nO!=2) ) {  // if not specified, and it's not 2-state 2-obs case, set to 1
        if(task.param_values_ub!=NULL) free(task.param_values_ub);  // if not specified, set to 1
        task.param_values_ub = Calloc(NUMBER, (size_t)( task.nS*(1+task.nS+task.nO) ) );
        for(int j=0; j<( task.nS*(1+task.nS+task.nO) ); j++)
            task.param_values_ub[j] = (NUMBER)1.0;
    }

    
}

bool read_and_structure_data(const char *filename, FILE *fid_console) {
    bool readok = true;
    if(task.binaryinput==0) {
        readok = InputUtilSt::readTxt(filename, &task);
    } else {
        readok = InputUtilSt::readBin(filename, &task);
    }
    if(! readok ) {
        return false;
    }
    
	NDAT t = 0;
	NCAT k;
    
	// Pass A, find first and last instance of null skill
	for(t=0; t<task.N; t++) {
        if(task.multiskill==0)
            k = task.dat_skill[t];//[t];
        else
            k = task.dat_skill_stacked[ task.dat_skill_rix[t] ]; // first skill of multi-skill

		// null skill : just count
		if( k < 0 ) {
            task.N_null++;
            if(task.first_null_skill==-1) task.first_null_skill = t; // update once
            task.last_null_skill = t; // update always
			continue;
		}
	}

    return true;
}

NUMBER cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//NUMBER cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    NPAR f;
    NCAT g;
    NDAT t;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(task.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(task.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)task.nG);//[task.nG];
    NDAT *fold_counts = Calloc(NDAT, (size_t)task.cv_folds);
    
    srand ( (unsigned int)time(NULL) );
    // folds file
    if(task.cv_folds_file[0] > 0) { // file is specified
        if(task.cv_inout_flag=='i') {
            fid_folds = fopen(task.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(task.cv_inout_flag=='o') {
            fid_folds = fopen(task.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(g=0; g<task.nG; g++) {
        if( task.cv_folds_file[0]==0 || (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='o') ) { // output or default
            folds[g] = (NPAR)(rand() % task.cv_folds);
            if(task.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[g]); // write out
        } else if( (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(task.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[g] = (NPAR)(atoi(ch));
        }
    }
    if(task.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(task.cv_inout_flag=='i')
            free(line);
    }
	task.metrics_target_obs = task.cv_target_obs;

    // count number of items in each fold
    for(t=0; t<task.N; t++) fold_counts[ folds[task.dat_group[t]] ]++;

    
    // create and fit multiple problems
    HMMProblemSt* hmms[task.cv_folds];
    int q = task.quiet;
//    task.quiet = 1;
    for(f=0; f<task.cv_folds; f++) {
        switch(task.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
//            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSt(&task);
                break;
            case STRUCTURE_ELO: // Conjugate Gradient Descent
//            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemEloSt(&task);
                break;
//            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedAB(&param);
//                break;
//            case STRUCTURE_SKAslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedA(&param);
//                break;
////            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
////                hmm = new HMMProblemPiG(&param);
////                break;
//            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                hmms[f] = new HMMProblemPiGK(&param);
//                break;
//            case STRUCTURE_PIgkww: // Gradient Descent, pLo=f(K,G), other by K
//                hmms[f] = new HMMProblemPiGKww(&param);
//                break;
//            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiAGK(&param);
//                break;
//            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemAGK(&param);
//                break;
//            case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiABGK(&param);
//                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
//                break;
            default:
                fprintf(stderr,"Solver specified (%d) is not supported!\n",task.structure);
                exit(1);
                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, (size_t)fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<task.N; t++) {
            if( folds[ task.dat_group[t]/*->get(t)*/ ] == f ) {
                saved_obs[count_saved++] = task.dat_obs[t];
                task.dat_obs[t] = -1;//->set(t, -1);
            }
        }


        // now compute
        tm0 = clock(); //SEQ
//        _tm0 = omp_get_wtime(); //PAR
        hmms[f]->fit();
        *(tm_fit) += (clock_t)(clock()- tm0);//SEQ
//        *(tm_fit) += omp_get_wtime()-_tm0;//PAR
        
        // write model
        char fname[1024];
        sprintf(fname,"%s_%i",model_file_name,f);
        hmms[f]->toFile(fname);
        
        // UN-block respective data
        count_saved = 0;
        for(t=0; t<task.N; t++)
            if( folds[ task.dat_group[t]/*->get(t)*/ ] == f )
                task.dat_obs[t]=saved_obs[count_saved++];
        free(saved_obs);
        printf("fold %d is done\n",f+1);
        if(task.duplicate_console==1) fprintf(fid_console, "fold %d is done\n",f+1);
    }
    task.quiet = (NPAR)q;
    
    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
	
	// new prediction
	//		create a general fold-identifying array
	NPAR *dat_fold = Calloc(NPAR, task.N);
	for(NDAT t=0; t<task.N; t++) dat_fold[t] = folds[ task.dat_group[t] ];
	//		predict
	hmms[0]->predict(metrics, filename,
                     //task.N, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix,
                     &task,
                     hmms, task.cv_folds/*nhmms*/, dat_fold);
	
    printf("Cross-validation is done wrong, if A or B matrices are 'sliced'\n");
    
//    if(task.structure==STRUCTURE_SKAslc) {
//        ((HMMProblemSlicedA *)hmm)->predict(metrics, predict_file, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix);
//    } else if(task.structure==STRUCTURE_SKABslc) {
//        ((HMMProblemSlicedAB *)hmm)->predict(metrics, predict_file, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix);
//    } else {
//        HMMProblem::predict(metrics, predict_file, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix, &hmm, 1, NULL);
//    }
    
    free(dat_fold);
	
	*(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
	
	
    // delete problems
    NUMBER n_par = 0;
    for(f=0; f<task.cv_folds; f++) {
        n_par += hmms[f]->getModelParamN();
        delete hmms[f];
    }
    n_par /= task.cv_folds;
	
    free(folds);
    free(fold_counts);
	
	return (n_par);
}

NUMBER cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//NUMBER cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    NPAR f;
    NCAT I; // item
    NDAT t;
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(task.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(task.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)task.nI);
    NDAT *fold_counts = Calloc(NDAT, (size_t)task.cv_folds);
    srand ( (unsigned int)time(NULL) ); // randomize

    // folds file
    if(task.cv_folds_file[0] > 0) { // file is specified
        if(task.cv_inout_flag=='i') {
            fid_folds = fopen(task.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(task.cv_inout_flag=='o') {
            fid_folds = fopen(task.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(I=0; I<task.nI; I++) {
        if( task.cv_folds_file[0]==0 || (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='o') ) { // output or default
            folds[I] = (NPAR)(rand() % task.cv_folds); // produce folds
            if(task.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[I]); // write out
        } else if( (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(task.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[I] = (NPAR)(atoi(ch));
        }
    }
    if(task.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(task.cv_inout_flag=='i')
            free(line);
    }
    
    // count number of items in each fold
    for(t=0; t<task.N; t++) fold_counts[ folds[task.dat_item[t]/*->get(t)*/] ]++;
    // create and fit multiple problems
    HMMProblemSt* hmms[task.cv_folds];
    int q = task.quiet;
//    task.quiet = 1;
    for(f=0; f<task.cv_folds; f++) {
        switch(task.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
//            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSt(&task);
                break;
            case STRUCTURE_ELO: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemEloSt(&task);
                break;
//            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
//                hmm = new HMMProblemPiG(&param);
//                break;
//            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedAB(&param);
//                break;
//            case STRUCTURE_SKAslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedA(&param);
//                break;
//            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                hmms[f] = new HMMProblemPiGK(&param);
//                break;
//            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiAGK(&param);
//                break;
//            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemAGK(&param);
//                break;
//            case STRUCTURE_PIABgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiABGK(&param);
//                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
//                break;
            default:
                fprintf(stderr,"Solver specified (%d) is not supported!\n",task.structure);
                exit(1);
                break;
        }
        
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, (size_t)fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<task.N; t++) {
            if( folds[ task.dat_item[t]/*->get(t)*/ ] == f ) {
                saved_obs[count_saved++] = task.dat_obs[t];
                task.dat_obs[t] = -1;//->set(t, -1);
            }
        }
        // now compute
        tm0 = clock(); //SEQ
//        _tm0 = omp_get_wtime(); //PAR
        
        hmms[f]->fit();
        *(tm_fit) += (clock_t)(clock()- tm0);//SEQ
//        *(tm_fit) += omp_get_wtime()-_tm0;//PAR
        
        // write model
        char fname[1024];
        sprintf(fname,"%s_%i",model_file_name,f);
        hmms[f]->toFile(fname);
        
        // UN-block respective data
        count_saved = 0;
        for(t=0; t<task.N; t++)
            if( folds[ task.dat_item[t]/*->get(t)*/ ] == f )
                task.dat_obs[t]=saved_obs[count_saved++];//->set(t, saved_obs[count_saved++]);
        free(saved_obs);
        printf("fold %d is done\n",f+1);
        if(task.duplicate_console==1) fprintf(fid_console, "fold %d is done\n",f+1);
    }
    free(fold_counts);
    task.quiet = (NPAR)q;

    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
	
	// new prediction
	//		create a general fold-identifying array
	NPAR *dat_fold = Calloc(NPAR, task.N);
	for(NDAT t=0; t<task.N; t++) dat_fold[t] = folds[ task.dat_item[t] ];
	//		predict
	HMMProblemSt::predict(metrics, filename,
                          //task.N, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix,
                          &task,
                          hmms, task.cv_folds/*nhmms*/, dat_fold);

    printf("Cross-validation is done wrong, if A or B matrices are 'sliced'\n");
    
    free(dat_fold);
	
	*(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<task.cv_folds; f++) {
        n_par += hmms[f]->getModelParamN();
        delete hmms[f];
    }
    n_par /= task.cv_folds;
    free(folds);

	return n_par;
}

NUMBER cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//NUMBER cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    NPAR f;
    NCAT U; // unstratified
    NDAT t;
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(task.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(task.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)task.N);
    NDAT *fold_counts = Calloc(NDAT, (size_t)task.cv_folds);
	
    srand ( (unsigned int)time(NULL) ); // randomize
	
    // folds file
    if(task.cv_folds_file[0] > 0) { // file is specified
        if(task.cv_inout_flag=='i') {
            fid_folds = fopen(task.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(task.cv_inout_flag=='o') {
            fid_folds = fopen(task.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(task.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(U=0; U<task.N; U++) {
        if( task.cv_folds_file[0]==0 || (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='o') ) { // output or default
            folds[U] = (NPAR)(rand() % task.cv_folds); // produce folds
            if(task.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[U]); // write out
        } else if( (task.cv_folds_file[0] > 0 && task.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(task.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[U] = (NPAR)(atoi(ch));
        }
    }
    if(task.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(task.cv_inout_flag=='i')
            free(line);
    }
    
    
    // count number of items in each fold
    for(t=0; t<task.N; t++)  fold_counts[ folds[task.dat_item[t]/*->get(t)*/] ]++;
    // create and fit multiple problems
    HMMProblemSt* hmms[task.cv_folds];
    int q = task.quiet;
//    task.quiet = 1;
    for(f=0; f<task.cv_folds; f++) {
        switch(task.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
//            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSt(&task);
                break;
            case STRUCTURE_ELO: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemEloSt(&task);
                break;
//            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedAB(&param);
//                break;
//            case STRUCTURE_SKAslc: // Conjugate Gradient Descent
//                hmms[f] = new HMMProblemSlicedA(&param);
//                break;
////            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
////                hmm = new HMMProblemPiG(&param);
////                break;
//            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                hmms[f] = new HMMProblemPiGK(&param);
//                break;
//            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiAGK(&param);
//                break;
//            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemAGK(&param);
//                break;
//            case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiABGK(&param);
//                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
//                break;
            default:
                fprintf(stderr,"Solver specified (%d) is not supported!\n",task.structure);
                exit(1);
                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, (size_t)fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<task.N; t++) {
            if( folds[ task.dat_item[t]/*->get(t)*/ ] == f ) {
                saved_obs[count_saved++] = task.dat_obs[t];
                task.dat_obs[t]=-1;//->set(t, -1);
            }
        }
        // now compute
        tm0 = clock(); //SEQ
//        _tm0 = omp_get_wtime(); //PAR
        hmms[f]->fit();
        *(tm_fit) += (clock_t)(clock()- tm0);//SEQ
//        *(tm_fit) += omp_get_wtime()-_tm0;//PAR
        
        // write model
        char fname[1024];
        sprintf(fname,"%s_%i",model_file_name,f);
        hmms[f]->toFile(fname);
        
        // UN-block respective data
        count_saved = 0;
        for(t=0; t<task.N; t++)
            if( folds[ task.dat_item[t]/*->get(t)*/ ] == f )
                task.dat_obs[t]=saved_obs[count_saved++];//->set(t, saved_obs[count_saved++]);
        free(saved_obs);
//        if(q == 0) {
            printf("fold %d is done\n",f+1);
            if(task.duplicate_console==1) fprintf(fid_console, "fold %d is done\n",f+1);
//        }
    }
    free(fold_counts);
    task.quiet = (NPAR)q;
    
    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
	
	// new prediction
	//		predict
	HMMProblemSt::predict(metrics, filename,
                          //task.N, task.dat_obs, task.dat_group, task.dat_skill, task.dat_skill_stacked, task.dat_skill_rcount, task.dat_skill_rix,
                          &task,
                          hmms, task.cv_folds/*nhmms*/, folds);
	
    printf("Cross-validation is done wrong, if A or B matrices are 'sliced'\n");
    
	*(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<task.cv_folds; f++) {
        n_par += hmms[f]->getModelParamN();
        delete hmms[f];
    }
    n_par /= task.cv_folds;
    free(folds);
	
	return n_par;
}
