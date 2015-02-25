/*
 
 Copyright (c) 2012-2014, Michael (Mikhail) Yudelson
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

#include <string.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <map>
#include <list>
#include "utils.h"
#include "InputUtil.h"
#include "HMMProblem.h"
//#include "HMMProblemPiG.h"
#include "HMMProblemPiGK.h"
#include "HMMProblemPiGKww.h"
#include "HMMProblemAGK.h"
//#include "HMMProblemAGKi.h"
#include "HMMProblemPiAGK.h"
#include "HMMProblemPiABGK.h"
//#include "HMMProblemKT.h"
#include "HMMProblemSliced.h"
#include "StripedArray.h"
//#include "SparseArray2D.h"
//#include <boost/numeric/ublas/matrix_sparse.hpp>//BOOST
//#include <boost/numeric/ublas/io.hpp>//BOOST
using namespace std;

struct param param;
void exit_with_help();

void parse_arguments_step1(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name, char *console_file_name); // things that do not need data file read
void parse_arguments_step2(int argc, char **argv, FILE *fid_console); // things that do need data file read, namely, number of observations

bool read_and_structure_data(const char *filename, FILE *fid_console);
void cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
void cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
void cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console);//SEQ
//void cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR
//void cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR
//void cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console);//PAR

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
void write_pLo_irt() {
    FILE *fid0 = fopen("uopx12_irt.txt","w");
    NPAR **group_skill_mask = init2D<NPAR>(param.nG, param.nK);
    NCAT g_k, g, k;
    NDAT t;
    data *dat;
    NPAR obs;
    for(g=0; g<param.nG; g++) {
        g_k = param.g_numk[g];
        for(k=0; k<g_k; k++) {
            dat = param.g_k_data[g][ k ];
            t = dat->ix[0];
            NCAT *ar;
            int n = 0;
            if(param.multiskill==0) {
                k = param.dat_skill[t];
                ar = &k;
                n = 1;
            } else {
//                ar = &param.dat_multiskill->get(t)[1];
//                n = param.dat_multiskill->get(t)[0];
                k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
                ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
                n = param.dat_skill_rcount[t];
                qsortNcat(ar, (NPAR)n);
            }
            obs = param.dat_obs[ dat->ix[0] ]; //->get( dat->ix[0] );
            NPAR count = 0; // 557687 >> 499117
            for(int l=0; l<n; l++)
                count = (NPAR)(count + (group_skill_mask[g][ ar[l] ] == 1));
                if(count<n) {
                    fprintf(fid0,"%s %u:1", ((1-obs)==0)?"-1":"+1",dat->g+1);
                    
                    for(int l=0; l<n; l++) {
                        fprintf(fid0, " %u:1",ar[l]+param.nG+1);
                        group_skill_mask[g][ ar[l] ] = 1;
                    }
                    fprintf(fid0,"\n");
                }
        }
    }
    fclose(fid0);
    free2D(group_skill_mask, param.nG);
}

int main (int argc, char ** argv) {
//    int array[] = {1, 2, 3, 4, 5, 6, 7};
//    NDAT lim = 0;
//    int v = 5;
//    NDAT ix = binsearch(&v, array, 7, sizeof(int)); // SparseArray2D<int>::binsearch(&v, array, 7, &lim);
//    fprintf(stderr,"%d found=%d, lim=%d\n",v, ix, lim);
//    
//    v = 8;
//    ix = SparseArray2D<int>::binsearch(&v, array, 7, &lim);
//    fprintf(stderr,"%d found=%d, lim=%d\n",v, ix, lim);
//    
//    v = 0;
//    ix = SparseArray2D<int>::binsearch(&v, array, 7, &lim);
//    fprintf(stderr,"%d found=%d, lim=%d\n",v, ix, lim);
//    
//    int c = (unsigned long)ceil((double)34600/20000);
    
	clock_t tm_all = clock();//overall time //SEQ
//    double _tm_all = omp_get_wtime(); //PAR
    
	char input_file[1024]; // data
	char output_file[1024]; // model
    char colsole_file[1024]; // console copy
	char predict_file[1024]; // predictions
    
	set_param_defaults(&param);
    
    
    // parse parameters, step 1
	parse_arguments_step1(argc, argv, input_file, output_file, predict_file, colsole_file);

    FILE *fid_console = NULL;
    if(param.duplicate_console==1)
        fid_console = fopen(colsole_file,"w");
    
    if(!param.quiet) {
        printf("trainhmm starting...\n");
        if(param.duplicate_console==1) fprintf(fid_console, "trainhmm starting...\n");
    }

    clock_t tm_read = clock();//overall time //SEQ
//    double _tm_read = omp_get_wtime(); //PAR
    int read_ok = read_and_structure_data(input_file, fid_console);
    tm_read = (clock_t)(clock()-tm_read);//SEQ
    
//    _tm_read = omp_get_wtime()-_tm_read;//PAR
    
//    // experimental
//    for(NCAT k=0; k<param.nK; k++) {
//        char fn [256];
//        sprintf (fn, "b89_kts_skill%d_A.txt", k);
//        InputUtil::writeInputMatrix(fn, &param, param.k_numg[k], param.k_g_data[k]);
//    }

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
//    param.item_complexity = Calloc(NUMBER, param.map_step_bwd->size());
//    line = (char *)malloc(max_line_length);// Malloc(char,max_line_length);
//    while( readline(fid)!=NULL) {
//		// Group
//		col = strtok(line,"\t\n\r");
//		item = string( col );
//		it = param.map_step_fwd->find(item);
//		if( it==param.map_step_fwd->end() ) { // not found
//            fprintf(stderr,"DID NOT FIND THE STEP!!\n");
//            return false;
//		}
//		else {
////            if( it->second > param.map_step_bwd->size()) {
////                int z = 0 ;
////            }
//            NUMBER v = atof( strtok(NULL,"\t\n\r") );
//            param.item_complexity[ it->second ] = v;
//		}
//    }
//    free(line);
//    fclose(fid);

    if(!param.quiet) {
        printf("input read, nO=%d, nG=%d, nK=%d, nI=%d, nZ=%d\n",param.nO, param.nG, param.nK, param.nI, param.nZ);
        if(param.duplicate_console==1) fprintf(fid_console, "input read, nO=%d, nG=%d, nK=%d, nI=%d, nZ=%d\n",param.nO, param.nG, param.nK, param.nI, param.nZ);
    }
    
    //    // write time
    //    if(param.time==1) {
    //        const char * fn = "a89_kts_times.txt";
    ////        const char * fn = "a89_uskts_times.txt";
    //        write_time_interval_data(&param, fn);
    //    }
    
    // erase blocking labels
    zeroLabels(&param);

    clock_t tm_fit; //SEQ
    clock_t tm_predict; //SEQ
//    double _tm_fit;//PAR
//    double _tm_predict;//PAR
    
    if(param.cv_folds==0) { // not cross-validation
        // create problem
        HMMProblem *hmm = NULL;
        switch(param.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmm = new HMMProblem(&param);
                break;
                //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
                //                hmm = new HMMProblemPiG(&param);
                //                break;
            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
                hmm = new HMMProblemSliced(&param);
                break;
            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
                hmm = new HMMProblemPiGK(&param);
                break;
            case STRUCTURE_PIgkww: // Gradient Descent, pLo=f(K,G), other by K
                hmm = new HMMProblemPiGKww(&param);
                break;
            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
                hmm = new HMMProblemPiAGK(&param);
                break;
            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
                hmm = new HMMProblemAGK(&param);
                break;
//                case STRUCTURE_Agki: // Gradient Descent, pT=f(K,G), other by K
//                    hmm = new HMMProblemAGKi(&param);
//                    break;
            case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
                hmm = new HMMProblemPiABGK(&param);
                break;
                //            case BKT_GD_T: // Gradient Descent with Transfer
                //                hmm = new HMMProblemKT(&param);
                //                break;
        }
        tm_fit = clock(); //SEQ
//        _tm_fit = omp_get_wtime(); //PAR
        hmm->fit();
        tm_fit = clock()-tm_fit;//SEQ
//        _tm_fit = omp_get_wtime()-_tm_fit;//PAR
        
        // write model
        hmm->toFile(output_file);
        
        if(param.metrics>0 || param.predictions>0) {
            NUMBER* metrics = Calloc(NUMBER, (size_t)7); // LL, AIC, BIC, RMSE, RMSEnonull, Acc, Acc_nonull;
            // takes care of predictions and metrics, writes predictions if param.predictions==1
            
            // temporary
            if(param.per_kc_rmse_acc) {
                param.kc_counts = Calloc(NDAT, (size_t)param.nK);
                param.kc_rmse = Calloc(NUMBER, (size_t)param.nK);
                param.kc_acc  = Calloc(NUMBER, (size_t)param.nK);
            }

            // NUMBER l1 = hmm->getSumLogPOPara(param.nSeq, param.k_data);
//            printf("hmm-style ll_no_null %15.7f\n",l1);
            
            tm_predict = clock(); //SEQ
//            _tm_predict = omp_get_wtime(); //PAR
//            hmm->predict(metrics, predict_file, param.dat_obs, param.dat_group, param.dat_skill, param.dat_multiskill, false/*all, not only unlabelled*/);
            hmm->predict(metrics, predict_file, param.dat_obs, param.dat_group, param.dat_skill, param.dat_skill_stacked, param.dat_skill_rcount, param.dat_skill_rix, false/*all, not only unlabelled*/);
            
            tm_predict = clock()-tm_predict;//SEQ
//            _tm_predict = omp_get_wtime()-_tm_predict;//PAR
            
            if( param.metrics>0 /*&& !param.quiet*/) {
                printf("trained model LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
                       metrics[0], metrics[1], // ll's
                       2*hmm->getNparams() + 2*metrics[0], hmm->getNparams()*safelog(param.N) + 2*metrics[0],
                       metrics[2], metrics[3], // rmse's
                       metrics[4], metrics[5]); // acc's
                if(param.duplicate_console==1) fprintf(fid_console, "trained model LL=%15.7f (%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",
                       metrics[0], metrics[1], // ll's
                       2*hmm->getNparams() + 2*metrics[0], hmm->getNparams()*safelog(param.N) + 2*metrics[0],
                       metrics[2], metrics[3], // rmse's
                       metrics[4], metrics[5]); // acc's

            }
            free(metrics);

            // temporary
            if(param.per_kc_rmse_acc) {
                for(NCAT i=0; i<param.nK; i++) {
                    printf("KC %4u RMSE=%8.6f Acc=%8.6f\n",i,param.kc_rmse[i],param.kc_acc[i]);
                    if(param.duplicate_console==1) fprintf(fid_console, "KC %4u RMSE=%8.6f Acc=%8.6f\n",i,param.kc_rmse[i],param.kc_acc[i]);
                }
                free(param.kc_counts);
                free(param.kc_rmse);
                free(param.kc_acc);
            }
        
        } // if predict or metrics
        
        delete hmm;
    } else { // cross-validation
        NUMBER* metrics = Calloc(NUMBER, (size_t)7); // AIC, BIC, RMSE, RMSE no null
        switch (param.cv_strat) {
            case CV_GROUP:
                cross_validate(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                cross_validate(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            case CV_ITEM:
                cross_validate_item(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                cross_validate_item(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            case CV_NSTR:
                cross_validate_nstrat(metrics, predict_file, output_file, &tm_fit, &tm_predict, fid_console);//SEQ
//                cross_validate_nstrat(metrics, predict_file, output_file, &_tm_fit, &_tm_predict, fid_console);//PAR
                break;
            default:
                
                break;
        }
//        if(!param.quiet) { // allow to report c-v result anyway
            printf("%d-fold cross-validation: LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",param.cv_folds, metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6]); //SEQ
            if(param.duplicate_console==1) fprintf(fid_console, "%d-fold cross-validation: LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",param.cv_folds, metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6]); //SEQ

//            printf("%d-fold cross-validation: LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",param.cv_folds, metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6]); //PAR
//            if(param.duplicate_console==1) fprintf(fid_console,"%d-fold cross-validation: LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",param.cv_folds, metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6]); //PAR
//        }
        
        free(metrics);
    }
	// free data
	destroy_input_data(&param);
	
//	if(param.quiet == 0) {
        printf("timing: overall %f seconds, read %f, fit %f, predict %f\n",(NUMBER)((clock()-tm_all)/CLOCKS_PER_SEC), (NUMBER)tm_read/CLOCKS_PER_SEC,  (NUMBER)tm_fit/CLOCKS_PER_SEC,  (NUMBER)tm_predict/CLOCKS_PER_SEC);//SEQ
        if(param.duplicate_console==1) fprintf(fid_console, "timing: overall %f seconds, read %f, fit %f, predict %f\n",(NUMBER)((clock()-tm_all)/CLOCKS_PER_SEC), (NUMBER)tm_read/CLOCKS_PER_SEC,  (NUMBER)tm_fit/CLOCKS_PER_SEC,  (NUMBER)tm_predict/CLOCKS_PER_SEC);//SEQ
//        if(param.duplicate_console==1) fprintf(fid_console, "timing: overall %lf sec, read %lf sec, fit %lf sec, predict %lf sec\n",omp_get_wtime()-_tm_all, _tm_read, _tm_fit, _tm_predict);//PAR
//        printf("timing: overall %lf sec, read %lf sec, fit %lf sec, predict %lf sec\n",omp_get_wtime()-_tm_all, _tm_read, _tm_fit, _tm_predict);//PAR
//    }
    
    if(param.duplicate_console==1)
        fclose(fid_console);
    return 0;
}

void exit_with_help() {
	printf(
		   "Usage: trainhmm [options] input_file [[output_file] predicted_response_file]\n"
           "options:\n"
           "-s : structure.solver[.solver setting], structures: 1-by skill, 2-by user;\n"
           "     solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent;\n"
           "     Conjugate Gradient Descent has 3 settings: 1-Polak-Ribiere,\n"
           "     2-Fletcher–Reeves, 3-Hestenes-Stiefel.\n"
           "     For example '-s 1.3.1' would be by skill structure (classical) with\n"
           "     Conjugate Gradient Descent and Hestenes-Stiefel formula, '-s 2.1' would be\n"
           "     by student structure fit using Baum-Welch method.\n"
           "-e : tolerance of termination criterion (0.01 default)\n"
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
           "     For standard BKT - 4 comma-separated numbers: C weight of the penalty and "
           "     centroids, for PI, A, and B matrices respectively. If used for iBKT with\n"
           "     student effects, 8 values will be used with 4 additional values for student\n"
           "     effect matrices. For example, '-c 1.0,0.5,0.5,0.0'."
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
				param.tol = atof(argv[i]);
				if(param.tol<0) {
					fprintf(stderr,"ERROR! Fitting tolerance cannot be negative\n");
					exit_with_help();
				}
				if(param.tol>10) {
					fprintf(stderr,"ERROR! Fitting tolerance cannot be >10\n");
					exit_with_help();
				}
				break;
			case 't':
				param.sliced = (NPAR)atoi(argv[i]);
				if(param.sliced!=0 && param.sliced!=1) {
					fprintf(stderr,"ERROR! Time parameter should be either 0 (off) or 1(om)\n");
					exit_with_help();
				}
				break;
			case 'i':
				param.maxiter = atoi(argv[i]);
				if(param.maxiter<10) {
					fprintf(stderr,"ERROR! Maximum iterations should be at least 10\n");
					exit_with_help();
				}
				break;
			case 'q':
				param.quiet = (NPAR)atoi(argv[i]);
				if(param.quiet!=0 && param.quiet!=1) {
					fprintf(stderr,"ERROR! Quiet param should be 0 or 1\n");
					exit_with_help();
				}
				break;
			case 'n':
				param.nS = (NPAR)atoi(argv[i]);
				if(param.nS<2) {
					fprintf(stderr,"ERROR! Number of hidden states should be at least 2\n");
					exit_with_help();
				}
                if(param.nS != 2) {
                    param.stat_specd_gt2 = true;
                }
				break;
			case 'S':
				param.scaled = (NPAR)atoi(argv[i]);
				if(param.scaled < 0 || param.scaled > 1) {
					fprintf(stderr,"ERROR! Scaling flag should be either 0 (off) or 1 (in)\n");
					exit_with_help();
				}
				break;
			case 's':
				param.structure = (NPAR)atoi( strtok(argv[i],".\t\n\r") );
                ch = strtok(NULL,".\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL)
                    param.solver = (NPAR)atoi(ch);
                ch = strtok(NULL,"\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL)
                    param.solver_setting = (NPAR)atoi(ch);
                if( param.structure != STRUCTURE_SKILL && param.structure != STRUCTURE_GROUP &&
                   param.structure != STRUCTURE_PIg   && param.structure != STRUCTURE_PIgk  &&
                   param.structure != STRUCTURE_PIAgk && param.structure != STRUCTURE_Agk &&
                   param.structure != STRUCTURE_PIABgk && param.structure != STRUCTURE_Agki &&
                   param.structure != STRUCTURE_PIgkww && param.structure != STRUCTURE_SKABslc) {
                    fprintf(stderr, "Model Structure specified (%d) is out of range of allowed values\n",param.structure);
					exit_with_help();
                }
                if( param.solver != METHOD_BW  && param.solver != METHOD_GD &&
                   param.solver != METHOD_CGD && param.solver != METHOD_GDL &&
                   param.solver != METHOD_GBB) {
                    fprintf(stderr, "Method specified (%d) is out of range of allowed values\n",param.solver);
					exit_with_help();
                }
                if( param.structure == STRUCTURE_SKABslc && param.solver == METHOD_BW) {
                    fprintf(stderr, "Method specified (%d) is not defined for this structure (%d) \n",param.solver,param.structure);
                    exit_with_help();
                }
                if( param.solver == METHOD_BW && ( param.solver != STRUCTURE_SKILL && param.solver != STRUCTURE_GROUP ) ) {
                    fprintf(stderr, "Baum-Welch solver does not support model structure specified (%d)\n",param.solver);
					exit_with_help();
                }
                if( param.solver == METHOD_CGD  &&
                   ( param.solver_setting != 1 && param.solver_setting != 2 &&
                    param.solver_setting != 3 )
                   ) {
                    fprintf(stderr, "Conjugate Gradient Descent setting specified (%d) is out of range of allowed values\n",param.solver_setting);
					exit_with_help();
                }
				break;
            case 'f':
                param.single_skill = (NPAR)atoi(argv[i]);
                break;
			case 'm':
                param.metrics = atoi( strtok(argv[i],",\t\n\r"));
                ch = strtok(NULL, "\t\n\r");
                if(ch!=NULL)
                    param.metrics_target_obs = atoi(ch)-1;
				if(param.metrics<0 || param.metrics>1) {
					fprintf(stderr,"value for -m should be either 0 or 1.\n");
					exit_with_help();
				}
				if(param.metrics_target_obs<0) {// || param.metrics_target_obs>(param.nO-1)) {
					fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",param.metrics_target_obs+1);
					exit_with_help();
				}
                break;
			case 'b':
                param.binaryinput = atoi( strtok(argv[i],"\t\n\r"));
                break;
			case 'v':
				param.cv_folds   = (NPAR)atoi( strtok(argv[i],",\t\n\r"));
                ch2 = strtok(NULL, ",\t\n\r");
                if(ch2!=NULL)
                    param.cv_strat = ch2[0];
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    param.cv_target_obs = (NPAR)(atoi(ch)-1);
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    strcpy(param.cv_folds_file, ch);
                ch = strtok(NULL, ",\t\n\r");
                if(ch!=NULL)
                    param.cv_inout_flag = ch[0];
                
				if(param.cv_folds<2) {
					fprintf(stderr,"number of cross-validation folds should be at least 2\n");

					exit_with_help();
				}
				if(param.cv_folds>10) {
					fprintf(stderr,"please keep number of cross-validation folds less than or equal to 10\n");
					exit_with_help();
				}
                if(param.cv_strat != CV_GROUP && param.cv_strat != CV_ITEM && param.cv_strat != CV_NSTR){
					fprintf(stderr,"cross-validation stratification parameter '%c' is illegal\n",param.cv_strat);
					exit_with_help();
                }
				if(param.cv_target_obs<0) {// || param.cv_target_obs>(param.nO-1)) {
					fprintf(stderr,"target observation to be cross-validated against cannot be '%d'\n",param.cv_target_obs+1);
					exit_with_help();
				}
                if( param.cv_inout_flag!='i' && param.cv_inout_flag!='o') {
					fprintf(stderr,"cross-validation folds input/output flag should be ither 'o' (out) or 'i' (in), while it is '%c'\n",param.cv_inout_flag);
					exit_with_help();
                }
                
				break;
            case  'p':
				param.predictions = atoi(argv[i]);
				if(param.predictions<0 || param.predictions>2) {
					fprintf(stderr,"a flag of whether to report predictions for training data (-p) should be 0, 1 or 2\n");
					exit_with_help();
				}
                break;
            case  'd':
				param.multiskill = argv[i][0]; // just grab first character (later, maybe several)
                break;
            case  'P':
				n = atoi(argv[i]);
                if(n!=0 && n!=1 && n!=2) {
					fprintf(stderr,"parallel processing flag (-P) should be 0 or 1\n");
					exit_with_help();
                }
                param.parallel = (NPAR)n;
                break;
			case 'r': // coordinate descend parameters
                // if two first_iteration_qualify,iterations_to_qualify
                // if one iterations_to_qualify (first_iteration_qualify==0)
				n = atoi( strtok(argv[i],",\t\n\r") );
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch==NULL) {// one parameter
                    param.first_iteration_qualify = 0;
                    param.iterations_to_qualify   = (NPAR)n;
                } else {
                    param.first_iteration_qualify = (NPAR)n;
                    param.iterations_to_qualify   = (NPAR)atoi(ch);
                }
				break;
            case 'c': {
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
                    param.Cslices = (NPAR) tmp_array->getSize() / 4;
                    param.Cw = Calloc(NUMBER, (size_t)param.Cslices);
                    param.Ccenters = Calloc(NUMBER, (size_t)(param.Cslices * 3) );
                    int c1 = 0, c2 = 0, i = 0;
                    for(int l=0; l<(int)tmp_array->getSize() / 4; l++) {
                        param.Cw[c1++] = tmp_array->get((NDAT)i++);
                        for(int j=0; j<3; j++)
                            param.Ccenters[c2++] = tmp_array->get((NDAT)i++);
                    }
                    delete tmp_array;
                }
                break;
            case 'o':
                param.duplicate_console = 1;
                strcpy(console_file_name,argv[i]);
                break;
            case '0':
                param.init_reset = true;
                break;
            case 'l': // just to keep it a valid option
                break;
            case 'u': // just to keep it a valid option
                break;
            case 'B': // just to keep it a valid option
                break;
			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
    // post-process checks
    // -v and -m collision
    if(param.cv_folds>0 && param.metrics>0) { // correct for 0-start coding
        fprintf(stderr,"values for -v and -m cannot be both non-zeros\n");
        exit_with_help();
    }
    // scaling
    if(param.scaled == 1 && param.solver != METHOD_BW) {
        param.scaled = 0;
        printf("Scaling can only be enabled for Baum-Welch method. Setting it to off\n");
    }
    // specifying >2 states via -n and mandatory specification of -0 (initial parameters)
    if(param.nS > 2 && !param.init_reset) {
        fprintf(stderr,"when >2 latent states specified via '-n', initial values of parameters have to be explicitly set via '-0'!\n");
        exit_with_help();
    }
    // STRUCTURE_SKABslc solver and -t 1 should be set together
    if( (param.sliced==1) != (param.structure == STRUCTURE_SKABslc) ) {
        fprintf(stderr,"Error! sliced parameter ('-t 1') and STRUCTURE_SKABslc structure should be either both set on or off.\n");
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
    int n;
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
                        strcpy(param.initfile, argv[i]);
                        break;
                    }
                }
                if(param.initfile[0]==0) { // init parameters parameters
                    // init params
                    if(param.init_params!=NULL) free(param.init_params);
                    param.init_params = Calloc(NUMBER, (size_t)n);
                    // read params and write to params
                    param.init_params[0] = atof( strtok(argv[i],",\t\n\r") );
                    for(int j=1; j<n; j++) {
                        param.init_params[j] = atof( strtok(NULL,",\t\n\r") );
                    }
                }
				break;
			case 'l': // lower boundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (NPAR)((argv[i][j]==',')?1:0);
				// init params
				if(param.param_lo!=NULL) free(param.param_lo);
				param.param_lo = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				param.param_lo[0] = atof( strtok(argv[i],",\t\n\r") );
                for(int j=1; j<n; j++) {
					param.param_lo[j] = atof( strtok(NULL,",\t\n\r") );
//                    if(param.param_lo[j] >0) {
//                        int a = 0;
//                    }
                }
                param.lo_lims_specd = true;
				break;
			case 'u': // upper boundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (argv[i][j]==',')?1:0;
				// init params
				if(param.param_hi!=NULL) free(param.param_hi);
				param.param_hi = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				param.param_hi[0] = atof( strtok(argv[i],",\t\n\r") );
                for(int j=1; j<n; j++) {
					param.param_hi[j] = atof( strtok(NULL,",\t\n\r") );
//                    if(param.param_hi[j] < 1) {
//                        int a = 0;
//                    }
                }
                param.hi_lims_specd = true;
				break;
			case 'B': // block fitting
                // first
				param.block_fitting[0] = (NPAR)atoi( strtok(argv[i],",\t\n\r") );
                if(param.block_fitting[0]!=0 && param.block_fitting[0]!=1) {
                    fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                    if(param.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                    exit_with_help();
                }
                // second
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL) {
                    param.block_fitting[1] = (NPAR)atoi(ch);
                    if(param.block_fitting[1]!=0 && param.block_fitting[1]!=1) {
                        fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        if(param.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        exit_with_help();
                    }
                }
                else {
                    fprintf(stderr,"There should be 3 blockig the fitting flags specified.\n");
                    if(param.duplicate_console==1) fprintf(fid_console,"There should be 3 blockig the fitting flags specified.\n");
                    exit_with_help();
                }
                // third
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch != NULL) {
                    param.block_fitting[2] = (NPAR)atoi(ch);
                    if(param.block_fitting[2]!=0 && param.block_fitting[2]!=1) {
                        fprintf(stderr,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        if(param.duplicate_console==1) fprintf(fid_console,"Values of blocking the fitting flags shuld only be 0 or 1.\n");
                        exit_with_help();
                    }
                }
                else {
                    fprintf(stderr,"There should be 3 blockig the fitting flags specified.\n");
                    if(param.duplicate_console==1) fprintf(fid_console,"There should be 3 blockig the fitting flags specified.\n");
                    exit_with_help();
                }
				break;
        } // end switch
    }// end for
    // post parse actions
    
    if(!param.lo_lims_specd && (param.nS!=2 || param.nO!=2) ) { // if not specified, and it's not 2-state 2-obs case, set to 0
        if(param.param_lo!=NULL) free(param.param_lo);
        param.param_lo = Calloc(NUMBER, (size_t)( param.nS*(1+param.nS+param.nO) ) );
    }
 
    if(!param.hi_lims_specd && (param.nS!=2 || param.nO!=2) ) {  // if not specified, and it's not 2-state 2-obs case, set to 1
        if(param.param_hi!=NULL) free(param.param_hi);  // if not specified, set to 1
        param.param_hi = Calloc(NUMBER, (size_t)( param.nS*(1+param.nS+param.nO) ) );
        for(int j=0; j<( param.nS*(1+param.nS+param.nO) ); j++)
            param.param_hi[j] = (NUMBER)1.0;
    }
    
    // post-argument checks - TODO - enable
    if( param.cv_target_obs>(param.nO-1)) {
        fprintf(stderr,"target observation to be cross-validated against cannot be '%d'\n",param.cv_target_obs+1);
        if(param.duplicate_console==1) fprintf(fid_console,"target observation to be cross-validated against cannot be '%d'\n",param.cv_target_obs+1);
        exit_with_help();
    }
    if(param.metrics_target_obs>(param.nO-1)) {
        fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",param.metrics_target_obs+1);
        if(param.duplicate_console==1) fprintf(fid_console,"target observation to compute metrics against cannot be '%d'\n",param.metrics_target_obs+1);
        exit_with_help();
    }
    
}

bool read_and_structure_data(const char *filename, FILE *fid_console) {
    bool readok = true;
    if(param.binaryinput==0)
        readok = InputUtil::readTxt(filename, &param);
    else
        readok = InputUtil::readBin(filename, &param);
    if(! readok )
        return false;
    
	//	2. distribute data into nK skill bins
	//		create
	//          skill_group_map[nK][nG] - explicit 'sparse' map of skills and groups, here 1 means done
	//			k_numg[nK]        - number of groups per skill                 RETAIN
	
	NDAT t = 0;
    NDAT t_stacked = 0;
	NCAT g, k;
//	NPAR o;
	NPAR **skill_group_map = init2D<NPAR>(param.nK, param.nG); // binary map of skills to groups
	param.k_numg = Calloc(NCAT, (size_t)param.nK);
	param.g_numk = Calloc(NCAT, (size_t)param.nG);
    NDAT *count_null_skill_group = Calloc(NDAT, (size_t)param.nG); // count null skill occurences per group
    NCAT *index_null_skill_group = Calloc(NCAT, (size_t)param.nG); // index of group in compressed array
    
	// Pass A
	for(t=0; t<param.N; t++) {
        if(param.multiskill==0)
            k = param.dat_skill[t];//[t];
        else
//            k = param.dat_multiskill->get(t)[1]; // #0 is count, #1 is first element
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ]; // first skill of multi-skill

        g = param.dat_group[t];//[t];
		// null skill : just count
		if( k < 0 ) {
            if(count_null_skill_group[g]==0) param.n_null_skill_group++;
            count_null_skill_group[g]++;
			continue;
		}
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            if( skill_group_map[k][g] == 0 ) {
                skill_group_map[k][g] = 1;
                param.k_numg[k]++;
                param.g_numk[g]++;
            }
        }
	}
    for(k=0; k<param.nK; k++) param.nSeq += param.k_numg[k];
    param.all_data = Calloc(struct data, (size_t)param.nSeq);
    
	// Section B
	param.k_g_data = Malloc(struct data **, (size_t)param.nK);
	param.k_data = Malloc(struct data *, (size_t)param.nSeq);
    //	for(k=0; k<param.nK; k++)
    //		param.k_g_data[k] = Calloc(struct data*, param.k_numg[k]);
	param.g_k_data = Calloc(struct data **, (size_t)param.nG);
	param.g_data = Malloc(struct data *, (size_t)param.nSeq);
    //	for(g=0; g<param.nG; g++)
    //		param.g_k_data[g] = Calloc(struct data*, param.g_numk[g]);
	param.null_skills = Calloc(struct data, (size_t)param.n_null_skill_group);
    // index compressed array of null-skill-BY-group
    NCAT idx = 0;
	for(g=0; g<param.nG; g++)
        if( count_null_skill_group[g] >0 ) index_null_skill_group[g] = idx++;
    
	// Pass C
	NDAT *k_countg = Calloc(NDAT, (size_t)param.nK); // track current group in skill
	NDAT *g_countk = Calloc(NDAT, (size_t)param.nG); // track current skill in group
    // set k_countg and g_countk pointers to relative positions
    NDAT sumk=0, sumg=0;
    for(k=0; k<param.nK; k++) {
        k_countg[k] = sumk;
        param.k_g_data[k] = &param.k_data[sumk];
        sumk += param.k_numg[k];
    }
    for(g=0; g<param.nG; g++) {
        g_countk[g] = sumg;
        param.g_k_data[g] = &param.g_data[sumg];
        sumg += param.g_numk[g];
    }
    
    NDAT n_all_data = 0;
	for(t=0; t<param.N; t++) {
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            g = param.dat_group[t];//[t];
            // now allocate space for the data
            if( k < 0 ) {
                NCAT gidx = index_null_skill_group[g];
                if( param.null_skills[gidx].ix != NULL) // was obs // check if we allocated it already
                    continue;
                param.null_skills[gidx].n = count_null_skill_group[g];
                param.null_skills[gidx].g = g;
                param.null_skills[gidx].cnt = 0;
                //                param.null_skills[gidx].obs = Calloc(NPAR, count_null_skill_group[g]);
                param.null_skills[gidx].ix = Calloc(NDAT, (size_t)count_null_skill_group[g]);
                // no time for null skills is necessary
                //                if(param.time)
                //                    param.null_skills[gidx].time = Calloc(int, (size_t)count_null_skill_group[g]);
                param.null_skills[gidx].alpha = NULL;
                param.null_skills[gidx].beta = NULL;
                param.null_skills[gidx].gamma = NULL;
                param.null_skills[gidx].xi = NULL;
                param.null_skills[gidx].c = NULL;
//                param.null_skills[gidx].time = NULL; // no longer used
                param.null_skills[gidx].p_O_param = 0.0;
                continue;
            }
            if( skill_group_map[k][g]==0) {
                fprintf(stderr, "ERROR! position [%d,%d] in skill_group_map should have been 1\n",k,g);
                if(param.duplicate_console==1) fprintf(fid_console,"ERROR! position [%d,%d] in skill_group_map should have been 1\n",k,g);
            } else if( skill_group_map[k][g]==1 ) { // insert new sequence and grab new data
                // link
                param.k_data[ k_countg[k] ] = &param.all_data[n_all_data]; // in linear array
                param.g_data[ g_countk[g] ] = &param.all_data[n_all_data]; // in linear array
                // param.k_g_data[k] and param.g_k_data[g] are already linked
                //                param.k_g_data[k][ k_countg[k] ] = Calloc(struct data, 1); // grab
                //                param.g_k_data[g][ g_countk[g] ] = param.k_g_data[k][ k_countg[k] ]; // relink
                //                param.k_g_data[k][ k_countg[k] ]->n = 1; // init data << VV
                //                param.k_g_data[k][ k_countg[k] ]->k = k; // init k
                //                param.k_g_data[k][ k_countg[k] ]->g = g; // init g
                //                param.k_g_data[k][ k_countg[k] ]->cnt = 0;
                //                param.k_g_data[k][ k_countg[k] ]->obs = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->alpha = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->beta = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->gamma = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->xi = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->c = NULL;
                //                param.k_g_data[k][ k_countg[k] ]->p_O_param = 0.0;
                //                param.k_g_data[k][ k_countg[k] ]->loglik = 0.0;
                param.all_data[n_all_data].n = 1; // init data << VV
                param.all_data[n_all_data].k = k; // init k
                param.all_data[n_all_data].g = g; // init g
                param.all_data[n_all_data].cnt = 0;
                //                param.all_data[n_all_data].obs = NULL;
                param.all_data[n_all_data].ix = NULL;
                param.all_data[n_all_data].ix_stacked = NULL;
                param.all_data[n_all_data].alpha = NULL;
                param.all_data[n_all_data].beta = NULL;
                param.all_data[n_all_data].gamma = NULL;
                param.all_data[n_all_data].xi = NULL;
                param.all_data[n_all_data].c = NULL;
                param.all_data[n_all_data].p_O_param = 0.0;
                param.all_data[n_all_data].loglik = 0.0;
                k_countg[k]++; // count
                g_countk[g]++; // count
                skill_group_map[k][g] = 2; // mark
                n_all_data++;
            } else if( skill_group_map[k][g]== 2) { // update data count, LINEAR SEARCH :(
                int gidx;
                //                for(gidx=(k_countg[k]-(NCAT)1); gidx>=0 && param.k_g_data[k][gidx]->g!=g; gidx--)
                //                    ;
                for(gidx=(k_countg[k]-1); gidx>=0 && param.k_data[gidx]->g!=g; gidx--)
                    ;
                if( param.k_data[gidx]->g==g)
                    param.k_data[gidx]->n++;
                else {
                    fprintf(stderr, "ERROR! position of group %d in skill %d not found\n",g,k);
                    if(param.duplicate_console==1) fprintf(fid_console,"ERROR! position of group %d in skill %d not found\n",g,k);
                }
            }
        }
	}
	// recycle
	free(k_countg);
	free(g_countk);
    free(count_null_skill_group);
    
    //	3 fill data
    //		pass A
    //			fill k_g_data.data (g_k_data is already linked)
    //				using skill_group_map as marker, 3 - data grabbed
	k_countg = Calloc(NDAT, (size_t)param.nK); // track current group in skill
	g_countk = Calloc(NDAT, (size_t)param.nG); // track current skill in group
	for(t=0; t<param.N; t++) {
		g = param.dat_group[t];
//		o = param.dat_obs->get(t);
//        if(param.sliced)
//            tm = param.dat_slice[t];
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            if( k < 0 ) {
                NCAT gidx = index_null_skill_group[g];
                //                param.null_skills[gidx].obs[ param.null_skills[gidx].cnt++ ] = o; // use .cnt as counter
                param.null_skills[gidx].ix[ param.null_skills[gidx].cnt++ ] = t; // use .cnt as counter
                continue;
            }
            if( skill_group_map[k][g]<2) {
                fprintf(stderr, "ERROR! position [%d,%d] in skill_group_map should have been 2\n",k,g);
                if(param.duplicate_console==1) fprintf(fid_console, "ERROR! position [%d,%d] in skill_group_map should have been 2\n",k,g);
            } else if( skill_group_map[k][g]==2 ) { // grab data and insert first dat point
                //                param.k_g_data[k][ k_countg[k] ]->obs = Calloc(NPAR, (size_t)param.k_g_data[k][ k_countg[k] ]->n); // grab
                //                param.k_g_data[k][ k_countg[k] ]->obs[0] = o; // insert
                param.k_g_data[k][ k_countg[k] ]->ix = Calloc(NDAT, (size_t)param.k_g_data[k][ k_countg[k] ]->n); // grab
                param.k_g_data[k][ k_countg[k] ]->ix[0] = t; // insert
                if(param.multiskill==1) {
                    param.k_g_data[k][ k_countg[k] ]->ix_stacked = Calloc(NDAT, (size_t)param.k_g_data[k][ k_countg[k] ]->n); // grab
                    param.k_g_data[k][ k_countg[k] ]->ix_stacked[0] = t_stacked; // first stacked index
                }
                param.k_g_data[k][ k_countg[k] ]->cnt++; // increase data counter
                k_countg[k]++; // count unique groups forward
                g_countk[g]++; // count unique skills forward
                skill_group_map[k][g] = 3; // mark
            } else if( skill_group_map[k][g]== 3) { // insert next data point and inc counter, LINEAR SEARCH :(
                NCAT gidx;
                //			for(gidx=0; gidx < k_countg[k] && param.k_g_data[k][gidx]->g!=g; gidx++)
                for(gidx=(k_countg[k]-(NCAT)1); gidx>=0 && param.k_g_data[k][gidx]->g!=g; gidx--)
                    ; // skip
                if( param.k_g_data[k][gidx]->g==g ) {
                    NDAT pos = param.k_g_data[k][ gidx ]->cnt; // copy position
                    //                    param.k_g_data[k][ gidx ]->obs[pos] = o; // insert
                    param.k_g_data[k][ gidx ]->ix[pos] = t; // insert
                    if(param.multiskill==1)
                        param.k_g_data[k][ gidx ]->ix_stacked[pos] = t_stacked; // insert
                    param.k_g_data[k][ gidx ]->cnt++; // increase data counter
                }
                else {
                    fprintf(stderr, "ERROR! position of group %d in skill %d not found\n",g,k);
                    if(param.duplicate_console==1) fprintf(fid_console, "ERROR! position of group %d in skill %d not found\n",g,k);
                }
            }
            t_stacked++;
        }// all skills on the row
    } // all t
	// recycle
	free(k_countg);
	free(g_countk);
    free(index_null_skill_group);
    // we delete linear column data in the end now
    //    if(param.cv_folds==0 && param.structure != STRUCTURE_SKILL_T) {
    //        if(param.multiskill==0) {
    //            delete param.dat_skill;//.clear();
    //            param.dat_skill = NULL;
    //        }
    //        else {
    //            delete param.dat_multiskill;
    //            param.dat_multiskill = NULL;
    //        }
    //        delete param.dat_group;
    //        delete param.dat_obs;
    //        delete param.dat_item;
    //        param.dat_group = NULL;
    //        param.dat_obs   = NULL;
    //        param.dat_item  = NULL;
    //    }
	free2D<NPAR>(skill_group_map, param.nK);
    // reset `cnt'
    for(g=0; g<param.nG; g++) // for all groups
        for(k=0; k<param.g_numk[g]; k++) // for all skills in it
            param.g_k_data[g][k]->cnt = 0;
    for(NCAT x=0; x<param.n_null_skill_group; x++)
        param.null_skills[x].cnt = 0;

    return true;
}

void cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//void cross_validate(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    NUMBER rmse = 0.0;
    NUMBER rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    NPAR f;
    NCAT g,k;
    NDAT count=0;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(param.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)param.nG);//[param.nG];
    srand ( (unsigned int)time(NULL) );
    // folds file
    if(param.cv_folds_file[0] > 0) { // file is specified
        if(param.cv_inout_flag=='i') {
            fid_folds = fopen(param.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(param.cv_inout_flag=='o') {
            fid_folds = fopen(param.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(g=0; g<param.nG; g++) {
        if( param.cv_folds_file[0]==0 || (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='o') ) { // output or default
            folds[g] = (NPAR)(rand() % param.cv_folds);
            if(param.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[g]); // write out
        } else if( (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(param.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[g] = (NPAR)(atoi(ch));
        }
    }
    if(param.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(param.cv_inout_flag=='i')
            free(line);
    }

    
    // create and fit multiple problems
    HMMProblem* hmms[param.cv_folds];
    int q = param.quiet;
    param.quiet = 1;
    for(f=0; f<param.cv_folds; f++) {
        switch(param.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblem(&param);
                break;
            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSliced(&param);
                break;
//            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
//                hmm = new HMMProblemPiG(&param);
//                break;
            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
                hmms[f] = new HMMProblemPiGK(&param);
                break;
            case STRUCTURE_PIgkww: // Gradient Descent, pLo=f(K,G), other by K
                hmms[f] = new HMMProblemPiGKww(&param);
                break;
            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
                hmms[f] = new HMMProblemPiAGK(&param);
                break;
            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
                hmms[f] = new HMMProblemAGK(&param);
                break;
            case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
                hmms[f] = new HMMProblemPiABGK(&param);
                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
//                break;
        }
        // block respective data - do not fit the data belonging to the fold
        for(g=0; g<param.nG; g++) // for all groups
            if(folds[g]==f) { // if in current fold
                for(k=0; k<param.g_numk[g]; k++) // for all skills in it
                    param.g_k_data[g][k]->cnt = 1; // block it
            }
        // block nulls
        for(NCAT x=0; x<param.n_null_skill_group; x++) {
            if( param.null_skills[x].g == f)
                param.null_skills[x].cnt = 1;
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
        for(g=0; g<param.nG; g++) // for all groups
            if(folds[g]==f) { // if in current fold
                for(k=0; k<param.g_numk[g]; k++) // for all skills in it
                    param.g_k_data[g][k]->cnt = 0; // UN-block it
            }
        // UN-block nulls
        for(NCAT x=0; x<param.n_null_skill_group; x++) {
            if( param.null_skills[x].g == f)
                param.null_skills[x].cnt = 0;
        }
        if(q == 0) {
            printf("fold %d is done\n",f+1);
            if(param.duplicate_console==1) fprintf(fid_console,"fold %d is done\n",f+1);
        }
    }
    param.quiet = (NPAR)q;
    
    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
    // go trhough original data and predict
	NDAT t;
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1D<NUMBER>(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3D<NUMBER>(param.nG, param.nK, param.nS); // knowledge states //UNBOOST
//   ::boost::numeric::ublas::mapped_matrix<NUMBER*> gsm (param.nG, param.nK);//BOOST
    NUMBER prob = 0, ll = 0;
    struct data dt;

	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs[t];//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group[t];//[t];
        dt.g = g;
        f = folds[g];
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        if(ar[0]<0) { // if no skill label
            rmse += pow(isTarget-hmms[f]->getNullSkillObs(param.cv_target_obs),2);
            accuracy += isTarget == (hmms[f]->getNullSkillObs(param.cv_target_obs)>=0.5);
            
            prob = safe0num(hmms[f]->getNullSkillObs(param.cv_target_obs));
            ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
            if(param.predictions>0) // write predictions file if it was opened
                for(m=0; m<param.nO; m++)
                    fprintf(fid,"%12.10f%s",hmms[f]->getNullSkillObs(m),(m<(param.nO-1))?"\t":"\n");
            continue;
        }

        // check if {g,k}'s were initialized
        for(int l=0; l<n; l++) {
            k = ar[l];
//          NUMBER *z = gsm(g,k); //BOOST
//          if( z==NULL )//BOOST
            if( group_skill_map[g][k][0]==0)//UNBOOST
            {
                dt.k = k;
//                NUMBER * pLbit = Calloc(NUMBER, param.nS);//BOOST
                
                for(i=0; i<param.nS; i++) {
                    group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//UNBOOST
//                    pLbit[i] = hmms[f]->getPI(&dt,i);//BOOST
                    count++;
                }
//              gsm(g,k) = pLbit; //BOOST
            }// pLo/pL not set
        }// for all skills at this transaction
        
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt); //UNBOOST
//        hmms[f]->producePCorrectBoost(&gsm, local_pred, ar, n, &dt); //BOOST

        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
//            NUMBER* pLbit = gsm(g,k); //BOOST
            if(o>-1) { // known observations
                // update p(L)
                pLe_denom = 0.0;
                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
                for(i=0; i<param.nS; i++)
                    pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);  ///// TODO: this is local_pred[o]!!!//UNBOOST
//                  pLe_denom += pLbit[i] * hmms[f]->getB(&dt,i,o); //BOOST
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //UNBOOST
//                  pLe[i] = pLbit[i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //BOOST
                // 2. L = (pLe'*A)';
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i]= 0.0; //UNBOOST
//                  pLbit[i]= 0.0; //BOOST
                for(j=0; j<param.nS; j++)
                    for(j=0; j<param.nS; j++)
                        for(i=0; i<param.nS; i++)
                            group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //UNBOOST
//                          pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //BOOST
            } else { // unknown observation
                // 2. L = (pL'*A)';
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i]; // copy first; //UNBOOST
//                  pLe[i] = pLbit[i]; // copy first; //BOOST
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i] = 0.0; // erase old value //UNBOOST
//                  pLbit[i] = 0.0; // erase old value //BOOST
                for(j=0; j<param.nS; j++)
                    for(i=0; i<param.nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//UNBOOST
//               pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//BOOST
            }// observations
        }
        // write prediction out (after update)
        if(param.predictions>0) { // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t": ((param.predictions==1)?"\n":"\t") );// if we print states of KCs, continut
            if(param.predictions==2) { // if we print out states of KC's as welll
                for(int l=0; l<n; l++) { // all KC here
                    fprintf(fid,"%12.10f%s",group_skill_map[g][ ar[l] ][0], (l==(n-1) && l==(n-1))?"\n":"\t"); // nnon boost // if end of all states: end line//UNBOOST
//           fprintf(fid,"%12.10f%s",gsm(g, ar[l] )[0], (l==(n-1) && l==(n-1))?"\n":"\t"); // if end of all states: end line //BOOST
                }
            }
        }
        
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
        *(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<param.cv_folds; f++) {
        n_par += hmms[f]->getNparams();
        delete hmms[f];
    }
    n_par /= f;
    free(folds);
    free(local_pred);
    free3D<NUMBER>(group_skill_map, param.nG, param.nK);//UNBOOST
//    gsm.clear();//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator1 it1_t;//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator2 it2_t;//BOOST
//    for (it1_t itgsm1 = gsm.begin1(); itgsm1 != gsm.end1(); itgsm1++)//BOOST
//       for (it2_t itgsm2 = itgsm1.begin(); itgsm2 != itgsm1.end(); itgsm2++)//BOOST
//           free( gsm( itgsm2.index1(), itgsm2.index2() ) );//BOOST
    
    
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);
}

void cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//void cross_validate_item(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    NUMBER rmse = 0.0, rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    NCAT I; // item
    NDAT t;
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(param.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)param.nI);
    NDAT *fold_counts = Calloc(NDAT, (size_t)param.cv_folds);
    srand ( (unsigned int)time(NULL) ); // randomize

    // folds file
    if(param.cv_folds_file[0] > 0) { // file is specified
        if(param.cv_inout_flag=='i') {
            fid_folds = fopen(param.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(param.cv_inout_flag=='o') {
            fid_folds = fopen(param.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(I=0; I<param.nI; I++) {
        if( param.cv_folds_file[0]==0 || (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='o') ) { // output or default
            folds[I] = (NPAR)(rand() % param.cv_folds); // produce folds
            if(param.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[I]); // write out
        } else if( (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(param.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[I] = (NPAR)(atoi(ch));
        }
    }
    if(param.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(param.cv_inout_flag=='i')
            free(line);
    }
    
    // count number of items in each fold
    //    for(I=0; I<param.nI; I++) fold_shortcounts[ folds[I] ]++; // produce folds
    for(t=0; t<param.N; t++) fold_counts[ folds[param.dat_item[t]/*->get(t)*/] ]++;
    // create and fit multiple problems
    HMMProblem* hmms[param.cv_folds];
    int q = param.quiet;
    param.quiet = 1;
    for(f=0; f<param.cv_folds; f++) {
        switch(param.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblem(&param);
                break;
                //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
                //                hmm = new HMMProblemPiG(&param);
                //                break;
            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSliced(&param);
                break;
            case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
                hmms[f] = new HMMProblemPiGK(&param);
                break;
            case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
                hmms[f] = new HMMProblemPiAGK(&param);
                break;
            case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
                hmms[f] = new HMMProblemAGK(&param);
                break;
            case STRUCTURE_PIABgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
                hmms[f] = new HMMProblemPiABGK(&param);
                break;
                //            case BKT_GD_T: // Gradient Descent with Transfer
                //                hmm = new HMMProblemKT(&param);
                //                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, (size_t)fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<param.N; t++) {
            if( folds[ param.dat_item[t]/*->get(t)*/ ] == f ) {
                saved_obs[count_saved++] = param.dat_obs[t];
                param.dat_obs[t] = -1;//->set(t, -1);
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
        for(t=0; t<param.N; t++)
            if( folds[ param.dat_item[t]/*->get(t)*/ ] == f )
                param.dat_obs[t]=saved_obs[count_saved++];//->set(t, saved_obs[count_saved++]);
        free(saved_obs);
        if(q == 0) {
            printf("fold %d is done\n",f+1);
            if(param.duplicate_console==1) fprintf(fid_console, "fold %d is done\n",f+1);
        }
    }
    free(fold_counts);
    param.quiet = (NPAR)q;

    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
    // go trhough original data and predict
	NPAR i, j, m, o, isTarget;
    NDAT count=0;
	NUMBER *local_pred = init1D<NUMBER>(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3D<NUMBER>(param.nG, param.nK, param.nS); // knowledge states //UNBOOST
//   ::boost::numeric::ublas::mapped_matrix<NUMBER*> gsm (param.nG, param.nK);//BOOST
    NUMBER prob = 0, ll = 0;
    struct data dt;

	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs[t];//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group[t];//
        dt.g = g;
        I = param.dat_item[t];
        f = folds[I];
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        if(ar[0]<0) { // if no skill label
            rmse += pow(isTarget-hmms[f]->getNullSkillObs(param.cv_target_obs),2);
            accuracy += isTarget == (hmms[f]->getNullSkillObs(param.cv_target_obs)>=0.5);
            
            prob = safe0num(hmms[f]->getNullSkillObs(param.cv_target_obs));
            ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
            if(param.predictions>0) // write predictions file if it was opened
                for(m=0; m<param.nO; m++)
                    fprintf(fid,"%12.10f%s",hmms[f]->getNullSkillObs(m),(m<(param.nO-1))?"\t":"\n");
            continue;
        }
        // check if {g,k}'s were initialized
        for(int l=0; l<n; l++) {
            k = ar[l];
//            NUMBER *z = gsm(g,k); //BOOST
//            if( z==NULL )//BOOST
            if( group_skill_map[g][k][0]==0)//UNBOOST
            {
                dt.k = k;
//                NUMBER * pLbit = Calloc(NUMBER, param.nS);//BOOST

                for(i=0; i<param.nS; i++) {
                    group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//UNBOOST
//                    pLbit[i] = hmms[f]->getPI(&dt,i);//BOOST
                    count++;
                }
//                gsm(g,k) = pLbit; //BOOST
            }// pLo/pL not set
        }// for all skills at this transaction
        
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt); //UNBOOST
//      hmms[f]->producePCorrectBoost(&gsm, local_pred, ar, n, &dt); //BOOST

        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
//          NUMBER* pLbit = gsm(g,k); //BOOST
            if(o>-1) { // known observations
                // update p(L)
                pLe_denom = 0.0;
                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
                for(i=0; i<param.nS; i++)
                    pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);  ///// TODO: this is local_pred[o]!!!//UNBOOST
//                  pLe_denom += pLbit[i] * hmms[f]->getB(&dt,i,o); //BOOST
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //UNBOOST
//                  pLe[i] = pLbit[i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //BOOST
                // 2. L = (pLe'*A)';
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i]= 0.0; //UNBOOST
//                  pLbit[i]= 0.0; //BOOST
                for(j=0; j<param.nS; j++)
                    for(j=0; j<param.nS; j++)
                        for(i=0; i<param.nS; i++)
                            group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //UNBOOST
//                          pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //BOOST
            } else { // unknown observation
                // 2. L = (pL'*A)';
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i]; // copy first; //UNBOOST
//                  pLe[i] = pLbit[i]; // copy first; //BOOST
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i] = 0.0; // erase old value //UNBOOST
//                  pLbit[i] = 0.0; // erase old value //BOOST
                for(j=0; j<param.nS; j++)
                    for(i=0; i<param.nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//UNBOOST
//               pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//BOOST
            }// observations
        }
        // write prediction out (after update)
        if(param.predictions>0) { // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t": ((param.predictions==1)?"\n":"\t") );// if we print states of KCs, continut
            if(param.predictions==2) { // if we print out states of KC's as welll
                for(int l=0; l<n; l++) { // all KC here
                    fprintf(fid,"%12.10f%s",group_skill_map[g][ ar[l] ][0], (l==(n-1) && l==(n-1))?"\n":"\t"); // nnon boost // if end of all states: end line//UNBOOST
//           fprintf(fid,"%12.10f%s",gsm(g, ar[l] )[0], (l==(n-1) && l==(n-1))?"\n":"\t"); // if end of all states: end line //BOOST
                }
            }
        }
        
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
        *(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<param.cv_folds; f++) {
        n_par += hmms[f]->getNparams();
        delete hmms[f];
    }
    n_par /= f;
    free(folds);
    //    free(fold_shortcounts);
    free(local_pred);
    free3D<NUMBER>(group_skill_map, param.nG, param.nK);//UNBOOST
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
//    gsm.clear();//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator1 it1_t;//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator2 it2_t;//BOOST
//    for (it1_t itgsm1 = gsm.begin1(); itgsm1 != gsm.end1(); itgsm1++)//BOOST
//       for (it2_t itgsm2 = itgsm1.begin(); itgsm2 != itgsm1.end(); itgsm2++)//BOOST
//           free( gsm( itgsm2.index1(), itgsm2.index2() ) );//BOOST
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);
}

void cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, clock_t *tm_fit, clock_t *tm_predict, FILE *fid_console) {//SEQ
//void cross_validate_nstrat(NUMBER* metrics, const char *filename, const char *model_file_name, double *tm_fit, double *tm_predict, FILE *fid_console) {//PAR
    NUMBER rmse = 0.0;
    NUMBER rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    NCAT U; // unstratified
    NDAT t;
    NDAT count = 0;
    clock_t tm0;//SEQ
//    double _tm0;//PAR
    char *ch;
    FILE *fid = NULL; // file for storing prediction should that be necessary
    FILE *fid_folds = NULL; // file for reading/writing folds
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr, "Can't write output model file %s\n",filename);
            if(param.duplicate_console==1) fprintf(fid_console, "Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, (size_t)param.N);
    NDAT *fold_counts = Calloc(NDAT, (size_t)param.cv_folds);
    //    NDAT *fold_shortcounts = Calloc(NDAT, (size_t)param.cv_folds);
    
    srand ( (unsigned int)time(NULL) ); // randomize
//    for(I=0; I<param.nI; I++) folds[I] = (NPAR)(rand() % param.cv_folds); // produce folds
    
    // folds file
    if(param.cv_folds_file[0] > 0) { // file is specified
        if(param.cv_inout_flag=='i') {
            fid_folds = fopen(param.cv_folds_file,"r");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for reading %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for reading %s\n",filename);
                exit(1);
            }
            max_line_length = 1024;
            line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
        } else if(param.cv_inout_flag=='o') {
            fid_folds = fopen(param.cv_folds_file,"w");
            if(fid_folds == NULL)
            {
                fprintf(stderr, "Can't open folds file for writing %s\n",filename);
                if(param.duplicate_console==1) fprintf(fid_console, "Can't open folds file for writing %s\n",filename);
                exit(1);
            }
        }
    }
    for(U=0; U<param.N; U++) {
        if( param.cv_folds_file[0]==0 || (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='o') ) { // output or default
            folds[U] = (NPAR)(rand() % param.cv_folds); // produce folds
            if(param.cv_folds_file[0] > 0 )
                fprintf(fid_folds,"%i\n",folds[U]); // write out
        } else if( (param.cv_folds_file[0] > 0 && param.cv_inout_flag=='i') ) {
            readline(fid_folds);//
            ch = strtok(line,"\t\n\r");
            if(ch == NULL) {
                fprintf(stderr, "Error reading input folds file (potentialy wrong number of rows)\n");
                if(param.duplicate_console==1) fprintf(fid_console, "Error reading input folds file (potentialy wrong number of rows)\n");
                exit(1);
            }
            folds[U] = (NPAR)(atoi(ch));
        }
    }
    if(param.cv_folds_file[0] > 0) { // file is specified
        fclose(fid_folds);
        if(param.cv_inout_flag=='i')
            free(line);
    }
    
    
    // count number of items in each fold
    //    for(I=0; I<param.nI; I++) fold_shortcounts[ folds[I] ]++; // produce folds
    for(t=0; t<param.N; t++)  fold_counts[ folds[param.dat_item[t]/*->get(t)*/] ]++;
    // create and fit multiple problems
    HMMProblem* hmms[param.cv_folds];
    int q = param.quiet;
    param.quiet = 1;
    for(f=0; f<param.cv_folds; f++) {
        switch(param.structure)
        {
            case STRUCTURE_SKILL: // Conjugate Gradient Descent
            case STRUCTURE_GROUP: // Conjugate Gradient Descent
                hmms[f] = new HMMProblem(&param);
                break;
            case STRUCTURE_SKABslc: // Conjugate Gradient Descent
                hmms[f] = new HMMProblemSliced(&param);
                break;
//                //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
//                //                hmm = new HMMProblemPiG(&param);
//                //                break;
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
//                //            case BKT_GD_T: // Gradient Descent with Transfer
//                //                hmm = new HMMProblemKT(&param);
//                //                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, (size_t)fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<param.N; t++) {
            if( folds[ param.dat_item[t]/*->get(t)*/ ] == f ) {
                saved_obs[count_saved++] = param.dat_obs[t];
                param.dat_obs[t]=-1;//->set(t, -1);
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
        for(t=0; t<param.N; t++)
            if( folds[ param.dat_item[t]/*->get(t)*/ ] == f )
                param.dat_obs[t]=saved_obs[count_saved++];//->set(t, saved_obs[count_saved++]);
        free(saved_obs);
        if(q == 0) {
            printf("fold %d is done\n",f+1);
            if(param.duplicate_console==1) fprintf(fid_console, "fold %d is done\n",f+1);
        }
    }
    free(fold_counts);
    param.quiet = (NPAR)q;
    
    tm0 = clock();//SEQ
//    _tm0 = omp_get_wtime();//PAR
    // go trhough original data and predict
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1D<NUMBER>(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3D<NUMBER>(param.nG, param.nK, param.nS); // knowledge states //UNBOOST
//    ::boost::numeric::ublas::mapped_matrix<NUMBER*> gsm (param.nG, param.nK);//BOOST
    NUMBER prob = 0, ll = 0;
    struct data dt;

	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs[t];//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group[t];
        dt.g = g;
        f = folds[t]; // folds are done on all N datapoints
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill[t];
            ar = &k;
            n = 1;
        } else {
//            ar = &param.dat_multiskill->get(t)[1];
//            n = param.dat_multiskill->get(t)[0];
            k = param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            ar = &param.dat_skill_stacked[ param.dat_skill_rix[t] ];
            n = param.dat_skill_rcount[t];
        }
        if(ar[0]<0) { // if no skill label
            rmse += pow(isTarget-hmms[f]->getNullSkillObs(param.cv_target_obs),2);
            accuracy += isTarget == (hmms[f]->getNullSkillObs(param.cv_target_obs)>=0.5);
            
            prob = safe0num(hmms[f]->getNullSkillObs(param.cv_target_obs));
            ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
            if(param.predictions>0) // write predictions file if it was opened
                for(m=0; m<param.nO; m++)
                    fprintf(fid,"%12.10f%s",hmms[f]->getNullSkillObs(m),(m<(param.nO-1))?"\t":"\n");
            continue;
        }
        
        // check if {g,k}'s were initialized
        for(int l=0; l<n; l++) {
            k = ar[l];
//          NUMBER *z = gsm(g,k); //BOOST
//          if( z==NULL )//BOOST
            if( group_skill_map[g][k][0]==0)//UNBOOST
            {
                dt.k = k;
//                NUMBER * pLbit = Calloc(NUMBER, param.nS);//BOOST
                for(i=0; i<param.nS; i++) {
                    group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//UNBOOST
//                    pLbit[i] = hmms[f]->getPI(&dt,i);//BOOST
                    count++;
                }
//              gsm(g,k) = pLbit; //BOOST
            }// pLo/pL not set
        }// for all skills at this transaction
        
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt); //UNBOOST
//        hmms[f]->producePCorrectBoost(&gsm, local_pred, ar, n, &dt); //BOOST
        
        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
//          NUMBER* pLbit = gsm(g,k); //BOOST
            if(o>-1) { // known observations
                // update p(L)
                pLe_denom = 0.0;
                // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
                for(i=0; i<param.nS; i++)
                    pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);  // TODO: this is local_pred[o]!!!//UNBOOST
//                  pLe_denom += pLbit[i] * hmms[f]->getB(&dt,i,o); //BOOST
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //UNBOOST
//                  pLe[i] = pLbit[i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom); //BOOST
                // 2. L = (pLe'*A)';
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i]= 0.0; //UNBOOST
//                  pLbit[i]= 0.0; //BOOST
                for(j=0; j<param.nS; j++)
                    for(j=0; j<param.nS; j++)
                        for(i=0; i<param.nS; i++)
                            group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //UNBOOST
//                          pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//A[i][j]; //BOOST
            } else { // unknown observation
                // 2. L = (pL'*A)';
                for(i=0; i<param.nS; i++)
                    pLe[i] = group_skill_map[g][k][i]; // copy first; //UNBOOST
//                  pLe[i] = pLbit[i]; // copy first; //BOOST
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][i] = 0.0; // erase old value //UNBOOST
//                  pLbit[i] = 0.0; // erase old value //BOOST
                for(j=0; j<param.nS; j++)
                    for(i=0; i<param.nS; i++)
                        group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);//UNBOOST
//               pLbit[j] += pLe[i] * hmms[f]->getA(&dt,i,j);//BOOST
            }// observations
        }
        // write prediction out (after update)
        if(param.predictions>0) { // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t": ((param.predictions==1)?"\n":"\t") );// if we print states of KCs, continut
            if(param.predictions==2) { // if we print out states of KC's as welll
                for(int l=0; l<n; l++) { // all KC here
                    fprintf(fid,"%12.10f%s",group_skill_map[g][ ar[l] ][0], (l==(n-1) && l==(n-1))?"\n":"\t"); // nnon boost // if end of all states: end line//UNBOOST
//           fprintf(fid,"%12.10f%s",gsm(g, ar[l] )[0], (l==(n-1) && l==(n-1))?"\n":"\t"); // if end of all states: end line //BOOST
                }
            }
        }
        
        
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
        *(tm_predict) += (clock_t)(clock()- tm0);//SEQ
//    *(tm_predict) += omp_get_wtime()-_tm0;//PAR
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<param.cv_folds; f++) {
        n_par += hmms[f]->getNparams();
        delete hmms[f];
    }
    n_par /= f;
    free(folds);
    //    free(fold_shortcounts);
    free(local_pred);
    free3D<NUMBER>(group_skill_map, param.nG, param.nK);//UNBOOST
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
//    gsm.clear();//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator1 it1_t;//BOOST
//    typedef boost::numeric::ublas::mapped_matrix<NUMBER *>::iterator2 it2_t;//BOOST
//    for (it1_t itgsm1 = gsm.begin1(); itgsm1 != gsm.end1(); itgsm1++)//BOOST
//       for (it2_t itgsm2 = itgsm1.begin(); itgsm2 != itgsm1.end(); itgsm2++)//BOOST
//           free( gsm( itgsm2.index1(), itgsm2.index2() ) );//BOOST
    
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);
}