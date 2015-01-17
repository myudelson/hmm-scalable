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
#include "HMMProblem.h"
////#include "HMMProblemPiG.h"
//#include "HMMProblemPiGK.h"
//#include "HMMProblemAGK.h"
//#include "HMMProblemPiAGK.h"
//#include "HMMProblemPiABGK.h"
////#include "HMMProblemKT.h"
#include "StripedArray.h"
using namespace std;

#define COLUMNS 4

struct param param;
static char *line = NULL;
static int max_line_length;

void exit_with_help();
void parse_arguments(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name);
bool structure_data(const char *filename);
bool read_data(const char *filename);
static char* readline(FILE *fid);
void cross_validate(NUMBER* metrics, const char *filename);
void cross_validate_item(NUMBER* metrics, const char *filename);
void cross_validate_nstrat(NUMBER* metrics, const char *filename);


int main (int argc, char ** argv) {
	clock_t tm0 = clock();
	char input_file[1024];
	char output_file[1024];
	char predict_file[1024];

	parse_arguments(argc, argv, input_file, output_file, predict_file);
    
    if(!param.quiet)
        printf("trainhmm starting...\n");
	bool is_data_read = structure_data(input_file);

    if( is_data_read ) {
        if(!param.quiet)
            printf("input read, nO=%d, nG=%d, nK=%d, nI=%d\n",param.nO, param.nG, param.nK, param.nI);
        clock_t tm;
        // erase blocking labels
        zeroTags(&param);
        // go
        if(param.cv_folds==0) { // not cross-validation
            // create problem
            HMMProblem *hmm;
            switch(param.structure)
            {
                case STRUCTURE_SKILL: // Conjugate Gradient Descent
                case STRUCTURE_GROUP: // Conjugate Gradient Descent
                    hmm = new HMMProblem(&param);
                    break;
//    //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
//    //                hmm = new HMMProblemPiG(&param);
//    //                break;
//                case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                    hmm = new HMMProblemPiGK(&param);
//                    break;
//                case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                    hmm = new HMMProblemPiAGK(&param);
//                    break;
//                case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
//                    hmm = new HMMProblemAGK(&param);
//                    break;
//                case STRUCTURE_PIABgk: // Gradient Descent, pT=f(K,G), other by K
//                    hmm = new HMMProblemPiABGK(&param);
//                    break;
//    //            case BKT_GD_T: // Gradient Descent with Transfer
//    //                hmm = new HMMProblemKT(&param);
//    //                break;
            }
            
            tm = clock();
            hmm->fit();
            
            //
            // this is "incremental" bit... will be removed later
            //
    //        HMMProblem *hmmO;
            
    //        // (1-2)
    //        hmmO = hmm;
    //        hmm = new HMMProblemPiGK(&param);
    //        param.structure = STRUCTURE_PIgk;
    //        for(NCAT k=0; k<param.nK; k++) { // copy the rest
    //            cpy1DNumber(hmmO->getPI(k), hmm->getPI(k), param.nS);
    //            cpy2DNumber(hmmO->getA(k), hmm->getA(k), param.nS, param.nS);
    //            cpy2DNumber(hmmO->getB(k), hmm->getB(k), param.nS, param.nO);
    //        }
    //        delete hmmO;
    //        hmm->fit();
            
    //        // (1-3)
    //        hmmO = hmm;
    //        hmm = new HMMProblemAGK(&param);
    //        param.structure = STRUCTURE_Agk;
    //        for(NCAT k=0; k<param.nK; k++) { // copy the rest
    //            cpy1DNumber(hmmO->getPI(k), hmm->getPI(k), param.nS);
    //            cpy2DNumber(hmmO->getA(k), hmm->getA(k), param.nS, param.nS);
    //            cpy2DNumber(hmmO->getB(k), hmm->getB(k), param.nS, param.nO);
    //        }
    //        delete hmmO;
    //        hmm->fit();

    //        // (1-4)
    //        hmmO = hmm;
    //        hmm = new HMMProblemPiAGK(&param);
    //        param.structure = STRUCTURE_PIAgk;
    //        for(NCAT k=0; k<param.nK; k++) { // copy the rest
    //            cpy1DNumber(hmmO->getPI(k), hmm->getPI(k), param.nS);
    //            cpy2DNumber(hmmO->getA(k), hmm->getA(k), param.nS, param.nS);
    //            cpy2DNumber(hmmO->getB(k), hmm->getB(k), param.nS, param.nO);
    //        }
    //        delete hmmO;
    //        hmm->fit();
            
    //        // (2-4,1-2-4)
    //        hmmO = hmm;
    //        hmm = new HMMProblemPiAGK(&param);
    //        param.structure = STRUCTURE_PIAgk;
    //        for(NCAT k=0; k<param.nK; k++) { // copy the rest
    //            cpy1DNumber(hmmO->getPI(k), hmm->getPI(k), param.nS);
    //            cpy2DNumber(hmmO->getA(k), hmm->getA(k), param.nS, param.nS);
    //            cpy2DNumber(hmmO->getB(k), hmm->getB(k), param.nS, param.nO);
    //        }
    //        for(NCAT g=0; g<param.nG; g++) { // copy the rest
    //            cpy1DNumber(((HMMProblemPiGK *)hmmO)->getPIg(g), ((HMMProblemPiAGK *)hmm)->getPIg(g), param.nS);
    //        }
    //        delete hmmO;
    //        hmm->fit();

    //        // (3-4, 1-3-4)
    //        hmmO = hmm;
    //        hmm = new HMMProblemPiAGK(&param);
    //        param.structure = STRUCTURE_PIAgk;
    //        for(NCAT k=0; k<param.nK; k++) { // copy the rest
    //            cpy1DNumber(hmmO->getPI(k), hmm->getPI(k), param.nS);
    //            cpy2DNumber(hmmO->getA(k), hmm->getA(k), param.nS, param.nS);
    //            cpy2DNumber(hmmO->getB(k), hmm->getB(k), param.nS, param.nO);
    //        }
    //        for(NCAT g=0; g<param.nG; g++) { // copy the rest
    //            cpy2DNumber(((HMMProblemAGK *)hmmO)->getAg(g), ((HMMProblemPiAGK *)hmm)->getAg(g), param.nS, param.nS);
    //        }
    //        delete hmmO;
    //        hmm->fit();
            
            
            if(param.quiet == 0)
                printf("fitting is done in %8.6f seconds\n",(NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
            
            // write model
            hmm->toFile(output_file);
            
            if(param.metrics>0 || param.predictions>0) {
                fprintf(stderr,"Prediction is not yet implemented");
//                NUMBER* metrics = Calloc(NUMBER, (size_t)7); // LL, AIC, BIC, RMSE, RMSEnonull, Acc, Acc_nonull;
//                if(param.predictions>0) {
//                    hmm->predict(metrics, predict_file, param.dat_obs, param.dat_group, param.dat_skill, param.dat_multiskill, false/*all, not only unlabelled*/);
//                } else {
//                    hmm->computeMetrics(metrics);
//                }
//                if( param.metrics>0 && !param.quiet) {
//                    printf("trained model LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f)\n",metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6]);
//                }
//                free(metrics);
            } // if predict or metrics
            
            delete hmm;
        } else { // cross-validation
            tm = clock();
            NUMBER* metrics = Calloc(NUMBER, (size_t)7); // Loglik, AIC, BIC, RMSE, RMSE no null, Accuracy, Accuracy no null
            switch (param.cv_strat) {
                case CV_GROUP:
                    cross_validate(metrics, predict_file);
                    break;
                case CV_ITEM:
                    cross_validate_item(metrics, predict_file);
                    break;
                case CV_NSTR:
                    cross_validate_nstrat(metrics, predict_file);
                    break;
                default:
                    
                    break;
            }
            if(!param.quiet)
                printf("%d-fold cross-validation: LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f), Acc=%8.6f (%8.6f) computed in %8.6f seconds\n",param.cv_folds, metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6], (NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
            free(metrics);
        }
    }
	// free data
	destroy_input_data(&param);
	
	if(param.quiet == 0)
		printf("overall time running is %8.6f seconds\n",(NUMBER)(clock()-tm0)/CLOCKS_PER_SEC);
    return 0;
}

void exit_with_help() {
	printf(
		   "Usage: trainhmm [options] input_file [[output_file] predicted_response_file]\n"
		   "options:\n"
		   "-s : structure.solver[.solver setting], structures: 1-by skill, 2-by user,\n"
           "     3-Pi by skill, A,B-by user, 4-Pi by skill and user, A,B-by user,\n"
           "     5-A by skill and user, Pi,B by skill, 6-Pi,A by skill and user, B by skill;\n"
           "     solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent;\n"
           "     Conjugate Gradient Descent has 3 setings: 1-Polak-Ribiere, 2-Fletcher–Reeves,\n"
           "     3-Hestenes-Stiefel.\n"
           "     For example '-s 1.3.1' would be by skill structure (classical) with Conjugate\n"
           "     Gradient Descent and Hestenes-Stiefel formula, '-s 2.1' would be by student structure\n"
           "     fit using Baum-Welch method.\n"
		   "-t : tolerance of termination criterion (0.01 default)\n"
		   "-i : maximum iterations (200 by default)\n"
		   "-q : quiet mode, without output, 0-no (default), or 1-yes\n"
		   "-n : number of hidden states, should be 2 or more (default 2)\n"
		   "-0 : initial params (comma-separated), {PI{n-1},A{(n*(n-1)},B{n*(observations-1)}}\n"
		   "     default 0.5,1.0,0.4,0.8,0.2\n"
		   "-l : lower boundaries for params (comma-separated), {PI{n},A{n*n},B{n*(observations)}}\n"
		   "     default 0,0,1,0,0,0,0,0,0,0\n"
		   "-u : upper boundaries for params (comma-separated), {PI{n},A{n*n},B{n*(observations)}}\n"
		   "     default 0,0,1,0,0,0,0,0,0,0\n"
		   "-c : weight of the L1 penalty, 0 (default) if not applicable, 1 if applicable\n"
		   "-f : fit as one skill, 0-no (default), 1 - fit as one skill and use params as starting\n"
           "     point for multi-skill, 2 - force one skill"
		   "-m : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To specify observation\n"
           "     for which metrics to be reported, list it after ';'. For example '-m 0', '-m 1' (by default, \n"
           "     observation 1 is assumed), '-m 1,2' (compute metrics for observation 2). Incompatible with-v option.\n"
		   "-v : cross-validation folds and target state to validate against, perform subject-stratified\n"
		   "     cross-validation, default 0 (no cross-validation),\n"
           "     examples '-v 5,2' - 5 fold, predict state 2, '-v 10' - 10-fold predict state 1 by default.\n"
		   "-p : report model predictions on the train set 0-no (default), 1-yes; works with any combination of\n"
           "     -v and -m params.\n"
		   "-d : multi-skill per observation delimiter 0-sinle skill per observation (default), [delimiter character].\n"
           "     For example, '-d ~'.\n"
		   "-b : treat input file as binary input file (lookup format specifications in help)\n"
		   );
	exit(1);
}

static char* readline(FILE *fid) {
	int length = 0;
	
	if(fgets(line,max_line_length,fid) == NULL)
		return NULL;
	
	while(strrchr(line,'\n') == NULL)
	{
		max_line_length *= 2;
		line = (char *) realloc(line,max_line_length);
		length = (int) strlen(line);
		if(fgets(line+length,max_line_length-length,fid) == NULL)
			break;
	}
	return line;
}

void parse_arguments(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name) {
    set_param_defaults(&param);
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
	int i;
    NPAR n;
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
				param.time = atof(argv[i]);
				if(param.time!=0 && param.time!=1) {
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
				param.quiet = atoi(argv[i]);
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
				//fprintf(stdout, "fit single skill=%d\n",param.quiet);
				break;
			case 's':
//				param.solver = (NPAR)atoi( strtok(argv[i],".\t\n\r") );
//                ch = strtok(NULL,"\t\n\r");
//                if(ch != NULL)
//                    param.solver_setting = (NPAR)atoi(ch);
//                if( param.solver != BKT_CGD      && param.solver != BKT_GD      &&
//                    param.solver != BKT_BW       && param.solver != BKT_GD_BW   &&
//                    param.solver != BKT_BW_GD    && param.solver != BKT_GD_G    &&
//                    param.solver != BKT_GD_PIg   && param.solver != BKT_GD_PIgk &&
//                    param.solver != BKT_GD_APIgk && param.solver != BKT_GD_Agk  &&
//                    param.solver != BKT_GD_T ) {
//                    fprintf(stderr, "Method specified (%d) is out of range of allowed values\n",param.solver);
//					exit_with_help();
//                }
//                // new
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
                    param.structure != STRUCTURE_PIABgk    ) {
                    fprintf(stderr, "Model Structure specified (%d) is out of range of allowed values\n",param.structure);
					exit_with_help();
                }
                if( param.solver != METHOD_BW  && param.solver != METHOD_GD &&
                    param.solver != METHOD_CGD ) {
                    fprintf(stderr, "Method specified (%d) is out of range of allowed values\n",param.solver);
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
			case '0': // init_params
				int len;
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (argv[i][j]==',')?(NPAR)1:(NPAR)0;
				// init params
				free(param.init_params);
				param.init_params = Calloc(NUMBER, n);
				// read params and write to params
				param.init_params[0] = atof( strtok(argv[i],",\t\n\r") );
				for(int j=1; j<n; j++)
					param.init_params[j] = atof( strtok(NULL,",\t\n\r") );
				break;
			case 'l': // lower poundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (argv[i][j]==',')?(NPAR)1:(NPAR)0;
				// init params
				free(param.param_lo);
				param.param_lo = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				param.param_lo[0] = atof( strtok(argv[i],",\t\n\r") );
				for(int j=1; j<n; j++)
					param.param_lo[j] = atof( strtok(NULL,",\t\n\r") );
				break;
			case 'u': // upper poundaries
				len = (int)strlen( argv[i] );
				// count delimiters
				n = 1; // start with 1
				for(int j=0;j<len;j++)
					n += (argv[i][j]==',')?(NPAR)1:(NPAR)0;
				// init params
				free(param.param_hi);
				param.param_hi = Calloc(NUMBER, (size_t)n);
				// read params and write to params
				param.param_hi[0] = atof( strtok(argv[i],",\t\n\r") );
				for(int j=1; j<n; j++)
					param.param_hi[j] = atof( strtok(NULL,",\t\n\r") );
				break;
			case 'c':
				param.C = atof(argv[i]);
				if(param.C < 0) {
					fprintf(stderr,"Regularization parameter C should be above 0.\n");
					exit_with_help();
				}
				if(param.C > 1000) {
					fprintf(stderr,"Regularization parameter C is _very) high and might be impractical(%f).\n", param.C);
				}
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
                ch = strtok(NULL, "\t\n\r");
                if(ch!=NULL)
                    param.cv_target_obs = (NPAR)(atoi(ch)-1);
                
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
				break;
            case  'p':
				param.predictions = atoi(argv[i]);
				if(param.predictions<0 || param.predictions>1) {
					fprintf(stderr,"a flag of whether to report predictions for training data (-p) should be 0 or 1\n");
					exit_with_help();
				}
                break;
            case  'd':
				param.multiskill = argv[i][0]; // just grab first character (later, maybe several)
                break;
			case 'r': // coordinate descend parameters
                // if two first_iteration_qualify,iterations_to_qualify
                // if one iterations_to_qualify (first_iteration_qualify==0)
				n = (NPAR)atoi( strtok(argv[i],",\t\n\r") );
                ch = strtok(NULL,",\t\n\r"); // could be NULL (default GD solver)
                if(ch==NULL) {// one parameter
                    param.first_iteration_qualify = 0;
                    param.iterations_to_qualify   = n;
                } else {
                    param.first_iteration_qualify = n;
                    param.iterations_to_qualify   = (NPAR)atoi(ch);
                }
				break;
			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
    // post-argument checks
//    if( param.cv_target_obs>(param.nO-1)) {
//        fprintf(stderr,"target observation to be cross-validated against cannot be '%d'\n",param.cv_target_obs+1);
//        exit_with_help();
//    }
//    if(param.metrics_target_obs>(param.nO-1)) {
//        fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",param.metrics_target_obs+1);
//        exit_with_help();
//    }
    if(param.cv_target_obs>0 && param.metrics>0) { // correct for 0-start coding
        fprintf(stderr,"values for -v and -m cannot be both non-zeros\n");
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

bool structure_data(const char *filename) {
    
    bool success;
    if(param.binaryinput==0) {
        success = InputUtil::readTxt(filename, &param);
    } else {
        success = InputUtil::readBin(filename, &param);
    }
    
    // set counts
	param.nG = (NCAT)param.map_group_fwd->size();
	param.nK = (NCAT)param.map_skill_fwd->size();
	param.nI = (NCAT)param.map_step_fwd->size();
	
	NDAT t = 0, q = 0;
//    int tm = 0; // time
	NCAT g, k;
    NCAT *ar;
    bool is_null = false;
	NPAR o;
    int n;
	NPAR **skill_group_map = init2D<NPAR>(param.nK+1, param.nG); // binary map of (skills+1) to groups, +1 for null skill groups
	param.k_nG = Calloc(NCAT, (size_t)param.nK + 1); //  +1 for null skill groups
	param.k_N = Calloc(NDAT, (size_t)param.nK+1); // +1 for null skills
//	param.g_numk = Calloc(NCAT, (size_t)param.nG);
//    NDAT *count_null_skill_group = Calloc(NDAT, (size_t)param.nG); // count null skill occurences per group
//    NCAT *index_null_skill_group = Calloc(NCAT, (size_t)param.nG); // index of group in compressed array
    
    //
	// Pass A: determine what skill (k) - group (g) combinations are there
    //          and how many groups (k) per each skill (k) are there
    //
	for(t=0; t<param.N; t++) {
        // first skill
        if(param.multiskill==0)
            k = param.dat_skill->get(t);
        else
            k = param.dat_multiskill->get(t)[1]; // #0 is count, #1 is first element
		g = param.dat_group->get(t);
		// null skill : just count
        // handle skills
        is_null = false;
		if( k < 0 ) { // Null skill
            is_null = true;
            k = param.nK; // last index
            ar = &k;
            n = 1;
		} else if(param.multiskill==0) { // single skill
            ar = &k; // array of skill ids, default &k
            n = 1; // number of skill ids, default 1
        } else { // multiple skill
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
        }
        for(int l=0; l<n; l++) { // for all skills
            k = ar[l];
            param.k_N[k]++; // count number of observations per skill
            if(!is_null)
                param.NN++; // number of skill observations (counting multi-skill codings multiple times) BUT NOT NULL
            if( skill_group_map[k][g] == 0 ) {
                skill_group_map[k][g] = 1;
                if(param.multiskill) { // set init flag - meaning start of the sequence
                    param.dat_multiskill->get(t)[n+1] = 1; // n'th skill +1 'cos first element is size
                }
                param.k_nG[k]++;
                if(is_null)
                    param.nSeqNull++;
                else
                    param.nSeq++;
            }
        }// all skills on the row
	}// for all rows
    // build k_ix1stSeq - ix of first sequence of the skill, account for nulls too
    // build k_Gs - array of arrays of Gs
    // build nK+1 array structure for t indexes per each skil, nK+1st is null
    NDAT idx = 0;
    NCAT **k_Gs = Calloc(NCAT *, (size_t)param.nK+1); // +1 for nulls
    param.k_t = Calloc(NDAT *, (size_t)param.nK+1); // +1 for nulls
    param.k_ix1stSeq = Calloc(NDAT, (size_t)param.nK+1); // +1 for nulls
    for(k=0; k<(param.nK+1); k++) {
        param.k_ix1stSeq[k] = idx;
        idx += param.k_nG[k];
        // k_Gs
        k_Gs[k] = Calloc(NCAT, (size_t)param.k_nG[k]);
        param.k_t[k]  = Calloc(NDAT, (size_t)param.k_N[k]);
        NCAT idxb = 0;
        for(g=0; g<param.nG; g++) // if not zero, record in the array of Gs for k
            if( skill_group_map[k][g]>0 )
                k_Gs[k][idxb++] = g;
    }
    
    
    // allocate all sequences
    param.all_seq = Calloc(struct data, (size_t)param.nSeq+param.nSeqNull);
    for(q=0; q<param.nSeq; q++)
        initDat( &param.all_seq[q] );
//	param.k_g_data = Malloc(struct data **, param.nK);
//	param.k_data = Malloc(struct data *, param.nSeq);
//	param.g_k_data = Calloc(struct data **, (size_t)param.nG);
//	param.g_data = Malloc(struct data *, param.nSeq);
//	param.null_skills = Malloc(struct data, param.n_null_skill_group);
    // index compressed array of null-skill-BY-group
//    idx = 0;
//	for(g=0; g<param.nG; g++)
//        if( count_null_skill_group[g] >0 ) index_null_skill_group[g] = idx++;
    if(param.multiskill)
        param.dat_multiskill_seq = Calloc(data **, (size_t)param.N);
    else
        param.dat_skill_seq = Calloc(data *, (size_t)param.N);
    // End of Pass A
    
    //
	// Pass B : link pointers to sequences in dat_[multi]skill_seq and count
    //          in each sequence
    //
    
//	NDAT *k_countg = Calloc(NDAT, (size_t)param.nK); // track current group in skill
////	NDAT *g_countk = Calloc(NDAT, (size_t)param.nG); // track current skill in group
//    // set k_countg and g_countk pointers to relative positions
//    NDAT sumk=0;//, sumg=0;
//    for(k=0; k<param.nK; k++) {
//        k_countg[k] = sumk;
//        param.k_g_data[k] = &param.k_data[sumk];
//        sumk += param.k_nG[k];
//    }
////    for(g=0; g<param.nG; g++) {
////        g_countk[g] = sumg;
////        param.g_k_data[g] = &param.g_data[sumg];
////        sumg += param.g_numk[g];
////    }
    
//    NDAT n_all_data = 0;
    NDAT *k_count = Calloc(NDAT, param.nK+1); // count t's in param.k_t arrays
	for(t=0; t<param.N; t++) {
        NCAT *ar;
        data **seq;
        int n;
        is_null = false;
        g = param.dat_group->get(t);//[t];
        if(param.multiskill==0) {
            k = param.dat_skill->get(t); // skills
            ar = &k; // array of skills
            n = 1; // size of aray
            seq = &param.dat_skill_seq[t]; // array of pointers to sequences
        } else {
            ar = &param.dat_multiskill->get(t)[1]; // array of skills
            n = param.dat_multiskill->get(t)[0]; // size of arrays
            param.dat_multiskill_seq[t] = Calloc(data*, n);
            seq = param.dat_multiskill_seq[t]; // array of pointers to sequences
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            // now allocate space for the data
            if( k < 0 ) {
                is_null = true;
                k = param.nK; // nK+1st position in the data
//                NCAT gidx = index_null_skill_group[g];
//                if( param.null_skills[gidx].ix != NULL) // was obs // check if we allocated it already
//                    continue;
//                param.null_skills[gidx].n = count_null_skill_group[g];
////                param.null_skills[gidx].g = g;
//                param.null_skills[gidx].tag = 0;
////                param.null_skills[gidx].obs = Calloc(NPAR, count_null_skill_group[g]);
//                param.null_skills[gidx].ix = Calloc(NDAT, count_null_skill_group[g]);
//                // no time for null skills is necessary
////                if(param.time)
////                    param.null_skills[gidx].time = Calloc(int, count_null_skill_group[g]);
//                param.null_skills[gidx].alpha = NULL;
//                param.null_skills[gidx].beta = NULL;
//                param.null_skills[gidx].gamma = NULL;
//                param.null_skills[gidx].xi = NULL;
//                param.null_skills[gidx].c = NULL;
//                param.null_skills[gidx].time = NULL;
//                param.null_skills[gidx].p_O_param = 0.0;
//                continue;
            }
            // fill in the k_t array
            param.k_t[k][ k_count[k]++ ] = t; // post increment of the counter
            
            if( skill_group_map[k][g]==0)
                printf("ERROR! position [%d,%d] in skill_group_map should have been 1\n",k,g);
            else {
                // get coordinates of the sequence
                NDAT baseK = param.k_ix1stSeq[k]; // first seq of the skill
                NCAT offsetG = (NCAT)binsearch(&g, k_Gs[k], param.k_nG[k]); // next is index which is the ix of Group
                if( offsetG<0 ) {
                    printf("ERROR! skill-group combination [%d,%d] is not found in the `k_Gs`\n",k,g);
                    continue;
                }
                NDAT ix = baseK + (NDAT)offsetG;
                seq[l] = &param.all_seq[ix];
                if( skill_group_map[k][g]==1 ) { // insert new sequence and grab new data
                    seq[l]->n = 1;
                    if(is_null)
                        seq[l]->id = -1;
                    else
                        seq[l]->id = k;
                    skill_group_map[k][g] = 2;
                } else if( skill_group_map[k][g]== 2) {
                    seq[l]->n++;
                }
            }
        }// all skills in step
	}// all rows in the data
    free(k_count);
    // set data->ix1st, just like we set param.k_ix1stSeq (not for null skills)
    idx = 0;
    for(q=0; q<(param.nSeq + param.nSeqNull); q++) { // for nulls as well
        param.all_seq[q].ix1st = idx;
        idx += param.all_seq[q].n;
    }
    if(idx != (param.NN+param.N_null)) {
        fprintf(stderr,"Final count of variable size should be the same as param.NN\n");
        exit(1);
    }
	free2D<NPAR>(skill_group_map, param.nK+1); // +1 is for nulls
    // End of Pass B
    
	// recycle
//	free(k_countg);
//	free(g_countk);
//    free(count_null_skill_group);
		
//	3 fill data
//		pass A
//			fill k_g_data.data (g_k_data is already linked)
//				using skill_group_map as marker, 3 - data grabbed	
//	k_countg = Calloc(NDAT, param.nK); // track current group in skill
//	g_countk = Calloc(NDAT, param.nG); // track current skill in group
    zeroTags(&param); // zero labels
    data * seq;
    NDAT ix;
    param.dat_tix = Calloc(NDAT, param.NN+param.N_null);
	for(t=0; t<param.N; t++) {
		g = param.dat_group->get(t);
		o = param.dat_obs->get(t);
//        if(param.time)
//            tm = param.dat_time->get(t);
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            if( k < 0 ) {
                k = param.nK; // last row for nulls
//                NCAT gidx = index_null_skill_group[g];
////                param.null_skills[gidx].obs[ param.null_skills[gidx].cnt++ ] = o; // use .cnt as counter
//                param.null_skills[gidx].ix[ param.null_skills[gidx].tag++ ] = t; // use .cnt as counter
//                continue;
            }
            NDAT baseK = param.k_ix1stSeq[k]; // first seq of the skill
            NCAT offsetG = (NCAT)binsearch(&g, k_Gs[k], param.k_nG[k]); // next is index which is the ix of Group
            seq = &param.all_seq[ baseK + (NDAT)offsetG ]; // grab the seq
            ix = seq->ix1st + (seq->cnt++);
            param.dat_tix[ix] = t;
            
//            if( skill_group_map[k][g]<2)
//                printf("ERROR! position [%d,%d] in skill_group_map should have been 2\n",k,g);
//            else if( skill_group_map[k][g]==2 ) { // grab data and insert first dat point
////                param.k_g_data[k][ k_countg[k] ]->obs = Calloc(NPAR, param.k_g_data[k][ k_countg[k] ]->n); // grab
////                param.k_g_data[k][ k_countg[k] ]->obs[0] = o; // insert
//                param.k_g_data[k][ k_countg[k] ]->ix = Calloc(NDAT, param.k_g_data[k][ k_countg[k] ]->n); // grab
//                if(param.time)
//                    param.k_g_data[k][ k_countg[k] ]->time = Calloc(int, param.k_g_data[k][ k_countg[k] ]->n); // grab
//                param.k_g_data[k][ k_countg[k] ]->ix[0] = t; // insert
//                if(param.time)
//                    param.k_g_data[k][ k_countg[k] ]->time[0] = tm; // insert
//                param.k_g_data[k][ k_countg[k] ]->tag++; // increase data counter
//                k_countg[k]++; // count unique groups forward
////                g_countk[g]++; // count unique skills forward
//                skill_group_map[k][g] = 3; // mark
//            } else if( skill_group_map[k][g]== 3) { // insert next data point and inc counter, LINEAR SEARCH :(
//                NCAT gidx;
//                //			for(gidx=0; gidx < k_countg[k] && param.k_g_data[k][gidx]->g!=g; gidx++)
//                for(gidx=(k_countg[k]-(NCAT)1); gidx>=0 && param.k_g_data[k][gidx]->g!=g; gidx--)
//                    ; // skip
//                if( param.k_g_data[k][gidx]->g==g ) {
//                    NDAT pos = param.k_g_data[k][ gidx ]->tag; // copy position
////                    param.k_g_data[k][ gidx ]->obs[pos] = o; // insert
//                    param.k_g_data[k][ gidx ]->ix[pos] = t; // insert
//                    if(param.time)
//                        param.k_g_data[k][ gidx ]->time[pos] = tm; // insert
//                    param.k_g_data[k][ gidx ]->tag++; // increase data counter
//                }
//                else
//                    printf("ERROR! position of group %d in skill %d not found\n",g,k);
//            }
        }
    }
    free2D<NCAT>(k_Gs, param.nK+1); // +1 is for null skill labels
    
	// recycle
//	free(k_countg);
//	free(g_countk);
//    free(index_null_skill_group);
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
    // reset `cnt'
    zeroTags(&param);
    return true;
}

void cross_validate(NUMBER* metrics, const char *filename) {
/*    NUMBER rmse = 0.0;
    NUMBER rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    FILE *fid; // file for storing prediction should that be necessary
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr,"Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, param.nG);//[param.nG];
    srand ( time(NULL) );
    for(g=0; g<param.nG; g++) folds[g] = rand() % param.cv_folds;
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
////            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
////                break;
        }
        // block respective data - do not fit the data belonging to the fold
        for(g=0; g<param.nG; g++) // for all groups
            if(folds[g]==f) { // if in current fold
                for(k=0; k<param.g_numk[g]; k++) // for all skills in it
                    param.g_k_data[g][k]->tag = 1; // block it
            }
        // block nulls
        for(NCAT x=0; x<param.n_null_skill_group; x++) {
            if( param.null_skills[x].g == f)
                param.null_skills[x].tag = 1;
        }
        // now compute
        hmms[f]->fit();
        // UN-block respective data
        for(g=0; g<param.nG; g++) // for all groups
            if(folds[g]==f) { // if in current fold
                for(k=0; k<param.g_numk[g]; k++) // for all skills in it
                    param.g_k_data[g][k]->tag = 0; // UN-block it
            }
        // UN-block nulls
        for(NCAT x=0; x<param.n_null_skill_group; x++) {
            if( param.null_skills[x].g == f)
                param.null_skills[x].tag = 0;
        }
        if(q == 0)
            printf("fold %d is done\n",f+1);
    }
    param.quiet = q;
    // go trhough original data and predict
	NDAT t;
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1DNumber(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(param.nG, param.nK, param.nS); // knowledge states
    NUMBER prob = 0, ll = 0;
    struct data dt;
	// initialize
	for(g=0; g<param.nG; g++) {
        dt.g = g;
        f = folds[g];
		for(k=0; k<param.nK; k++) {
            dt.k = k;
			for(i=0; i<param.nO; i++)
                group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//PI[i];
		}
    }
	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs->get(t);//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group->get(t);//[t];
        dt.g = g;
        f = folds[g];
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
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
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt);
        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
//            // produce prediction and copy to result
//            for(m=0; m<param.nO; m++)
//                for(i=0; i<param.nS; i++)
//                    local_pred[m] += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,m);
            // update p(L)
            pLe_denom = 0.0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<param.nS; i++) pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);
            for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0;
            for(j=0; j<param.nS; j++)
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);
        }
//        for(m=0; m<param.nO; m++)
//            local_pred[m] /= n;
        if(param.predictions>0) // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
    
    // delete problems
    NCAT n_par = 0;
    for(f=0; f<param.cv_folds; f++) {
        n_par += hmms[f]->getNparams();
        delete hmms[f];
    }
    n_par /= f;
    free(folds);
    free(local_pred);
    free3DNumber(group_skill_map, param.nG, param.nK);
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);*/
}

void cross_validate_item(NUMBER* metrics, const char *filename) {
/*    NUMBER rmse = 0.0, rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    NCAT I; // item
    NDAT t;
    FILE *fid; // file for storing prediction should that be necessary
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr,"Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, param.nI);
    NDAT *fold_counts = Calloc(NDAT, param.cv_folds);
//    NDAT *fold_shortcounts = Calloc(NDAT, param.cv_folds);
    srand ( time(NULL) ); // randomize
    for(I=0; I<param.nI; I++) folds[I] = rand() % param.cv_folds; // produce folds
    // count number of items in each fold
//    for(I=0; I<param.nI; I++) fold_shortcounts[ folds[I] ]++; // produce folds
    for(t=0; t<param.N; t++) fold_counts[ folds[param.dat_item->get(t)] ]++;
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
//            case STRUCTURE_PIABgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiABGK(&param);
//                break;
////            case BKT_GD_T: // Gradient Descent with Transfer
////                hmm = new HMMProblemKT(&param);
////                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<param.N; t++) {
            if( folds[ param.dat_item->get(t) ] == f ) {
                saved_obs[count_saved++] = param.dat_obs->get(t);
                param.dat_obs->set(t, -1);
            }
        }
        // now compute
        hmms[f]->fit();
        
        // UN-block respective data
        count_saved = 0;
        for(t=0; t<param.N; t++)
            if( folds[ param.dat_item->get(t) ] == f )
                param.dat_obs->set(t, saved_obs[count_saved++]);
        free(saved_obs);
        if(q == 0)
            printf("fold %d is done\n",f+1);
    }
    free(fold_counts);
    param.quiet = q;
    // go trhough original data and predict
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1DNumber(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(param.nG, param.nK, param.nS); // knowledge states
    NUMBER prob = 0, ll = 0;
    struct data dt;
	// initialize
	for(g=0; g<param.nG; g++) {
        dt.g = g;
        f = folds[g];
		for(k=0; k<param.nK; k++) {
            dt.k = k;
			for(i=0; i<param.nO; i++)
                group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//PI[i];
		}
    }
	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs->get(t);//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group->get(t);//[t];
        dt.g = g;
        f = folds[g];
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
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
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt);
        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
            //            // produce prediction and copy to result
            //            for(m=0; m<param.nO; m++)
            //                for(i=0; i<param.nS; i++)
            //                    local_pred[m] += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,m);
            // update p(L)
            pLe_denom = 0.0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<param.nS; i++) pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);
            for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0;
            for(j=0; j<param.nS; j++)
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);
        }
        //        for(m=0; m<param.nO; m++)
        //            local_pred[m] /= n;
        if(param.predictions>0) // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
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
    free3DNumber(group_skill_map, param.nG, param.nK);
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);*/
}

void cross_validate_nstrat(NUMBER* metrics, const char *filename) {
/*    NUMBER rmse = 0.0;
    NUMBER rmse_no_null = 0.0, accuracy = 0.0, accuracy_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    NCAT I; // item
    NDAT t;
    FILE *fid; // file for storing prediction should that be necessary
    if(param.predictions>0) {  // if we have to write the predictions file
        fid = fopen(filename,"w");
        if(fid == NULL)
        {
            fprintf(stderr,"Can't write output model file %s\n",filename);
            exit(1);
        }
    }
    // produce folds
    NPAR *folds = Calloc(NPAR, param.nI);
    NDAT *fold_counts = Calloc(NDAT, param.cv_folds);
    //    NDAT *fold_shortcounts = Calloc(NDAT, param.cv_folds);
    srand ( time(NULL) ); // randomize
    for(I=0; I<param.nI; I++) folds[I] = rand() % param.cv_folds; // produce folds
    // count number of items in each fold
    //    for(I=0; I<param.nI; I++) fold_shortcounts[ folds[I] ]++; // produce folds
    for(t=0; t<param.N; t++)  fold_counts[ folds[param.dat_item->get(t)] ]++;
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
////            case BKT_GD_T: // Gradient Descent with Transfer
////                hmm = new HMMProblemKT(&param);
////                break;
        }
        // block respective data - do not fit the data belonging to the fold
        NPAR *saved_obs = Calloc(NPAR, fold_counts[f]);
        NDAT count_saved = 0;
        for(t=0; t<param.N; t++) {
            if( folds[ param.dat_item->get(t) ] == f ) {
                saved_obs[count_saved++] = param.dat_obs->get(t);
                param.dat_obs->set(t, -1);
            }
        }
        // now compute
        hmms[f]->fit();
        
        // UN-block respective data
        count_saved = 0;
        for(t=0; t<param.N; t++)
            if( folds[ param.dat_item->get(t) ] == f )
                param.dat_obs->set(t, saved_obs[count_saved++]);
        free(saved_obs);
        if(q == 0)
            printf("fold %d is done\n",f+1);
    }
    free(fold_counts);
    param.quiet = q;
    // go trhough original data and predict
	NPAR i, j, m, o, isTarget;
	NUMBER *local_pred = init1DNumber(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(param.nG, param.nK, param.nS); // knowledge states
    NUMBER prob = 0, ll = 0;
    struct data dt;
	// initialize
	for(g=0; g<param.nG; g++) {
        dt.g = g;
        f = folds[g];
		for(k=0; k<param.nK; k++) {
            dt.k = k;
			for(i=0; i<param.nO; i++)
                group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//PI[i];
		}
    }
	// deal with null skill
	for(t=0; t<param.N; t++) {
		o = param.dat_obs->get(t);//[t]; correct: obs 1 (0 code), incorect obs 2 (1 code), hence 1-code is the conversion
        isTarget = (NPAR)(param.cv_target_obs == o);
		g = param.dat_group->get(t);//[t];
        dt.g = g;
        f = folds[g];
        
        NCAT *ar;
        int n;
        if(param.multiskill==0) {
            k = param.dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
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
        hmms[f]->producePCorrect(group_skill_map, local_pred, ar, n, &dt);
        for(int l=0; l<n; l++) {
            k = ar[l];
            dt.k = k;
            //            // produce prediction and copy to result
            //            for(m=0; m<param.nO; m++)
            //                for(i=0; i<param.nS; i++)
            //                    local_pred[m] += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,m);
            // update p(L)
            pLe_denom = 0.0;
            // 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
            for(i=0; i<param.nS; i++) pLe_denom += group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o);
            for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i] * hmms[f]->getB(&dt,i,o) / safe0num(pLe_denom);
            // 2. L = (pLe'*A)';
            for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0;
            for(j=0; j<param.nS; j++)
                for(i=0; i<param.nS; i++)
                    group_skill_map[g][k][j] += pLe[i] * hmms[f]->getA(&dt,i,j);
        }
        //        for(m=0; m<param.nO; m++)
        //            local_pred[m] /= n;
        if(param.predictions>0) // write predictions file if it was opened
            for(m=0; m<param.nO; m++)
                fprintf(fid,"%12.10f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        accuracy += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        accuracy_no_null += isTarget == (local_pred[param.cv_target_obs]>=0.5);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
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
    free3DNumber(group_skill_map, param.nG, param.nK);
    if(param.predictions>0) // close predictions file if it was opened
        fclose(fid);
    metrics[0] = ll;
    metrics[1] = 2*(n_par) + 2*ll;
    metrics[2] = n_par*safelog(param.N) + 2*ll;
    metrics[3] = rmse;
    metrics[4] = rmse_no_null;
    metrics[5] = accuracy / param.N;
    metrics[6] = accuracy_no_null / (param.N - param.N_null);*/
}