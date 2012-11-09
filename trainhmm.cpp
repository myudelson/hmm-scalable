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
//#include "HMMProblemPiG.h"
#include "HMMProblemPiGK.h"
#include "HMMProblemAGK.h"
//#include "HMMProblemPiAGK.h"
//#include "HMMProblemKT.h"
#include "StripedArray.h"
using namespace std;

#define COLUMNS 4

struct param param;
static char *line = NULL;
static int max_line_length;
//StripedArray<NPAR> *dat_obs = NULL;
//StripedArray<NCAT> *dat_group = NULL;
//StripedArray<NCAT> *dat_skill = NULL;
//StripedArray< NCAT* > *dat_multiskill = NULL;

void exit_with_help();
void parse_arguments(int argc, char **argv, char *input_file_name, char *output_file_name, char *predict_file_name);
void read_train_data(const char *filename);
void read_predict_data(const char *filename);
void destroy_input_data();
static char* readline(FILE *fid);
void cross_validate(NUMBER* metrics, const char *filename);

int main (int argc, char ** argv) {
	clock_t tm0 = clock();
	
	printf("trainhmm starting...\n");
	set_param_defaults(&param);
	
	char input_file[1024];
	char output_file[1024];
	char predict_file[1024];
	
	parse_arguments(argc, argv, input_file, output_file, predict_file);
	read_train_data(input_file);

	if(param.quiet == 0)
		printf("input read, nO=%d, nG=%d, nK=%d\n",param.nO, param.nG, param.nK);

    clock_t tm;
    // erase blocking labels
    zeroLabels(&param);
    // go
    if(param.cv_folds==0) {
        // create problem
        HMMProblem *hmm;
        switch(param.solver)
        {
            case BKT_CGD: // Conjugate Gradient Descent
            case BKT_GD: // Gradient Descent
            case BKT_BW: // Expectation Maximization (Baum-Welch)
            case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
            case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            case BKT_GD_G: // Gradient Descent by group
                hmm = new HMMProblem(&param);
                break;
//            case BKT_GD_PIg: // Gradient Descent: PI by group, A,B by skill
//                hmm = new HMMProblemPiG(&param);
//                break;
            case BKT_GD_PIgk: // Gradient Descent, pLo=f(K,G), other by K
                hmm = new HMMProblemPiGK(&param);
                break;
//            case BKT_GD_APIgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmm = new HMMProblemPiAGK(&param);
//                break;
            case BKT_GD_Agk: // Gradient Descent, pT=f(K,G), other by K
                hmm = new HMMProblemAGK(&param);
                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmm = new HMMProblemKT(&param);
//                break;
        }
        
        tm = clock();
        hmm->fit();
        if(param.quiet == 0)
            printf("fitting is done in %8.6f seconds\n",(NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
        
        // write model
        hmm->toFile(output_file);
        
        if(param.metrics>0 || param.predictions>0) {
            NUMBER* metrics = Calloc(NUMBER, 5); // AIC, BIC, RMSE
            if(param.predictions>0) {
                read_predict_data(input_file); // re-read the data
                hmm->predict(metrics, predict_file, param.dat_obs, param.dat_group, param.dat_skill, param.dat_multiskill);
                // now remove the data read
                if(param.multiskill==0) {
                    delete param.dat_skill;//.clear();
                    param.dat_skill = NULL;
                }
                else {
                    delete param.dat_multiskill;
                    param.dat_multiskill = NULL;
                }
                delete param.dat_group;
                delete param.dat_obs;
                param.dat_group = NULL;
                param.dat_obs = NULL;
            } else {
                hmm->computeMetrics(metrics);
            }
            if( param.metrics>0 ) {
                printf("trained model LL=%15.7f(%15.7f), AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f)\n",metrics[0], hmm->getLogLik(), metrics[1], metrics[2], metrics[3], metrics[4]);
            }
            free(metrics);
        } // if predict or metrics
        
        delete hmm;
    } else {
        tm = clock();
        NUMBER* metrics = Calloc(NUMBER, 4); // AIC, BIC, RMSE, RMSE no null
        cross_validate(metrics, predict_file);
        printf("%d-fold cross-validation: AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f) computed in %8.6f seconds\n",param.cv_folds, metrics[0],metrics[1],metrics[2], metrics[3], (NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
        free(metrics);
        // now delete source data
        if(param.multiskill==0) {
            delete param.dat_skill;//.clear();
            param.dat_skill = NULL;
        }
        else {
            delete param.dat_multiskill;
            param.dat_multiskill = NULL;
        }
        delete param.dat_group;
        delete param.dat_obs;
        param.dat_group = NULL;
        param.dat_obs = NULL;
    }

	// free data
	destroy_input_data();
	
	if(param.quiet == 0)
		printf("overall time running is %8.6f seconds\n",(NUMBER)(clock()-tm0)/CLOCKS_PER_SEC);
    return 0;
}

void exit_with_help() {
	printf(
		   "Usage: trainhmm [options] input_file [[output_file] predicted_response_file]\n"
		   "options:\n"
		   "(-s) : solver, 0 - brute-force, 1 - grid-search, 2 - gradient descent (2 default)\n"
           "       3 - expectation maximization (Baum-Welch)\n"
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
		   "(-c) : weight of the L1 penalty, 0 (default) if not applicable, 1 if applicable\n"
		   "(-f) : fit as one skill, 0-no (default), 1 - fit as one skill and use params as starting\n"
           "       point for multi-skill, 2 - force one skill"
		   "-m : report model fitting metrics (AIC, BIC, RMSE) 0-no (default), 1-yes. To specify observation\n"
           "     for which metrics to be reported, list it after ';'. For example '-r 0', '-r 1' (by default, \n"
           "     observation 1 is assumed), '-r1:2' (compute metrics for observation 2). Incompatible with-v option.\n"
		   "-v : cross-validation folds and target state to validate against, perform subject-stratified\n"
		   "     cross-validation, default 0 (no cross-validation),\n"
           "     examples '-v 5;2' - 5 fold, predict state 2, '-v 10' - 10-fold predict state 1 by default.\n"
           "     Incompatible with -r option.\n"
		   "-p : report model predictions on the train set 0-no (default), 1-yes; works with any combination of\n"
           "     -v and -m params.\n"
		   "-d : multi-skill per observation delimiter 0-sinle skill per observation (default), [delimiter character].\n"
           "     For example '-d ~'.\n"
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
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
	int i;
    NPAR n;
    char * ch;
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break; // end of options stop parsing
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 't':
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
				param.solver = (NPAR)atoi( strtok(argv[i],".\t\n\r") );
                ch = strtok(NULL,"\t\n\r");
                if(ch != NULL)
                    param.solver_settting = (NPAR)atoi(ch);
                if( param.solver != BKT_CGD      && param.solver != BKT_GD      &&
                    param.solver != BKT_BW       && param.solver != BKT_GD_BW   &&
                    param.solver != BKT_BW_GD    && param.solver != BKT_GD_G    &&
                    param.solver != BKT_GD_PIg   && param.solver != BKT_GD_PIgk &&
                    param.solver != BKT_GD_APIgk && param.solver != BKT_GD_Agk  &&
                    param.solver != BKT_GD_T ) {
                    fprintf(stderr, "Method specified (%d) is out of range of allowed values\n",param.solver);
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
				param.param_lo = Calloc(NUMBER, n);
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
				param.param_hi = Calloc(NUMBER, n);
				// read params and write to params
				param.param_hi[0] = atof( strtok(argv[i],",\t\n\r") );
				for(int j=1; j<n; j++)
					param.param_hi[j] = atof( strtok(NULL,",\t\n\r") );
				break;
			case 'c':
				param.C = atof(argv[i]);
				break;
			case 'm':
                param.metrics = atoi( strtok(argv[i],";\t\n\r"));
                ch = strtok(NULL, "\t\n\r");
                if(ch!=NULL)
                    param.metrics_target_obs = atoi(ch)-1;
				if(param.metrics<0 || param.metrics>1) {
					fprintf(stderr,"value for -r should be either 0 or 1.\n");
					exit_with_help();
				}
				if(param.metrics_target_obs<0) {// || param.metrics_target_obs>(param.nO-1)) {
					fprintf(stderr,"target observation to compute metrics against cannot be '%d'\n",param.metrics_target_obs+1);
					exit_with_help();
				}
                break;
			case 'v':
				param.cv_folds   = (NPAR)atoi( strtok(argv[i],";\t\n\r"));
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
        fprintf(stderr,"values for -v and -r cannot be both non-zeros\n");
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

void read_train_data(const char *filename) {
	read_predict_data(filename);
	param.nG = (NCAT)param.map_group_fwd->size();
	param.nK = (NCAT)param.map_skill_fwd->size();
	
	//	2. distribute data into nK skill bins
	//		create
	//          skill_group_map[nK][nG] - explicit 'sparse' map of skills and groups, here 1 means done
	//			k_numg[nK]        - number of groups per skill                 RETAIN
	
	NDAT t = 0;
	NCAT g, k;
	NPAR o;
	NPAR **skill_group_map = init2DNCat(param.nK, param.nG); // binary map of skills to groups
	param.k_numg = Calloc(NCAT, param.nK);
	param.g_numk = Calloc(NCAT, param.nG);
    NDAT *count_null_skill_group = Calloc(NDAT, param.nG); // count null skill occurences per group
    NCAT *index_null_skill_group = Calloc(NCAT, param.nG); // index of group in compressed array

	// Pass A
	for(t=0; t<param.N; t++) {
        if(param.multiskill==0)
            k = param.dat_skill->get(t);//[t];
        else
            k = param.dat_multiskill->get(t)[1]; // #0 is count, #1 is first element
		g = param.dat_group->get(t);//[t];
		// null skill : just count
		if( k < 0 ) {
            if(count_null_skill_group[g]==0) param.n_null_skill_group++;
            count_null_skill_group[g]++;
			continue;
		}
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
            if( skill_group_map[k][g] == 0 ) {
                skill_group_map[k][g] = 1;
                param.k_numg[k]++;
                param.g_numk[g]++;
            }
        }
	}
    for(k=0; k<param.nK; k++) param.ndata += param.k_numg[k];
    param.all_data = Calloc(struct data, param.ndata);
		
	// Section B
	param.k_g_data = Malloc(struct data **, param.nK);
	param.k_data = Malloc(struct data *, param.ndata);
//	for(k=0; k<param.nK; k++)
//		param.k_g_data[k] = Calloc(struct data*, param.k_numg[k]);
	param.g_k_data = Calloc(struct data **, param.nG);
	param.g_data = Malloc(struct data *, param.ndata);
//	for(g=0; g<param.nG; g++)
//		param.g_k_data[g] = Calloc(struct data*, param.g_numk[g]);
	param.null_skills = Malloc(struct data, param.n_null_skill_group);
    // index compressed array of null-skill-BY-group
    NCAT idx = 0;
	for(g=0; g<param.nG; g++)
        if( count_null_skill_group[g] >0 ) index_null_skill_group[g] = idx++;
    
	// Pass C
	NDAT *k_countg = Calloc(NDAT, param.nK); // track current group in skill
	NDAT *g_countk = Calloc(NDAT, param.nG); // track current skill in group
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
            k = param.dat_skill->get(t);
            ar = &k;
            n = 1;
        } else {
            ar = &param.dat_multiskill->get(t)[1];
            n = param.dat_multiskill->get(t)[0];
        }
        for(int l=0; l<n; l++) {
            k = ar[l];
            g = param.dat_group->get(t);//[t];
            // now allocate space for the data
            if( k < 0 ) {
                NCAT gidx = index_null_skill_group[g];
                if( param.null_skills[gidx].obs != NULL) // check if we allocated it already
                    continue;
                param.null_skills[gidx].ndat = count_null_skill_group[g];
                param.null_skills[gidx].g = g;
                param.null_skills[gidx].cnt = 0;
                param.null_skills[gidx].obs = Calloc(NPAR, count_null_skill_group[g]);
                param.null_skills[gidx].alpha = NULL;
                param.null_skills[gidx].beta = NULL;
                param.null_skills[gidx].gamma = NULL;
                param.null_skills[gidx].xi = NULL;
                param.null_skills[gidx].c = NULL;
                param.null_skills[gidx].p_O_param = 0.0;
                continue;
            }
            if( skill_group_map[k][g]==0)
                printf("ERROR! position [%d,%d] in skill_group_map in should have been 1\n",k,g);
            else if( skill_group_map[k][g]==1 ) { // insert new sequence and grab new data
                // link
                param.k_data[ k_countg[k] ] = &param.all_data[n_all_data]; // in linear array
                param.g_data[ g_countk[g] ] = &param.all_data[n_all_data]; // in linear array
                // param.k_g_data[k] and param.g_k_data[g] are already linked
//                param.k_g_data[k][ k_countg[k] ] = Calloc(struct data, 1); // grab
//                param.g_k_data[g][ g_countk[g] ] = param.k_g_data[k][ k_countg[k] ]; // relink
//                param.k_g_data[k][ k_countg[k] ]->ndat = 1; // init data << VV
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
                param.all_data[n_all_data].ndat = 1; // init data << VV
                param.all_data[n_all_data].k = k; // init k
                param.all_data[n_all_data].g = g; // init g
                param.all_data[n_all_data].cnt = 0;
                param.all_data[n_all_data].obs = NULL;
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
                    param.k_data[gidx]->ndat++;
                else
                    printf("ERROR! position of group %d in skill %d not found\n",g,k);
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
	k_countg = Calloc(NDAT, param.nK); // track current group in skill
	g_countk = Calloc(NDAT, param.nG); // track current skill in group
	for(t=0; t<param.N; t++) {
		g = param.dat_group->get(t);//[t];
		o = param.dat_obs->get(t);//[t];
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
                NCAT gidx = index_null_skill_group[g];
                param.null_skills[gidx].obs[ param.null_skills[gidx].cnt++ ] = o; // use .cnt as counter
                continue;
            }
            if( skill_group_map[k][g]<2)
                printf("ERROR! position [%d,%d] in skill_group_map should have been 2\n",k,g);
            else if( skill_group_map[k][g]==2 ) { // grab data and insert first dat point
                param.k_g_data[k][ k_countg[k] ]->obs = Calloc(NPAR, param.k_g_data[k][ k_countg[k] ]->ndat); // grab
                param.k_g_data[k][ k_countg[k] ]->obs[0] = o; // insert
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
                    param.k_g_data[k][ gidx ]->obs[pos] = o; // insert
                    param.k_g_data[k][ gidx ]->cnt++; // increase data counter
                }
                else
                    printf("ERROR! position of group %d in skill %d not found\n",g,k);
            }
        }
    }
	// recycle
	free(k_countg);
	free(g_countk);
    free(index_null_skill_group);
    if(param.cv_folds==0 && param.solver != BKT_GD_T) {
        if(param.multiskill==0) {
            delete param.dat_skill;//.clear();
            param.dat_skill = NULL;
        }
        else {
            delete param.dat_multiskill;
            param.dat_multiskill = NULL;
        }
        delete param.dat_group;
        delete param.dat_obs;
        param.dat_group = NULL;
        param.dat_obs = NULL;
    }
	free2DNCat(skill_group_map, param.nK);
    // reset `cnt'
    for(g=0; g<param.nG; g++) // for all groups
            for(k=0; k<param.g_numk[g]; k++) // for all skills in it
                param.g_k_data[g][k]->cnt = 0;
    for(NCAT x=0; x<param.n_null_skill_group; x++)
            param.null_skills[x].cnt = 0;
}

void destroy_input_data() {
	if(param.init_params != NULL) free(param.init_params);
	if(param.param_lo != NULL) free(param.param_lo);
	if(param.param_hi != NULL) free(param.param_hi);
	
    // data
	if(param.dat_obs != NULL) free(param.dat_obs);
	if(param.dat_group != NULL) free(param.dat_group);
	if(param.dat_skill != NULL) free(param.dat_skill);
	if(param.dat_multiskill != NULL) free(param.dat_multiskill);

    free(param.all_data); // ndat of them
    free(param.k_data); // ndat of them (reordered by k)
    free(param.g_data); // ndat of them (reordered by g)
    free(param.k_g_data); // nK of them
    free(param.g_k_data); // nG of them
//	NCAT k,g;
//	NDAT t;
//	for(k=0; k<param.nK; k++) {
//		for(g=0; g<param.k_numg[k]; g++) {
//			free(param.k_g_data[k][g]->obs); // only free data here
//			if( param.k_g_data[k][g]->c != NULL )
//				free(param.k_g_data[k][g]->c);  // only free data here
//			if( param.k_g_data[k][g]->alpha != NULL )
//				free2DNumber(param.k_g_data[k][g]->alpha,  param.k_g_data[k][g]->ndat);  // only free data here
//			if( param.k_g_data[k][g]->beta != NULL )
//				free2DNumber(param.k_g_data[k][g]->beta,  param.k_g_data[k][g]->ndat); // only free data here
//			if( param.k_g_data[k][g]->gamma != NULL )
//				free2DNumber(param.k_g_data[k][g]->gamma,  param.k_g_data[k][g]->ndat); // only free data here
//			if( param.k_g_data[k][g]->xi != NULL )
//				for(t=0;t<param.k_g_data[k][g]->ndat; t++)
//					free2DNumber(param.k_g_data[k][g]->xi[t],  param.nS); // only free data here
//			free(param.k_g_data[k][g]);
//		}
//		free(param.k_g_data[k]);
//	}
//	for(g=0; g<param.nG; g++)
//		free(param.g_k_data[g]);
    
	free(param.k_numg);
	free(param.g_numk);
    // null skills
    for(NCAT g=0;g<param.n_null_skill_group; g++)
        free(param.null_skills[g].obs);
    free(param.null_skills);
    // vocabularies
    delete param.map_group_fwd;
    delete param.map_group_bwd;
    delete param.map_step;
    delete param.map_skill_fwd;
    delete param.map_skill_bwd;
}

void read_predict_data(const char *filename) {
	FILE *fid = fopen(filename,"r");
	int number_columns = 0;
	max_line_length = 1024;
	char *col;
	
    
	// count lines and check for number of columns
	line = (char *)malloc(max_line_length);// Malloc(char,max_line_length);
	
	// grab memory and read all data
	param.dat_obs   = new StripedArray<NPAR>();//Malloc(NPAR, global_N);
	param.dat_group = new StripedArray<NCAT>();//Malloc(NCAT, global_N);
    if(param.multiskill==0)
        param.dat_skill = new StripedArray<NCAT>();//Malloc(NCAT, global_N);
    else
        param.dat_multiskill = new StripedArray< NCAT* >(true);
    param.map_group_fwd = new map<string,NCAT>();
    param.map_group_bwd = new map<NCAT,string>();
    param.map_skill_fwd = new map<string,NCAT>();
    param.map_skill_bwd = new map<NCAT,string>();
	string s_group, s_step, s_skill;
	map<string,NCAT>::iterator it;
	bool wrong_no_columns = false;
    param.N = 0;
    param.N_null = 0;
	while( readline(fid)!=NULL && !wrong_no_columns) {
//    while( NULL!=fgets(line,max_line_length,fid) && !wrong_no_columns) {
		number_columns = 0;
		// Observation
		col = strtok(line,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		NPAR obs = (NPAR)(atoi( col )-1);
		if(obs==NPAR_MAX) {
			fprintf(stderr,"Number of observtions exceeds allowed maximum of %d.\n",NPAR_MAX);
			exit(1);
		}
		param.dat_obs->add((NPAR)obs); // dat_obs[t] = (NPAR)obs;
		if( (obs >= 0) && ((param.nO-1) < obs) )
			param.nO = obs + (NPAR)1; // obs[t] + 1;
		// Group
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_group = string( col );
		it = param.map_group_fwd->find(s_group);
		if( it==param.map_group_fwd->end() ) { // not found
			if(param.map_group_fwd->size()==NCAT_MAX) {
				fprintf(stderr,"Number of unique groups exceeds allowed maximum of %d.\n",NCAT_MAX);
				exit(1);
			}
			NCAT newg = param.map_group_fwd->size();
			param.dat_group->add(newg); //[t] = param.map_group_fwd.size();
			param.map_group_fwd->insert(pair<string,NCAT>(s_group, newg));
			param.map_group_bwd->insert(pair<NCAT,string>(newg, s_group));
		}
		else
			param.dat_group->add(it->second); // [t] = it->second;
		
		
		// Step
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_step = string( col );
		//		it = param.map_step.find(s_step);
		//		if( it==param.map_step.end() ) { // not found
		//			if(param.map_step.size()==NCAT_MAX) {
		//				fprintf(stderr,"Number of unique steps exceeds allowed maximum of %d.\n",NCAT_MAX);
		//				exit(1);
		//			}
		//			dat_step[t] = param.map_step.size();
		//			param.map_step.insert(pair<string,NCAT>(s_step, param.map_step.size()));
		//		}
		//		else
		//			dat_step[t] = it->second;
		
		// Skill
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_skill = string( col );
		if( (s_skill.empty() || ( s_skill.size()==1 && (s_skill[0]=='.' || s_skill[0]==' ') ) ) ) { // null skill
            param.N_null++;
            if(param.multiskill == 0) {
                param.dat_skill->add(-1); // [t] = -1;
            }
            else {
                NCAT* a_skills = Malloc(NCAT, 2);
                a_skills[0] = 1; // count
                a_skills[1] = -1; // value
                param.dat_multiskill->add(a_skills);
            }
		}
		else {
            // multiskill
            if(param.multiskill != 0) {
                list<NCAT> a_skills;//
                char* a_kc;
                a_kc  = strtok(col,"~\n\r");
                string s_kc;
                while(a_kc != NULL) {
                    s_kc = string(a_kc);
                    // adding vvvv
                    it = param.map_skill_fwd->find(s_kc);
                    if( it==param.map_skill_fwd->end() ) { // not found
                        if(param.map_skill_fwd->size()==NCAT_MAX) {
                            fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                            exit(1);
                        }
                        a_skills.insert(a_skills.end(), param.map_skill_fwd->size()); //dat_skill->add(param.map_skill_fwd->size());
                        param.map_skill_fwd->insert(pair<string,NCAT>(s_kc, param.map_skill_fwd->size()));
                        param.map_skill_bwd->insert(pair<NCAT,string>(param.map_skill_bwd->size(),s_kc));
                    }
                    else
                        a_skills.insert(a_skills.end(), it->second); //dat_skill->add(it->second); //[t] = it->second;
                    // adding ^^^^
                    a_kc  = strtok(NULL,"~\n\r");
                }
                NCAT *b_skills = Malloc(NCAT, a_skills.size()+1);
                b_skills[0] = a_skills.size();
                int count = 0;
                for(list<NCAT>::iterator it=a_skills.begin(); it!=a_skills.end(); it++)
                    b_skills[++count] = *it;
                param.dat_multiskill->add(b_skills);
                // multi skill
            } else {
                // single skill
                it = param.map_skill_fwd->find(s_skill);
                if( it==param.map_skill_fwd->end() ) { // not found
                    if(param.map_skill_fwd->size()==NCAT_MAX) {
                        fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                        exit(1);
                    }
                    param.dat_skill->add(param.map_skill_fwd->size()); //[t] = param.map_skill_fwd.size();
                    param.map_skill_fwd->insert(pair<string,NCAT>(s_skill, param.map_skill_fwd->size()));
                    param.map_skill_bwd->insert(pair<NCAT,string>(param.map_skill_bwd->size(),s_skill));
                }
                else
                    param.dat_skill->add(it->second); //[t] = it->second;
            } // single skill
		}
		// count lines
		param.N++;	// increase line count
//        fprintf(stdout,"Line %d\n",param.N);
	}// reading loop
	if(wrong_no_columns) {
		fprintf(stderr,"Wrong number of columns in line %u. Expected %d, found %d\n",param.N+1,COLUMNS, number_columns);
        if(param.multiskill==0) {
            delete param.dat_skill;//.clear();
            param.dat_skill = NULL;
        }
        else {
            delete param.dat_multiskill;
            param.dat_multiskill = NULL;
        }
        delete param.dat_group;
        delete param.dat_obs;
        param.dat_group = NULL;
        param.dat_obs = NULL;
		free(line);
		fclose(fid);
		exit(1);
	}
	fclose(fid);
	free(line);
//    if(param.cv_folds>0) {
//        delete param.map_group_fwd;
//        delete param.map_group_bwd;
//        delete param.map_step;
//        delete param.map_skill_fwd;
//        delete param.map_skill_bwd;
//    }
}

void cross_validate(NUMBER* metrics, const char *filename) {
    NUMBER rmse = 0.0;
    NUMBER rmse_no_null = 0.0;
    NPAR f;
    NCAT g,k;
    FILE *fid; // file for storing prediction should that be necessary
    if(param.predictions>0) {
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
        switch(param.solver)
        {
            case BKT_CGD: // Conjugate Gradient Descent
            case BKT_GD: // Gradient Descent
            case BKT_BW: // Expectation Maximization (Baum-Welch)
            case BKT_GD_BW: // Gradient Descent then Expectation Maximization (Baum-Welch)
            case BKT_BW_GD: // Expectation Maximization (Baum-Welch) then Gradient Descent
            case BKT_GD_G: // Gradient Descent by group
                hmms[f] = new HMMProblem(&param);
                break;
//            case BKT_GD_PIg: // Gradient Descent: PI by group, A,B by skill
//                hmms[f] = new HMMProblemPiG(&param);
//                break;
//            case BKT_GD_PIgk: // Gradient Descent, pLo=f(K,G), other by K
//                hmms[f] = new HMMProblemPiGK(&param);
//                break;
//            case BKT_GD_APIgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
//                hmms[f] = new HMMProblemPiAGK(&param);
//                break;
//            case BKT_GD_Agk: // Gradient Descent, pT=f(K,G), other by K
//                hmms[f] = new HMMProblemAGK(&param);
//                break;
//            case BKT_GD_T: // Gradient Descent with Transfer
//                hmms[f] = new HMMProblemKT(&param);
//                break;
        }
        // block respective data
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
        hmms[f]->fit();
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
//    NUMBER* result = Calloc(NUMBER, param.N); // just one observation, assuming it is first one - correct (index 0)
	for(g=0; g<param.nG; g++)
		for(k=0; k<param.nK; k++) {
            f = folds[g];
            dt.k = k;
            dt.g = g;
			for(i=0; i<param.nO; i++)
                group_skill_map[g][k][i] = hmms[f]->getPI(&dt,i);//PI[i];
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
            prob = safe0num(hmms[f]->getNullSkillObs(param.cv_target_obs));
            ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
            if(param.predictions>0) // write predictions file if it was opened
                for(m=0; m<param.nO; m++)
                    fprintf(fid,"%10.8f%s",hmms[f]->getNullSkillObs(m),(m<(param.nO-1))?"\t":"\n");
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
                fprintf(fid,"%10.8f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
        rmse += pow(isTarget-local_pred[param.cv_target_obs],2);
        rmse_no_null += pow(isTarget-local_pred[param.cv_target_obs],2);
        prob = safe01num(local_pred[param.metrics_target_obs]);
        ll -= safelog(  prob)*   isTarget  +  safelog(1-prob)*(1-isTarget);
	} // for all data
    rmse = sqrt( rmse / param.N );
    rmse_no_null = sqrt( rmse_no_null / (param.N - param.N_null) );
    // delete data
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
    metrics[0] = 2*(n_par) + 2*ll;
    metrics[1] = n_par*safelog(param.N) + 2*ll;
    metrics[2] = rmse;
    metrics[3] = rmse_no_null;
}