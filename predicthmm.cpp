/*
 *  predicthmm.cpp
 *  HMM
 *
 *  Created by Mikhail Yudelson on 6/6/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
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
//#include "HMMProblemPiG.h"
#include "HMMProblemPiGK.h"
#include "HMMProblemAGK.h"
#include "HMMProblemPiAGK.h"
//#include "HMMProblemKT.h"
using namespace std;

#define COLUMNS 4

struct param param;
static char *line = NULL;
static int max_line_length;
static NDAT global_predict_N;
HMMProblem *hmm;
NUMBER* metrics;
map<string,NCAT> data_map_group_fwd;
map<NCAT,string> data_map_group_bwd;
map<string,NCAT> map_step;
map<string,NCAT> model_map_skill_fwd;
map<NCAT,string> model_map_skill_bwd;
void exit_with_help();
void parse_arguments(int argc, char **argv, char *input_file_name, char *model_file_name, char *predict_file_name);
void read_predict_data(const char *filename);
void read_model(const char *filename);
//void predict(/*NUMBER **result, */const char *input_file, const char *predict_file);
//void write_prediction(const char *filename, NUMBER **result);
static char* readline(FILE *fid);

int main (int argc, char ** argv) {
	clock_t tm0 = clock();
	printf("predicthmm starting...\n");
	set_param_defaults(&param);
	
	char input_file[1024];
	char model_file[1024];
	char predict_file[1024];
	
	parse_arguments(argc, argv, input_file, model_file, predict_file);
    param.predictions = 1; // force it on, since we, you know, predictinng :)
	read_model(model_file);
    
	read_predict_data(input_file);
	
	if(param.quiet == 0)
		printf("input read, nO=%d, nG=%d, nK=%d\n",param.nO, param.nG, param.nK);
	
    //	// default params
    //	defPI = init1DNumber(param.nS); defPI[0] = 0.5; defPI[1] = 0.5;
    //	defA = init2DNumber(param.nS,param.nS); defA[0][0] = 1.0; defA[0][1] = 0.0; defA[1][0] = 0.4; defA[1][1] = 0.6;
    //	defB = init2DNumber(param.nS,param.nO); defB[0][0] = 0.8; defB[0][1] = 0.2; defB[1][0] = 0.2; defB[1][1] = 0.8;
    //
	clock_t tm = clock();
    if(param.metrics>0 || param.predictions>0) {
        metrics = Calloc(NUMBER, 5); // AIC, BIC, RMSE
    }
    //	predict(input_file, predict_file);
    hmm->predict(metrics, predict_file, param.dat_obs, param.dat_group, param.dat_skill, param.dat_multiskill, true/*only unlabelled*/);
	if(param.quiet == 0)
		printf("predicting is done in %8.6f seconds\n",(NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
    if( param.metrics>0 ) {
        printf("trained model LL=%15.7f, AIC=%8.6f, BIC=%8.6f, RMSE=%8.6f (%8.6f)\n",metrics[0], metrics[1], metrics[2], metrics[3], metrics[4]);
    }
    free(metrics);
    
	destroy_input_data(&param);
	
    delete hmm;
	if(param.quiet == 0)
		printf("overall time running is %8.6f seconds\n",(NUMBER)(clock()-tm0)/CLOCKS_PER_SEC);
    return 0;
}

void exit_with_help() {
	printf(
		   "Usage: trainhmm [options] input_file [[output_file] predicted_response_file]\n"
		   "options:\n"
		   "(-s) : structure.solver[.solver setting], structures: 1-by skill, 2-by user,\n"
           "       3-Pi by skill, A,B-by user, 4-Pi by skill and user, A,B-by user,\n"
           "       5-A by skill and user, Pi,B by skill, 6-Pi,A by skill and user, B by skill;\n"
           "       solvers: 1-Baum-Welch, 2-Gradient Descent, 3-Conjugate Gradient Descent;\n"
           "       Conjugate Gradient Descent has 3 setings: 1-Polak-Ribiere, 2-Fletcher–Reeves,\n"
           "       3-Hestenes-Stiefel.\n"
           "       For example '-s 1.3.1' would be by skill structure (classical) with Conjugate\n"
           "       Gradient Descent and Hestenes-Stiefel formula, '-s 2.1' would be by student structure\n"
           "       fit using Baum-Welch method.\n"
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

void parse_arguments(int argc, char **argv, char *input_file_name, char *model_file_name, char *predict_file_name) {
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
	int i;
    char * ch;
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break; // end of options stop parsing
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 'q':
				param.quiet = atoi(argv[i]);
				if(param.quiet!=0 && param.quiet!=1) {
					fprintf(stderr,"ERROR! Quiet param should be 0 or 1\n");
					exit_with_help();
				}
				break;
                //			case 'n':
                //				param.nS = (NPAR)atoi(argv[i]);
                //				if(param.nS<2) {
                //					fprintf(stderr,"ERROR! Number of hidden states should be at least 2\n");
                //					exit_with_help();
                //				}
                //				//fprintf(stdout, "fit single skill=%d\n",param.quiet);
                //				break;
			case 'm':
                param.metrics = atoi( strtok(argv[i],";\t\n\r"));
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
            case  'd':
				param.multiskill = argv[i][0]; // just grab first character (later, maybe several)
                break;
			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
    if(param.cv_target_obs>0 && param.metrics>0) { // correct for 0-start coding
        fprintf(stderr,"values for -v and -m cannot be both non-zeros\n");
        exit_with_help();
    }
	
	// next argument should be input (training) file name
	if(i>=argc) // if not
		exit_with_help(); // leave
	
	strcpy(input_file_name, argv[i++]); // copy and advance
	
	if(i>=argc) { // no model file name specified
        fprintf(stderr,"no model file specified\n");
		exit_with_help(); // leave
	}
	else {
		strcpy(model_file_name,argv[i++]); // copy and advance
		if(i>=argc) // no prediction file name specified
			strcpy(predict_file_name,"predict_hmm.txt"); // the predict file too
		else
			strcpy(predict_file_name,argv[i]);
	}
}

void read_predict_data(const char *filename) {
	FILE *fid = fopen(filename,"r");
    //	global_N = 0;
    //	global_predict_N = 0;
	int number_columns = 0;
	max_line_length = 1024;
	char *col;
	if(fid == NULL)
	{
		fprintf(stderr,"Can't open input file %s\n",filename);
		exit(1);
	}
    
	// grab memory and read all data
	line = (char *)malloc(max_line_length);// Malloc(char,max_line_length);
	param.dat_obs   = new StripedArray<NPAR>();//Malloc(NPAR, global_N);
	param.dat_group = new StripedArray<NCAT>();//Malloc(NCAT, global_N);
    if(param.multiskill==0)
        param.dat_skill = new StripedArray<NCAT>();//Malloc(NCAT, global_N);
    else
        param.dat_multiskill = new StripedArray< NCAT* >(true);
	
	string s_obs, s_group, s_step, s_skill;
	map<string,NCAT>::iterator it_k;
	map<string,NCAT>::iterator it_g;
    NPAR obs = 0;
    
	bool wrong_no_columns = false;
    param.N = 0;
    param.N_null = 0;
	while( readline(fid)!=NULL && !wrong_no_columns) {
		number_columns = 0;
        
		// Observation
		col = strtok(line,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
        s_obs = string( col );
		if( (s_obs.empty() || ( s_obs.size()==1 && (s_obs[0]=='.' || s_obs[0]==' '/**/) ) ) ) {
			// is emply as labelled!
			obs = -1;//dat_obs[t] = -1;
			global_predict_N++;
		} else {
			obs = atoi( s_obs.c_str() )-1;
			if(obs==NPAR_MAX) {
				fprintf(stderr,"Number of observtions exceeds allowed maximum of %d.\n",NPAR_MAX);
				exit(1);
			}
			if(obs>(param.nO-1)) {
				fprintf(stderr,"Number of observtions exceeds preset of %d.\n",param.nO);
				exit(1);
			}
		}
        param.dat_obs->add(obs);
        
		// Group
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_group = string( col );
		it_g = param.map_group_fwd->find(s_group);
		if( it_g==param.map_group_fwd->end() ) { // not found
            //			param.dat_group->add(NCAT(-1));
			if(param.map_group_fwd->size()==NCAT_MAX) {
				fprintf(stderr,"Number of unique groups exceeds allowed maximum of %d.\n",NCAT_MAX);
				exit(1);
			}
			NCAT newg = param.map_group_fwd->size();
			param.dat_group->add(newg); //[t] = param.map_group_fwd.size();
			param.map_group_fwd->insert(pair<string,NCAT>(s_group, newg));
			param.map_group_bwd->insert(pair<NCAT,string>(newg, s_group));
		}
		else                                     // found
			param.dat_group->add(it_g->second);
        
		// Step
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_step = string( col );
        
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
                param.dat_skill->add((NCAT)-1); // [t] = -1;
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
                    it_k = param.map_skill_fwd->find(s_kc);
                    if( it_k==param.map_skill_fwd->end() ){         // not found
                        a_skills.insert(a_skills.end(), (NCAT)-1);
                        param.N_null++;
                    }
                    else                                            // found
                        a_skills.insert(a_skills.end(), it_k->second);
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
                it_k = param.map_skill_fwd->find(s_skill);
                if( it_k==param.map_skill_fwd->end() ){ // not found
                    param.dat_skill->add( (NCAT)-1 );
                    param.N_null++;
                }
                else                                    // found
                    param.dat_skill->add(it_k->second);
            } // single skill
		}
		// count lines
		param.N++;	// increase line count
        //        fprintf(stdout,"Line %d\n",param.N);
	}// reading loop
    
	free(line);
	fclose(fid);
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

void read_model(const char *filename) {
	FILE *fid = fopen(filename,"r");
	if(fid == NULL)
	{
		fprintf(stderr,"Can't read model file %s\n",filename);
		exit(1);
	}
	max_line_length = 1024;
	line = Malloc(char,max_line_length);
	NDAT line_no = 0;
    //
    // read solver info
    //
    readSolverInfo(fid, &param, &line_no);
    
    //
    // create hmm Object
    //
    switch(param.structure)
    {
        case STRUCTURE_SKILL: // Conjugate Gradient Descent
        case STRUCTURE_GROUP: // Conjugate Gradient Descent
            hmm = new HMMProblem(&param);
            break;
            //            case STRUCTURE_PIg: // Gradient Descent: PI by group, A,B by skill
            //                hmm = new HMMProblemPiG(&param);
            //                break;
        case STRUCTURE_PIgk: // Gradient Descent, pLo=f(K,G), other by K
            hmm = new HMMProblemPiGK(&param);
            break;
        case STRUCTURE_PIAgk: // Gradient Descent, pLo=f(K,G), pT=f(K,G), other by K
            hmm = new HMMProblemPiAGK(&param);
            break;
        case STRUCTURE_Agk: // Gradient Descent, pT=f(K,G), other by K
            hmm = new HMMProblemAGK(&param);
            break;
            //            case BKT_GD_T: // Gradient Descent with Transfer
            //                hmm = new HMMProblemKT(&param);
            //                break;
    }
    //
    // read model
    //
    hmm->readModel(fid, &line_no);
    
	
    //	k=0;
    //	map<NCAT,string>::iterator it;
    //	for(k=0; k<param.nK; k++) {
    //		it	= model_map_skill_bwd.find(k);
    //		printf("%d %d %s \n", k, it->first, it->second.c_str());
    //	}
	
	fclose(fid);
	free(line);
}

//void predict(/*NUMBER **result, */const char *input_file, const char *predict_file) {
//	FILE *fid = fopen(input_file,"r");
//	FILE *fidP = fopen(predict_file,"w");
//	if(fid == NULL)
//	{
//		fprintf(stderr,"Can't open input file %s\n",input_file);
//		exit(1);
//	}
//	if(fidP == NULL)
//	{
//		fprintf(stderr,"Can't open prediction output file %s\n",input_file);
//		exit(1);
//	}
//
//	NDAT t;
//	NCAT g, k;
//	NPAR i, j, m, o;
//	NUMBER *local_pred = init1DNumber(param.nO); // local prediction
//	NUMBER pLe[param.nS];// p(L|evidence);
//	NUMBER pLe_denom; // p(L|evidence) denominator
//	NUMBER ***group_skill_map = init3DNumber(param.nG, param.nK, param.nS);
//	// initialize
//
//	for(g=0; g<param.nG; g++)
//		for(k=0; k<param.nK; k++) {
//			for(i=0; i<param.nO; i++)
//					group_skill_map[g][k][i] = PI[k][i];
//
//		}
//	NDAT predict_idx = 0;
//	string s_obs, s_group, s_step, s_skill;
//	map<string,NCAT>::iterator it;
//
//	max_line_length = 1024;
//	line = Malloc(char,max_line_length);
//	for(t=0; t<global_N; t++) {
//		readline(fid);// read file not
//		s_obs = strtok(line,"\t\n\r");
//		s_group = string( strtok(NULL,"\t\n\r") );
//		s_step = string( strtok(NULL,"\t\n\r") );
//		s_skill = string( strtok(NULL,"\t\n\r") );
//
//		if( (s_obs.empty() || ( s_obs.size()==1 && (s_obs[0]=='.' || s_obs[0]==' '/**/) ) ) )
//			o = -1;
//		else
//			o = atoi(s_obs.c_str())-1;
//
//		it = model_map_skill_fwd.find(s_skill);
//		if( it==model_map_skill_fwd.end() )
//			k = -1;
//		else
//			k = it->second;
//		it = data_map_group_fwd.find(s_group);
//		if( it==data_map_group_fwd.end() ) {
//			fprintf(stderr,"Group not found in line %d while predicting.\n",t);
//			exit(1);
//		}
//		else
//			g = it->second;
//
//
//		// produce prediction and copy to result
//		if(k<0) { // if no skill label
//			//for(m=0; m<param.nO; m++)
//			//	result[t][m] = param.null_obs_ratio[m];
//			if(o==-1) {// if output
//				for(m=0; m<param.nO; m++)
//					fprintf(fidP,"%10.8f%s",param.null_obs_ratio[m],(m<(param.nO-1))?"\t":"\n");
//				predict_idx++;
//			}
//			continue;
//		}
//		NUMBER **a_A = A[k];
//		NUMBER **a_B = B[k];
//		if(o==-1) { // if we need to predict
//			for(m=0; m<param.nO; m++) local_pred[m] = 0.0;
//			for(m=0; m<param.nO; m++)
//				for(i=0; i<param.nS; i++)
//					local_pred[m] += group_skill_map[g][k][i] * a_B[i][m];
//			//cpy1DNumber(local_pred, result[predict_idx], param.nO);
//			for(m=0; m<param.nO; m++)
//				fprintf(fidP,"%10.8f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
//			predict_idx++;
//		}
//
//		// update p(L)
//		if(o>=0) { // known observation
//			pLe_denom = 0.0;
//			// 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
//			for(i=0; i<param.nS; i++) pLe_denom += group_skill_map[g][k][i] * a_B[i][o];
//			for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i] * a_B[i][o] / safe0num(pLe_denom);
//			// 2. L = (pLe'*A)';
//			for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0;
//			for(j=0; j<param.nS; j++)
//				for(i=0; i<param.nS; i++)
//					group_skill_map[g][k][j] += pLe[i] * a_A[i][j];
//		} else { // unknown observation
//			// 2. L = (pL'*A)';
//			for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i]; // copy first;
//			for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0; // erase old value
//			for(j=0; j<param.nS; j++)
//				for(i=0; i<param.nS; i++)
//					group_skill_map[g][k][j] += pLe[i] * a_A[i][j];
//		}
//	} // for all data
//	free(line);
//	free(local_pred);
//	fclose(fid);
//	fclose(fidP);
//}
