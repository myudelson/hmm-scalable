/*
 *  predicthmm.cpp
 *  HMM
 *
 *  Created by Mikhail Yudelson on 6/6/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <map>
#include "utils.h"
using namespace std;

#define COLUMNS 4

struct param param;
static char *line = NULL;
static int max_line_length;
static NDAT global_N;
static NDAT global_predict_N;
//static NPAR *dat_obs;
//static NCAT *dat_group;
//static NCAT *dat_skill;
map<string,NCAT> data_map_group_fwd;
map<NCAT,string> data_map_group_bwd;
map<string,NCAT> map_step;
map<string,NCAT> model_map_skill_fwd;
map<NCAT,string> model_map_skill_bwd;
NUMBER *defPI;
NUMBER **defA;
NUMBER **defB;

NUMBER **PI;
NUMBER ***A;
NUMBER ***B;

void exit_with_help();
void parse_arguments(int argc, char **argv, char *input_file_name, char *model_file_name, char *predict_file_name);
void read_predict_data(const char *filename);
void read_model(const char *filename);
void predict(/*NUMBER **result, */const char *input_file, const char *predict_file);
void write_prediction(const char *filename, NUMBER **result);
bool readline(FILE *fid);

int main (int argc, char ** argv) {
	clock_t tm0 = clock();
	printf("predicthmm starting...\n");
	set_param_defaults(&param);
	
	char input_file[1024];
	char model_file[1024];
	char predict_file[1024];
	
	parse_arguments(argc, argv, input_file, model_file, predict_file);
	read_model(model_file);
	read_predict_data(input_file);
	
	if(param.quiet == 0)
		printf("input read, nO=%d, nG=%d, nK=%d\n",param.nO, param.nG, param.nK);
	
	// default params
	defPI = init1DNumber(param.nS); defPI[0] = 0.5; defPI[1] = 0.5;
	defA = init2DNumber(param.nS,param.nS); defA[0][0] = 1.0; defA[0][1] = 0.0; defA[1][0] = 0.4; defA[1][1] = 0.6;
	defB = init2DNumber(param.nS,param.nO); defB[0][0] = 0.8; defB[0][1] = 0.2; defB[1][0] = 0.2; defB[1][1] = 0.8;
	
	clock_t tm = clock();
	predict(input_file, predict_file);
	if(param.quiet == 0)
		printf("predicting is done in %8.6f seconds\n",(NUMBER)(clock()-tm)/CLOCKS_PER_SEC);
	
	// free data
	free2DNumber(PI, param.nK);
	free3DNumber(A, param.nK, param.nS);
	free3DNumber(B, param.nK, param.nS);
	free(defPI);
	free2DNumber(defA, param.nS);
	free2DNumber(defB, param.nS);
	//free(dat_skill);
	//free(dat_group);
	//free(dat_obs);
	data_map_group_fwd.clear();
	data_map_group_bwd.clear();
	map_step.clear();
	model_map_skill_bwd.clear();
	model_map_skill_fwd.clear();
	
	if(param.quiet == 0)
		printf("overall time running is %8.6f seconds\n",(NUMBER)(clock()-tm0)/CLOCKS_PER_SEC);
    return 0;
}

void exit_with_help(){
	printf(
		   "Usage: trainhmm [options] input_file [[output_file] predicted_response_file]\n"
		   "options:\n"
		   "-q : quiet mode, withour output, 0-no (default), or 1-yes\n"
		   "-s : number of hidden states, should be 2 or more (default 2)\n"
		   );
	exit(1);
}
 
void parse_arguments(int argc, char **argv, char *input_file_name, char *model_file_name, char *predict_file_name) {	
	// parse command line options, starting from 1 (0 is path to executable)
	// go in pairs, looking at whether first in pair starts with '-', if not, stop parsing arguments
	int i;
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break; // end of options stop parsing
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
			case 'q':
				param.quiet = atoi(argv[i]);
				//fprintf(stdout, "quiet=%d\n",param.quiet);
				break;
			case 'n':
				param.nS = atoi(argv[i]);
				if(param.nS<2) {
					printf("ERROR! Number of hidden states should be at least 2\n");
					exit(1);
				}
				//fprintf(stdout, "fit single skill=%d\n",param.quiet);
				break;
			default:
				fprintf(stderr,"unknown option: -%c\n", argv[i-1][1]);
				exit_with_help();
				break;
		}
	}
	
	// next argument should be input file name
	if(i>=argc) // if not
		exit_with_help(); // leave 
	
	strcpy(input_file_name, argv[i++]); // copy and advance
	
	if(i>=argc) { // no output file name specified
		fprintf(stderr,"Error! no model file specified\n");
		exit_with_help();
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
	global_N = 0;
	global_predict_N = 0;
	int number_columns = 0;
	max_line_length = 1024;
	char *col;
	
	if(fid == NULL)
	{
		fprintf(stderr,"Can't open input file %s\n",filename);
		exit(1);
	}
	
	// count lines and check for number of columns
	line = Malloc(char,max_line_length);
	
	while( readline(fid)/*, line, max_line_length)!=NULL*/) {
		number_columns = 0;
		col = strtok(line,"\t\n\r");
		// columns
		while(col != NULL)
		{
			number_columns++;
			col = strtok(NULL,"\t\n\r");
		}
		if(number_columns != COLUMNS) {
			fprintf(stderr,"Wrong number of columns in line %d. Expected %d, found %d\n",global_N+1,COLUMNS, number_columns);
			exit(1);
		}
		global_N++;	// increase line count
	}
	rewind(fid);
	
	// grab memory and read all data
	//dat_obs   = Malloc(NPAR, global_N);
	//dat_group = Malloc(NCAT, global_N);
	//NCAT *dat_step  = Malloc(NCAT, global_N);
	//dat_skill = Malloc(NCAT, global_N);
	
	//	1. read raw data 
	//		create
	//			dat_group[N] - all group data
	//			dat_skill[N] - all skill data
	//			dat_obs[N]   - all observation data
	//			nG - number of unique groups                                    RETAIN
	//			nK - number of unique skills                                    RETAIN
	//			nO - number of unique observations                              RETAIN
	//		purge
	//			N/A
	string s_obs, s_group, s_step, s_skill; 
	map<string,NCAT>::iterator it;
	for(NDAT t=0; t<global_N; t++) {
		readline(fid/*, line, max_line_length*/);
		s_obs = strtok(line,"\t\n\r");
		if( (s_obs.empty() || ( s_obs.size()==1 && (s_obs[0]=='.' || s_obs[0]==' '/**/) ) ) ) {
			// is emply as labelled!
			//dat_obs[t] = -1;
			global_predict_N++;
		} else {
			int obs = atoi( s_obs.c_str() )-1;
			if(obs==NPAR_MAX) {
				fprintf(stderr,"Number of observtions exceeds allowed maximum of %d.\n",NPAR_MAX);
				exit(1);
			}
			if(obs>(param.nO-1)) {
				fprintf(stderr,"Number of observtions exceeds preset of %d.\n",param.nO);
				exit(1);
			}
			//dat_obs[t] = (NPAR)obs;
		}
		
		s_group = string( strtok(NULL,"\t\n\r") );
		it = data_map_group_fwd.find(s_group);
		if( it==data_map_group_fwd.end() ) { // not found
			if(data_map_group_fwd.size()==NCAT_MAX) {
				fprintf(stderr,"Number of unique groups exceeds allowed maximum of %d.\n",NCAT_MAX);
				exit(1);
			}
			//dat_group[t] = data_map_group_fwd.size();
			data_map_group_fwd.insert(pair<string,NCAT>(s_group, data_map_group_fwd.size()));
			data_map_group_bwd.insert(pair<NCAT,string>(data_map_group_bwd.size(),s_group));
		}
	}//main reading loop
	param.nG = data_map_group_fwd.size();
	free(line);
	fclose(fid);
}

bool readline(FILE *fid) {
	int length = 0;
	
	if(fgets(line,max_line_length,fid) == NULL)
		return false;
	
	while(strrchr(line,'\n') == NULL)
	{
		max_line_length *= 2;
		line = (char *) realloc(line,max_line_length);
		length = (int) strlen(line);
		if(fgets(line+length,max_line_length-length,fid) == NULL)
			break;
	}
	return true;
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
	char *col;
	NPAR num_col = 0;
	NDAT line_no = 0;
	//
	// count null skill ratios
	//
	if(readline(fid)/*, line, max_line_length)!=NULL*/) {
		line_no++;
		col = strtok(line,"\t\n\r");
		// columns
		while(col != NULL)
		{
			num_col++;
			col = strtok(NULL,"\t\n\r");
		}
		param.nO = num_col-1; // first instance is text
	}
	else {
		fprintf(stderr,"Input file is empty\n");
		exit(1);
	}
	NPAR position = 0; // 0 - skill name, 1 - PI, 2 - A, 3 - B; 
	NCAT num_skill = 0;
	//
	// count skills
	//
	while( readline(fid)/*, line, max_line_length)!=NULL*/) {
		line_no++;
		num_col = 0;
		col = strtok(line,"\t\n\r"); // ID of skill or param text
		// columns
		while(col != NULL)
		{
			num_col++;
			col = strtok(NULL,"\t\n\r");
		}
		switch (position) {
			case 0: // skill name
				if(num_col != 2) {
					fprintf(stderr,"Error in model file format, line %d\n",line_no);
					exit(1);
				}
				break;
			case 1: // PI
				if(num_col != (param.nS+1)) {
					fprintf(stderr,"Error in model file format, line %d\n",line_no);
					exit(1);
				}
				break;
			case 2: // A
				if(num_col != (param.nS*param.nS+1)) {
					fprintf(stderr,"Error in model file format, line %d\n",line_no);
					exit(1);
				}
				break;
			case 3: // B
				if(num_col != (param.nS*param.nO+1)) {
					fprintf(stderr,"Error in model file format, line %d\n",line_no);
					exit(1);
				}
				else
					num_skill++;
				break;
		}
		position = (position+1) % 4;
	}
	if(position!=0) {
		fprintf(stderr,"Model is incomplete%d\n",line_no);
		exit(1);
	}
	rewind(fid);
	param.nK = num_skill;
	NPAR i,j,m;
	string s;
	//
	// read null skill ratios
	//
	param.null_obs_ratio = Calloc(NUMBER, param.nO);
	readline(fid/*, line, max_line_length*/);
	col = strtok(line,"\t\n\r");
	for(i=0; i<param.nO; i++) {
		param.null_obs_ratio[i] = atof( strtok(NULL,"\t\n\r") );
	}
	//
	// read skills
	//
	PI = init2DNumber(param.nK, param.nS);
	A = init3DNumber(param.nK, param.nS, param.nS);
	B = init3DNumber(param.nK, param.nS, param.nO);	
	NCAT k = 0;
	for(k=0; k<param.nK; k++) {
		// read label
		readline(fid/*, line, max_line_length*/);
		strtok(line,"\t\n\r"); // ID
		s = string( strtok(NULL,"\t\n\r") );
		model_map_skill_fwd.insert(pair<string,NCAT>(s, model_map_skill_fwd.size()));
		model_map_skill_bwd.insert(pair<NCAT,string>(model_map_skill_bwd.size(), s));
		// read PI
		readline(fid/*, line, max_line_length*/);
		s = string( strtok(line,"\t\n\r") ); // "PI"
		for(i=0; i<param.nS; i++) {
			s = string( strtok(NULL,"\t\n\r") );
			PI[k][i] = atof( s.c_str() );
		}
		// read A
		readline(fid/*, line, max_line_length*/);
		s = string( strtok(line,"\t\n\r") ); // "A"
		for(i=0; i<param.nS; i++)
			for(j=0; j<param.nS; j++) {
				s = string( strtok(NULL,"\t\n\r") );
				A[k][i][j] = atof( s.c_str() );
			}
		// read B
		readline(fid/*, line, max_line_length*/);
		col = strtok(line,"\t\n\r"); // "B"
		for(i=0; i<param.nS; i++)
			for(m=0; m<param.nO; m++) {
				s = string( strtok(NULL,"\t\n\r") );
				B[k][i][m] = atof( s.c_str() );
			}
	}
	
//	k=0;
//	map<NCAT,string>::iterator it;
//	for(k=0; k<param.nK; k++) {
//		it	= model_map_skill_bwd.find(k);
//		printf("%d %d %s \n", k, it->first, it->second.c_str());
//	}
	
	fclose(fid);
	free(line);
}

void predict(/*NUMBER **result, */const char *input_file, const char *predict_file) {
	FILE *fid = fopen(input_file,"r");
	FILE *fidP = fopen(predict_file,"w");
	if(fid == NULL)
	{
		fprintf(stderr,"Can't open input file %s\n",input_file);
		exit(1);
	}
	if(fidP == NULL)
	{
		fprintf(stderr,"Can't open prediction output file %s\n",input_file);
		exit(1);
	}
	
	NDAT t;
	NCAT g, k;
	NPAR i, j, m, o;
	NUMBER *local_pred = init1DNumber(param.nO); // local prediction
	NUMBER pLe[param.nS];// p(L|evidence);
	NUMBER pLe_denom; // p(L|evidence) denominator
	NUMBER ***group_skill_map = init3DNumber(param.nG, param.nK, param.nS);
	// initialize
	
	for(g=0; g<param.nG; g++)
		for(k=0; k<param.nK; k++) {
			for(i=0; i<param.nO; i++)
					group_skill_map[g][k][i] = PI[k][i];
			
		}
	NDAT predict_idx = 0;
	string s_obs, s_group, s_step, s_skill; 
	map<string,NCAT>::iterator it;

	max_line_length = 1024;
	line = Malloc(char,max_line_length);
	for(t=0; t<global_N; t++) {
		readline(fid);// read file not
		s_obs = strtok(line,"\t\n\r");
		s_group = string( strtok(NULL,"\t\n\r") );
		s_step = string( strtok(NULL,"\t\n\r") );
		s_skill = string( strtok(NULL,"\t\n\r") );
		
		if( (s_obs.empty() || ( s_obs.size()==1 && (s_obs[0]=='.' || s_obs[0]==' '/**/) ) ) )
			o = -1;
		else
			o = atoi(s_obs.c_str())-1;
		
		it = model_map_skill_fwd.find(s_skill);
		if( it==model_map_skill_fwd.end() )
			k = -1;
		else
			k = it->second;
		it = data_map_group_fwd.find(s_group);
		if( it==data_map_group_fwd.end() ) {
			fprintf(stderr,"Group not found in line %d while predicting.\n",t);
			exit(1);
		}
		else
			g = it->second;

		
		// produce prediction and copy to result
		if(k<0) { // if no skill label
			//for(m=0; m<param.nO; m++)
			//	result[t][m] = param.null_obs_ratio[m];
			if(o==-1) {// if output 
				for(m=0; m<param.nO; m++)
					fprintf(fidP,"%10.8f%s",param.null_obs_ratio[m],(m<(param.nO-1))?"\t":"\n");
				predict_idx++;
			}
			continue;
		}
		NUMBER **a_A = A[k];
		NUMBER **a_B = B[k];
		if(o==-1) { // if we need to predict
			for(m=0; m<param.nO; m++) local_pred[m] = 0.0;
			for(m=0; m<param.nO; m++)
				for(i=0; i<param.nS; i++)
					local_pred[m] += group_skill_map[g][k][i] * a_B[i][m];
			//cpy1DNumber(local_pred, result[predict_idx], param.nO);
			for(m=0; m<param.nO; m++)
				fprintf(fidP,"%10.8f%s",local_pred[m],(m<(param.nO-1))?"\t":"\n");
			predict_idx++;
		}
		
		// update p(L)
		if(o>=0) { // known observation
			pLe_denom = 0.0;
			// 1. pLe =  (L .* B(:,o)) ./ ( L'*B(:,o)+1e-8 );
			for(i=0; i<param.nS; i++) pLe_denom += group_skill_map[g][k][i] * a_B[i][o];
			for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i] * a_B[i][o] / safe0num(pLe_denom);
			// 2. L = (pLe'*A)';
			for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0;
			for(j=0; j<param.nS; j++)
				for(i=0; i<param.nS; i++)
					group_skill_map[g][k][j] += pLe[i] * a_A[i][j];
		} else { // unknown observation
			// 2. L = (pL'*A)';
			for(i=0; i<param.nS; i++) pLe[i] = group_skill_map[g][k][i]; // copy first;
			for(i=0; i<param.nS; i++) group_skill_map[g][k][i] = 0.0; // erase old value
			for(j=0; j<param.nS; j++)
				for(i=0; i<param.nS; i++)
					group_skill_map[g][k][j] += pLe[i] * a_A[i][j];
		}
	} // for all data
	free(line);
	free(local_pred);
	fclose(fid);
	fclose(fidP);
}
