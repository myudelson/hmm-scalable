//
//  InputUtil.cpp
//  HMM
//
//  Created by Yudelson, Michael on 7/17/13.
//
//

#include "InputUtil.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <list>

static int max_line_length;
static char * line;

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

void InputUtil::writeString(FILE *f, string str) {
    // Create char pointer from string.
    char* text = const_cast<char*>(str.c_str());
    // Find the length of the string.
    unsigned int size = (unsigned int)str.size();
    // Write the string's size to the file.
    fwrite(&size, sizeof(unsigned int), 1, f);
    // Followed by the string itself.
    fwrite(text, 1, size, f);
}

unsigned long InputUtil::writeMultiSkill(FILE *f, struct param * param) {
    if(param->multiskill == 0) {
        fprintf(stderr,"Error: multiskill flag should not be 0\n.");
        return 0;
    }
    NCAT * ar;
    unsigned long all_nwrit = 0;
    unsigned long nwrit     = 0;
    for(unsigned int t=0; t<param->N; t++) {
        ar = param->dat_multiskill->get(t);
        nwrit = (unsigned int)fwrite (ar, sizeof(NCAT), ar[0]+1, f);
        if(nwrit != ar[0]+1) {
            fprintf(stderr,"Errr while writing element %d of the multi-skill data. Expected %d, written %lu\n.",t,ar[0]+1,nwrit);
            return 0;
        }
        all_nwrit += nwrit;
    }
    return all_nwrit;
}

unsigned long InputUtil::readMultiSkill(FILE *f, struct param * param) {
    if(param->multiskill == 0) {
        fprintf(stderr,"Error: multiskill flag should not be 0\n.");
        return 0;
    }
    NCAT * ar;
    NCAT n;
    unsigned long all_nread = 0;
    unsigned long nread     = 0;
    for(unsigned int t=0; t<param->N; t++) {
        // read count
        nread = fread(&n, sizeof(NCAT), 1, f);
        if( nread!= 1) {
            fprintf(stderr,"Error: read a wrong number of datapoints\n.");
            return 0;
        }
        ar = Calloc(NCAT, n+1);
        ar[0] = n;
        // read data
        nread = fread(&ar[1], sizeof(NCAT), n, f);
        if( nread!= n) {
            fprintf(stderr,"Error: read a wrong number of datapoints\n.");
            return 0;
        }
        all_nread += (nread+1);
        // place
        param->dat_multiskill->set(t, ar);
    }
    return all_nread;
}

string InputUtil::readString(FILE *f) {
    // Create new string object to store the retrieved text and to return to the calling function.
    string str;
    // UInt for storing the string's size.
    unsigned int size;
    // Read the size of the string from the file and store it in size.
    fread(&size, sizeof(unsigned int), 1, f);
    // Create a char pointer for temporary storage.
    char* text = new char[size];
    // Read [size] number of characters from the string and store them in text.
    fread(text, 1, size, f);
    // Store the contents of text in str.
    str = text;
    // Resize str to match the size else we get extra cruft (line endings methinks).
    str.resize(size);
    // Finally, return the string to the calling function.
    delete[] text;
    return str;
}

bool InputUtil::readTxt(const char *fn, struct param * param) {
	FILE *fid = fopen(fn,"r");
    if( fid == NULL) {
        fprintf(stderr,"Could not read input file (%s).\n",fn);
        return false;
    }
    
	int number_columns = 0;
	max_line_length = 1024;
	char *col;
    
	// count lines and check for number of columns
	line = (char *)malloc(max_line_length);// Malloc(char,max_line_length);
	
	// grab memory and read all data
	param->dat_obs   = new StripedArray<NPAR>();
	param->dat_group = new StripedArray<NCAT>();
    if(param->multiskill==0)
        param->dat_skill = new StripedArray<NCAT>();
    else
        param->dat_multiskill = new StripedArray< NCAT* >(true);
	param->dat_item = new StripedArray<NCAT2>();
    if(param->time)
        param->dat_time = new StripedArray<int>();
    param->map_group_fwd = new map<string,NCAT>();
    param->map_group_bwd = new map<NCAT,string>();
    param->map_skill_fwd = new map<string,NCAT>();
    param->map_skill_bwd = new map<NCAT,string>();
    param->map_step_fwd = new map<string,NCAT2>();
    param->map_step_bwd = new map<NCAT2,string>();
	string s_group, s_step, s_skill;
    int time = 0;
	map<string,NCAT>::iterator it;
	map<string,NCAT2>::iterator it2;
	bool wrong_no_columns = false;
    param->N = 0;
    param->N_null = 0;
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
			return false;
		}
		param->dat_obs->add((NPAR)obs); // dat_obs[t] = (NPAR)obs;
		if( (obs >= 0) && ((param->nO-1) < obs) )
			param->nO = obs + (NPAR)1; // obs[t] + 1;
		// Group
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_group = string( col );
		it = param->map_group_fwd->find(s_group);
		if( it==param->map_group_fwd->end() ) { // not found
			if(param->map_group_fwd->size()==NCAT_MAX) {
				fprintf(stderr,"Number of unique groups exceeds allowed maximum of %d.\n",NCAT_MAX);
				return false;
			}
			NCAT newg = param->map_group_fwd->size();
			param->dat_group->add(newg); //[t] = param->map_group_fwd.size();
			param->map_group_fwd->insert(pair<string,NCAT>(s_group, newg));
			param->map_group_bwd->insert(pair<NCAT,string>(newg, s_group));
		}
		else
			param->dat_group->add(it->second); // [t] = it->second;
		
		
		// Step
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_step = string( col );
		it2 = param->map_step_fwd->find(s_step);
		if( it2==param->map_step_fwd->end() ) { // not found
			if(param->map_step_fwd->size()==NCAT2_MAX) {
				fprintf(stderr,"Number of unique steps exceeds allowed maximum of %d.\n",NCAT2_MAX);
				return false;
			}
            NCAT2 news = param->map_step_fwd->size();
			param->dat_item->add(news); //[t] = param->map_group_fwd.size();
			param->map_step_fwd->insert(pair<string,NCAT2>(s_step, news));
			param->map_step_bwd->insert(pair<NCAT2,string>(news, s_step));
		}
		else
            param->dat_item->add(it2->second);
		
        
		// Skill
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_skill = string( col );
        // process skills later
        // Time
        if(param->time) {
            col = strtok(NULL,"\t\n\r");
            if(col == NULL) {
                wrong_no_columns = true;
                break;
            }
            number_columns++;
            time = atoi( col );
            if( time<=0 ) {
                fprintf(stderr,"Time cannot be negative or zero (line %d).\n",param->N+1);
				return false;
            }
            param->dat_time->add(time);
        } // time
        
        // back to skill processing
		if( (s_skill.empty() || ( s_skill.size()==1 && (s_skill[0]=='.' || s_skill[0]==' ') ) ) ) { // null skill
            param->N_null++;
            if(param->multiskill == 0) {
                param->dat_skill->add(-1); // [t] = -1;
            }
            else {
                NCAT* a_skills = Malloc(NCAT, 2);
                a_skills[0] = 1; // count
                a_skills[1] = -1; // value
                param->dat_multiskill->add(a_skills);
            }
		} // empty skill
		else { // non-empty skill
            // multiskill
            if(param->multiskill != 0) {
                list<NCAT> a_skills;//
                char* a_kc;
                col = &s_skill[0];//. c_str();
                a_kc  = strtok(col, "~\n\r");
                string s_kc;
                while(a_kc != NULL) {
                    s_kc = string(a_kc);
                    // adding vvvv
                    it = param->map_skill_fwd->find(s_kc);
                    if( it==param->map_skill_fwd->end() ) { // not found
                        if(param->map_skill_fwd->size()==NCAT_MAX) {
                            fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                            return false;
                        }
                        a_skills.insert(a_skills.end(), param->map_skill_fwd->size()); //dat_skill->add(param->map_skill_fwd->size());
                        param->map_skill_fwd->insert(pair<string,NCAT>(s_kc, param->map_skill_fwd->size()));
                        param->map_skill_bwd->insert(pair<NCAT,string>(param->map_skill_bwd->size(),s_kc));
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
                param->dat_multiskill->add(b_skills);
                // multi skill
            } else {
                // single skill
                it = param->map_skill_fwd->find(s_skill);
                if( it==param->map_skill_fwd->end() ) { // not found
                    if(param->map_skill_fwd->size()==NCAT_MAX) {
                        fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                        return false;
                    }
                    param->dat_skill->add(param->map_skill_fwd->size()); //[t] = param->map_skill_fwd.size();
                    param->map_skill_fwd->insert(pair<string,NCAT>(s_skill, param->map_skill_fwd->size()));
                    param->map_skill_bwd->insert(pair<NCAT,string>(param->map_skill_bwd->size(),s_skill));
                }
                else
                    param->dat_skill->add(it->second); //[t] = it->second;
            } // single skill
		} // non empty skill
        
		// count lines
		param->N++;	// increase line count
        //        fprintf(stdout,"Line %d\n",param->N);
	}// reading loop
	if(wrong_no_columns) {
		fprintf(stderr,"Wrong number of columns in line %u. Expected %d, found %d\n",param->N+1,COLUMNS+param->time, number_columns);
		free(line);
		fclose(fid);
        return false;
	}
	param->nG = (NCAT)param->map_group_fwd->size();
	param->nK = (NCAT)param->map_skill_fwd->size();
	param->nI = (NCAT2)param->map_step_fwd->size();
	fclose(fid);
	free(line);
    return true;
}

bool InputUtil::readBin(const char *fn, struct param * param) {
    char c;
    unsigned int i;
    unsigned long nread;
    FILE *fid = fopen(fn,"rb");
    
    // version
    nread = fread (&c, sizeof(char), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading version data from %s\n",fn);
        return false;
    }
    if( c != bin_input_file_verstion) {
        fprintf(stderr,"Wrong version of the data file. Expected %c, actual %c\n",(int)bin_input_file_verstion, c);
        return false;
    }
    
    // N
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading N from %s\n",fn);
        return false;
    }
    param->N = i;
    
    // N_null
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading N_null from %s\n",fn);
        return false;
    }
    param->N_null = i;
    
    // nO
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading nO from %s\n",fn);
        return false;
    }
    param->nO = i;
    
    // nG
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading nG from %s\n",fn);
        return false;
    }
    param->nG = i;
    
    // nI
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading nI from %s\n",fn);
        return false;
    }
    param->nI = i;
    
    // nK
    nread = fread (&i, sizeof(unsigned int), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading nK from %s\n",fn);
        return false;
    }
    param->nK = i;
    
    // multiskill
    nread = fread (&c, sizeof(char), 1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading multiskill flag from %s\n",fn);
        return false;
    }
    param->multiskill = c;
    
    // dat_obs
    param->dat_obs = new StripedArray<NPAR>(fid, param->N);
    // dat_group
    param->dat_group = new StripedArray<NCAT>(fid, param->N);
    // dat_skill
    if(param->multiskill == 0)
        param->dat_skill = new StripedArray<NCAT>(fid, param->N);
    else {
        param->dat_multiskill = new StripedArray< NCAT* >(param->N,true);
        nread = readMultiSkill(fid, param);
    }
    // dat_item
    param->dat_item = new StripedArray<NCAT2>(fid, param->N);
    
    
    string str;
    param->map_group_fwd = new map<string,NCAT>();
    param->map_group_bwd = new map<NCAT,string>();
    param->map_skill_fwd = new map<string,NCAT>();
    param->map_skill_bwd = new map<NCAT,string>();
    param->map_step_fwd = new map<string,NCAT2>();
    param->map_step_bwd = new map<NCAT2,string>();
    // voc_group
    for(NCAT g=0; g<param->nG; g++) {
        str = readString(fid);
        param->map_group_fwd->insert(pair<string,NCAT>(str, g));
        param->map_group_bwd->insert(pair<NCAT,string>(g, str));
    }
    
    // voc_skill
    for(NCAT k=0; k<param->nK; k++) {
        str = readString(fid);
        param->map_skill_fwd->insert(pair<string,NCAT>(str, k));
        param->map_skill_bwd->insert(pair<NCAT,string>(k, str));
    }
    
    // voc_item
    for(NCAT2 i=0; i<param->nI; i++) {
        str = readString(fid);
        param->map_step_fwd->insert(pair<string,NCAT2>(str, i));
        param->map_step_bwd->insert(pair<NCAT2,string>(i, str));
    }
    map<NCAT2,string>::iterator it2;
    for (it2 =  param->map_step_bwd->begin(); it2 != param->map_step_bwd->end(); ++it2) {
        // fprintf(stderr,"key: %u, value: %s\n",it2->first, it2->second.c_str());
        fwrite (it2->second.c_str() , sizeof(std::string), 1, fid);
    }
    
    fclose(fid);
    return true;
}

/*
 * File format:
 *  - version number: char 1:...
 *  - N : unsigned int 1...
 *  - N_null : unsigned int 0...
 *  - nO : unsigned int 1...
 *  - nG : unsigned int 1...
 *  - nI : unsigned int 1...
 *  - nK : unsigned int 1...
 *  - multiskill : char (0, 1)
 *  - dat_obs : char * N
 *  - dat_group : unsigned int * N
 *  - dat_skill :
 *      a) single skill: unsigned int * N
 *      b) multiple skill: unsigned int * N * skills' (' variable)
 *  - dat_item : unsigned int * N
 *  - voc_group : string * nG : ordered by 1:nG
 *  - voc_skill : string * nK : ordered by 1:nK
 *  - voc_item  : string * nI : ordered by 1:nI
 */

bool InputUtil::toBin(struct param * param, const char *fn) {
    char c;
    unsigned int i;
    FILE *fid = fopen(fn,"wb");

    // version
    c = bin_input_file_verstion;
    fwrite (&c , sizeof(char), 1, fid);
    
    // N
    i = param->N;
    fwrite (&i , sizeof(unsigned int), 1, fid);

    // N_null
    i = param->N_null;
    fwrite (&i , sizeof(unsigned int), 1, fid);
    
    // nO
    i = param->nO;
    fwrite (&i , sizeof(unsigned int), 1, fid);
    
    // nG
    i = param->nG;
    fwrite (&i , sizeof(unsigned int), 1, fid);
    
    // nI
    i = param->nI;
    fwrite (&i , sizeof(unsigned int), 1, fid);
    
    // nK
    i = param->nK;
    fwrite (&i , sizeof(unsigned int), 1, fid);
    
    // multiskill
    c = param->multiskill;
    fwrite (&c , sizeof(char), 1, fid);
    
    unsigned long nwrit;
    // dat_obs
    nwrit = param->dat_obs->toBinFile(fid);
    // dat_group
    nwrit = param->dat_group->toBinFile(fid);
    // dat_skill
    if(param->multiskill == 0)
        nwrit = param->dat_skill->toBinFile(fid);
    else
        nwrit = writeMultiSkill(fid, param);
    // dat_item
    nwrit = param->dat_item->toBinFile(fid);

    map<NCAT,string>::iterator it;
    // voc_group
    for (it = param->map_group_bwd->begin(); it != param->map_group_bwd->end(); ++it) {
        // fprintf(stderr,"key: %u, value: %s\n",it->first, it->second.c_str());
        writeString(fid, it->second);
    }

    // voc_skill
    for (it = param->map_skill_bwd->begin(); it != param->map_skill_bwd->end(); ++it) {
        // fprintf(stderr,"key: %u, value: %s\n",it->first, it->second.c_str());
        writeString(fid, it->second);
    }
    
    // voc_item
    map<NCAT2,string>::iterator it2;
    for (it2 =  param->map_step_bwd->begin(); it2 != param->map_step_bwd->end(); ++it2) {
        // fprintf(stderr,"key: %u, value: %s\n",it2->first, it2->second.c_str());
        writeString(fid, it->second);
    }
    
    fclose(fid);
    return true;
}
