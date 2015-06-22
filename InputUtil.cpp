/*
 
 Copyright (c) 2012-2015, Michael (Mikhail) Yudelson
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

#include "InputUtil.h"
#include "utils.h"
#include <stdio.h>
#include <map>
#include <list>

#include <fstream>
#include <string>
#include <sstream>
#include <iostream>


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

void InputUtil::writeString(FILE *f, string str) {
    // Create char pointer from string.
    char* text = const_cast<char*>(str.c_str());
    // Find the length of the string.
    NDAT size = (NDAT)str.size();
    // Write the string's size to the file.
    fwrite(&size, sizeof(NDAT), 1, f);
    // Followed by the string itself.
    fwrite(text, 1, (size_t)size, f);
}

//NDAT InputUtil::writeMultiSkill(FILE *f, struct param * param) {
//    if(param->multiskill == 0) {
//        fprintf(stderr,"Error: multiskill flag should not be 0\n.");
//        return 0;
//    }
//    NCAT * ar;
//    NDAT all_nwrit = 0;
//    NDAT nwrit     = 0;
//    for(NDAT t=0; t<param->N; t++) {
//        ar = param->dat_multiskill->get(t);
//        nwrit = (NDAT)fwrite (ar, sizeof(NCAT), (size_t)ar[0]+1, f);
//        if(nwrit != ar[0]+1) {
//            fprintf(stderr,"Errr while writing element %d of the multi-skill data. Expected %d, written %u\n.",t,ar[0]+1,nwrit);
//            return 0;
//        }
//        all_nwrit += nwrit;
//    }
//    return all_nwrit;
//}
//
//NDAT InputUtil::readMultiSkill(FILE *f, struct param * param, char version) {
//    if(param->multiskill == 0) {
//        fprintf(stderr,"Error: multiskill flag should not be 0\n.");
//        return 0;
//    }
//    NCAT * ar;
//    short * arv1;
//    NCAT n;
//    NDAT all_nread = 0;
//    NDAT nread     = 0;
//    for(NDAT t=0; t<param->N; t++) {
//        // read count
//        if(version == 1)
//            nread = (NCAT)fread(&n, sizeof(short), 1, f);
//        else
//            nread = (NDAT)fread(&n, sizeof(NCAT), (size_t)1, f);
//        if( nread!= 1) {
//            fprintf(stderr,"Error: read a wrong number of datapoints\n.");
//            return 0;
//        }
//        ar = Calloc(NCAT, (size_t)n+1);
//        ar[0] = n;
//        if(version==1)
//            arv1 = Calloc(short, (size_t)n+1);;
//        // read data
//        if(version == 1) {
//            nread = (NCAT)fread(&arv1[1], sizeof(short), (size_t)n, f);
//            for(int i=1;i<(n+1);i++)
//                ar[i] = (NCAT)arv1[i];
//        }
//        else
//            nread = (NDAT)fread(&ar[1], sizeof(NCAT), (size_t)n, f);
//        if( nread!= n) {
//            fprintf(stderr,"Error: read a wrong number of datapoints\n.");
//            return 0;
//        }
//        all_nread += (nread+1);
//        // place
//        param->dat_multiskill->set(t, ar);
//    }
//    return all_nread;
//}

string InputUtil::readString(FILE *f) {
    // Create new string object to store the retrieved text and to return to the calling function.
    string str;
    // UInt for storing the string's size.
    NDAT size;
    // Read the size of the string from the file and store it in size.
    fread(&size, sizeof(NDAT), 1, f);
    // Create a char pointer for temporary storage.
    char* text = new char[size];
    // Read [size] number of characters from the string and store them in text.
    fread(text, 1, (size_t)size, f);
    // Store the contents of text in str.
    str = text;
    // Resize str to match the size else we get extra cruft (line endings methinks).
    str.resize((size_t)size);
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
	line = (char *)malloc((size_t)max_line_length);// Malloc(char,max_line_length);
	
	// grab memory and read all data
	StripedArray<NPAR> *striped_dat_obs = new StripedArray<NPAR>();
	StripedArray<NCAT> *striped_dat_group = new StripedArray<NCAT>();
    StripedArray<NCAT> *striped_dat_skill = NULL;
    StripedArray<NCAT> *striped_dat_skill_stacked = NULL;
    StripedArray<NCAT> *striped_dat_skill_rcount  = NULL;
    StripedArray<NCAT> *striped_dat_skill_rix     = NULL;
    NDAT current_stacked = 0;
    if(param->multiskill==0)
        striped_dat_skill = new StripedArray<NCAT>();
    else {
//        param->dat_multiskill = new StripedArray< NCAT* >(true);
        striped_dat_skill_stacked = new StripedArray<NCAT>();
        striped_dat_skill_rcount  = new StripedArray<NCAT>();
        striped_dat_skill_rix     = new StripedArray<NCAT>();
    }
	StripedArray<NCAT> * striped_dat_item = new StripedArray<NCAT>();
    StripedArray<NPAR> *striped_dat_slice = NULL;
    if(param->sliced)
        striped_dat_slice = new StripedArray<NPAR>();
    param->map_group_fwd = new map<string,NCAT>();
    param->map_group_bwd = new map<NCAT,string>();
    param->map_skill_fwd = new map<string,NCAT>();
    param->map_skill_bwd = new map<NCAT,string>();
    param->map_step_fwd = new map<string,NCAT>();
    param->map_step_bwd = new map<NCAT,string>();
	string s_group, s_step, s_skill;
    NPAR slice = 0;
	map<string,NCAT>::iterator it;
	map<string,NCAT>::iterator it2;
	bool wrong_no_columns = false;
    param->N = 0;
    param->Nstacked = 0;
    NDAT Nstacked_alt = 0;
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
		striped_dat_obs->add(obs); // dat_obs[t] = (NPAR)obs;
		if( (obs >= 0) && ((param->nO-1) < obs) )
			param->nO = (NPAR)(obs + 1); // obs[t] + 1;

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
			NCAT newg = (NCAT)param->map_group_fwd->size();
			striped_dat_group->add(newg); //[t] = param->map_group_fwd.size();
			param->map_group_fwd->insert(pair<string,NCAT>(s_group, newg));
			param->map_group_bwd->insert(pair<NCAT,string>(newg, s_group));
		}
		else
			striped_dat_group->add(it->second); // [t] = it->second;
		
		
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
			if(param->map_step_fwd->size()==NCAT_MAX) {
				fprintf(stderr,"Number of unique steps exceeds allowed maximum of %d.\n",NCAT_MAX);
				return false;
			}
            NCAT news = (NCAT)param->map_step_fwd->size();
			striped_dat_item->add(news); //[t] = param->map_group_fwd.size();
			param->map_step_fwd->insert(pair<string,NCAT>(s_step, news));
			param->map_step_bwd->insert(pair<NCAT,string>(news, s_step));
		}
		else
            striped_dat_item->add(it2->second);
		
        
		// Skill
		col = strtok(NULL,"\t\n\r");
		if(col == NULL) {
			wrong_no_columns = true;
			break;
		}
		number_columns++;
		s_skill = string( col );
        // process skills later
        
        // Slice
        if(param->sliced) {
            col = strtok(NULL,"\t\n\r");
            if(col == NULL) {
                wrong_no_columns = true;
                break;
            }
            number_columns++;
            slice = (NPAR) atoi( col );
            if( slice<0 ) {
                fprintf(stderr,"Slice cannot be negative (line %d).\n",param->N+1);
				return false;
            }
            if( (slice >= 0) && ((param->nZ-1) < slice) ) // update slice count
                param->nZ = (NPAR)(slice + 1);
            striped_dat_slice->add(slice);
        } // time
        
        // back to skill processing
		if( (s_skill.empty() || ( s_skill.size()==1 && (s_skill[0]=='.' || s_skill[0]==' ') ) ) ) { // null skill
            param->N_null++;
            param->Nstacked++;	// increase stackd count too
            if(param->multiskill == 0) {
                striped_dat_skill->add(-1); // [t] = -1;
            }
            else {
                NCAT* a_skills = Malloc(NCAT, 2);
                a_skills[0] = 1; // count
                a_skills[1] = -1; // value
//                param->dat_multiskill->add(a_skills);
                // add whole multi-skill row
                striped_dat_skill_rcount->add(1);
                striped_dat_skill_rix->add(current_stacked++); // first add pointer, then increase
                striped_dat_skill_stacked->add(-1);
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
                // stacked
                NDAT skill_count = 0;
                striped_dat_skill_rix->add(current_stacked); // add stacked pointer once (but increase it with every next skill
                while(a_kc != NULL) {
                    s_kc = string(a_kc);
                    // stacked
                    skill_count++;
                    current_stacked++;
                    param->Nstacked++;	// increase line count
                    // adding vvvv
                    it = param->map_skill_fwd->find(s_kc);
                    if( it==param->map_skill_fwd->end() ) { // not found
                        if(param->map_skill_fwd->size()==NCAT_MAX) {
                            fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                            return false;
                        }
                        a_skills.insert(a_skills.end(), (NCAT)param->map_skill_fwd->size()); //dat_skill->add(param->map_skill_fwd->size());
                        param->map_skill_fwd->insert(pair<string,NCAT>(s_kc, (NCAT)param->map_skill_fwd->size()));
                        param->map_skill_bwd->insert(pair<NCAT,string>((NCAT)param->map_skill_bwd->size(),s_kc));
                        // add stacked skill
                        striped_dat_skill_stacked->add((NCAT)param->map_skill_fwd->size()-1); // -1 because after adding
                    }
                    else {
                        a_skills.insert(a_skills.end(), it->second); //dat_skill->add(it->second); //[t] = it->second;
                        // add stacked skill
                        striped_dat_skill_stacked->add(it->second);
                    }
                    // adding ^^^^
                    a_kc  = strtok(NULL,"~\n\r");
                }
//                NCAT *b_skills = Malloc(NCAT, a_skills.size()+1);
//                b_skills[0] = (NCAT)a_skills.size();
//                int count = 0;
//                for(list<NCAT>::iterator it=a_skills.begin(); it!=a_skills.end(); it++)
//                    b_skills[++count] = *it;
//                param->dat_multiskill->add(b_skills);
                // add stacked count
                striped_dat_skill_rcount->add(skill_count);
                Nstacked_alt+=skill_count;
                // multi skill
            } else {
                // single skill
                it = param->map_skill_fwd->find(s_skill);
                if( it==param->map_skill_fwd->end() ) { // not found
                    if(param->map_skill_fwd->size()==NCAT_MAX) {
                        fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                        return false;
                    }
                    striped_dat_skill->add((NCAT)param->map_skill_fwd->size()); //[t] = param->map_skill_fwd.size();
                    param->map_skill_fwd->insert(pair<string,NCAT>(s_skill, (NCAT)param->map_skill_fwd->size()));
                    param->map_skill_bwd->insert(pair<NCAT,string>((NCAT)param->map_skill_bwd->size(),s_skill));
                }
                else
                    striped_dat_skill->add(it->second); //[t] = it->second;
            } // single skill
		} // non empty skill
        
		// count lines
		param->N++;	// increase line count
	}// reading loop
	if(wrong_no_columns) {
		fprintf(stderr,"Wrong number of columns in line %u. Expected %d, found %d\n",param->N+1,COLUMNS+param->sliced, number_columns);
		free(line);
		fclose(fid);
        return false;
	}
	param->nG = (NCAT)param->map_group_fwd->size();
	param->nK = (NCAT)param->map_skill_fwd->size();
	param->nI = (NCAT)param->map_step_fwd->size();
    
    // copy striped to lined
    param->dat_obs = striped_dat_obs->toArray();
  
    delete striped_dat_obs;
    param->dat_group = striped_dat_group->toArray();
    delete striped_dat_group;
    if(param->multiskill==0) {
        param->dat_skill = striped_dat_skill->toArray();
        delete striped_dat_skill;
    } else {
        param->dat_skill_stacked = striped_dat_skill_stacked->toArray();
        param->dat_skill_rcount  = striped_dat_skill_rcount->toArray();
        param->dat_skill_rix     = striped_dat_skill_rix->toArray();
        delete striped_dat_skill_stacked;
        delete striped_dat_skill_rcount;
        delete striped_dat_skill_rix;
    }
    param->dat_item = striped_dat_item->toArray();
    delete striped_dat_item;
    if(param->sliced) {
        param->dat_slice = striped_dat_slice->toArray();
        delete striped_dat_slice;
    }

	fclose(fid);
	free(line);
    return true;
}

bool InputUtil::readBin(const char *fn, struct param * param) {
    char c, v/*version*/;
    NDAT i;
    NDAT nread;
    FILE *fid = fopen(fn,"rb");
    
    // version
    nread = (NDAT)fread (&v, sizeof(char), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading version data from %s\n",fn);
        return false;
    }
    if( v > bin_input_file_verstion) { // we will handle earlier versions
        fprintf(stderr,"Wrong version of the data file. Expected %d, actual %d\n",(int)bin_input_file_verstion, (int)v);
        return false;
    }
    
    // N
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of rows from %s\n",fn);
        return false;
    }
    param->N = (NDAT)i;
    
    // Nstacked
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of rows (stacked) from %s\n",fn);
        return false;
    }
    param->Nstacked = (NDAT)i;
    
    // N_null
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of rows with null skills from %s\n",fn);
        return false;
    }
    param->N_null = (NDAT)i;
    
    // nO
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of observations from %s\n",fn);
        return false;
    }
    param->nO = (NPAR)i;
    
    // nG
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of groups (students) from %s\n",fn);
        return false;
    }
    param->nG = (NCAT)i;
    
    // nI
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of items from %s\n",fn);
        return false;
    }
    param->nI = (NCAT)i;
    
    // nK
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of skills from %s\n",fn);
        return false;
    }
    param->nK = (NCAT)i;
    
    // nZ
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of slices from %s\n",fn);
        return false;
    }
    param->nZ = (NPAR)i;
    if(param->nZ<1) {
        fprintf(stderr,"Number of slices should be at least 1\n");
        return true;
    }
    
    // multiskill
    nread = (NDAT)fread (&c, sizeof(char), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading multiskill flag from %s\n",fn);
        return false;
    }
    param->multiskill = (NPAR)c;
    
    
    // dat_obs
    StripedArray<NPAR> *striped_dat_obs = new StripedArray<NPAR>(fid, param->N);
    param->dat_obs = striped_dat_obs->toArray();
    delete striped_dat_obs;

//    if(v==1) { // older NCAT of unsigned short
//        StripedArray<short> *dat_group_legacy = NULL;
//        StripedArray<short> *dat_skill_legacy = NULL;
//        StripedArray<short*> *dat_multiskill_legacy = NULL;
//        dat_group_legacy = new StripedArray<short>(fid, param->N);
//        for(t=0; t<param->N; t++)
//            param->dat_group[t] = (NCAT)dat_group_legacy->get(t);
//        if(param->multiskill == 0) {
//            dat_skill_legacy = new StripedArray<short>(fid, param->N);
//            for(t=0; t<param->N; t++)
//                param->dat_skill[t] = (NCAT)dat_skill_legacy->get(t);
//            delete dat_skill_legacy;
//        } else {
//            param->dat_multiskill = new StripedArray<NCAT*>(param->N,true);
//            nread = readMultiSkill(fid, param, v);
//            delete dat_multiskill_legacy;
//        }
//        delete dat_group_legacy;
//    } else {
        // dat_group
        StripedArray<NCAT> *striped_dat_group = new StripedArray<NCAT>(fid, param->N);
        param->dat_group = striped_dat_group->toArray();
        delete striped_dat_group;
//    }
    
    // dat_item
    StripedArray<NCAT> *striped_dat_item = new StripedArray<NCAT>(fid, param->N);
    param->dat_item = striped_dat_item->toArray();
    delete striped_dat_item;
    
    // dat_skill
    if(param->multiskill == 0) {
        StripedArray<NCAT> *striped_dat_skill = new StripedArray<NCAT>(fid, param->N);
        param->dat_skill = striped_dat_skill->toArray();
        delete striped_dat_skill;
    } else {
//            param->dat_multiskill = new StripedArray< NCAT* >(param->N,true);
//            nread = readMultiSkill(fid, param, v);
        StripedArray<NCAT> *striped_dat_skill_stacked = new StripedArray<NCAT>(fid, param->Nstacked);
        param->dat_skill_stacked = striped_dat_skill_stacked->toArray();
        StripedArray<NCAT> *striped_dat_skill_rcount = new StripedArray<NCAT>(fid, param->N);
        param->dat_skill_rcount = striped_dat_skill_rcount->toArray();
        StripedArray<NCAT> *striped_dat_skill_rix = new StripedArray<NCAT>(fid, param->N);
        param->dat_skill_rix = striped_dat_skill_rix->toArray();
        delete striped_dat_skill_stacked;
        delete striped_dat_skill_rcount;
        delete striped_dat_skill_rix;
    }
    // dat_slices, only of nZ > 1
    if(param->nZ > 1) {
        NDAT szZ = (param->multiskill == 0)?param->N:param->Nstacked;
        StripedArray<NPAR> *striped_dat_slice = new StripedArray<NPAR>(fid, szZ);
        param->dat_slice = striped_dat_slice->toArray();
        delete striped_dat_slice;
    }
        
    string str;
    param->map_group_fwd = new map<string,NCAT>();
    param->map_group_bwd = new map<NCAT,string>();
    param->map_skill_fwd = new map<string,NCAT>();
    param->map_skill_bwd = new map<NCAT,string>();
    param->map_step_fwd = new map<string,NCAT>();
    param->map_step_bwd = new map<NCAT,string>();
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
    for(NCAT i=0; i<param->nI; i++) {
        str = readString(fid);
        param->map_step_fwd->insert(pair<string,NCAT>(str, i));
        param->map_step_bwd->insert(pair<NCAT,string>(i, str));
    }
    
    fclose(fid);
    return true;
}

/*
 * File format:
 *  - version number: char 1:...
 *  - N : NDAT 1...
 *  - N_stacked : NDAT 1...
 *  - N_null : NDAT 0...
 *  - nO : NDAT 1...
 *  - nG : NDAT 1...
 *  - nI : NDAT 1...
 *  - nK : NDAT 1...
 *  - multiskill : char (0, 1)
 *  - dat_obs : char * N
 *  - dat_group : NDAT * N
 *  - dat_skill :
 *      a) single skill: NDAT * N
 *      b) multiple skill: NDAT * N * skills' (' variable)
 *  - dat_item : NDAT * N
 *  - voc_group : string * nG : ordered by 1:nG
 *  - voc_skill : string * nK : ordered by 1:nK
 *  - voc_item  : string * nI : ordered by 1:nI
 */

bool InputUtil::toBin(struct param * param, const char *fn) {
    char c;
    NDAT i;
    FILE *fid = fopen(fn,"wb");

    // version
    c = bin_input_file_verstion;
    fwrite (&c , sizeof(char), 1, fid);
    
    // N
    i = param->N;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // Nstacked
    i = param->Nstacked;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // N_null
    i = param->N_null;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nO
    i = param->nO;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nG
    i = param->nG;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nI
    i = param->nI;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nK
    i = param->nK;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nZ
    i = param->nZ;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // multiskill
    c = param->multiskill;
    fwrite (&c , sizeof(char), 1, fid);
    
//    NDAT nwrit;
    // dat_obs
    /*nwrit = */StripedArray<NPAR>::arrayToBinFile(param->dat_obs, param->N, fid);//param->dat_obs->toBinFile(fid);
    // dat_group
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_group, param->N, fid);//->toBinFile(fid);
    // dat_item
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_item, param->N, fid);//->toBinFile(fid);
    // dat_skill
    if(param->multiskill == 0)
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_skill, param->N, fid);//->toBinFile(fid);
    else {
        //        /*nwrit = */writeMultiSkill(fid, param);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_skill_stacked, param->Nstacked, fid);//->toBinFile(fid);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_skill_rcount , param->N,        fid);//->toBinFile(fid);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(param->dat_skill_rix    , param->N,        fid);//->toBinFile(fid);
    }
    // dat_slices, only of nZ > 1
    if(param->nZ > 1) {
        NDAT szZ = (param->multiskill == 0)?param->N:param->Nstacked;
        StripedArray<NPAR>::arrayToBinFile(param->dat_slice, szZ, fid);
    }

    map<NCAT,string>::iterator it;
    // voc_group
    for (it = param->map_group_bwd->begin(); it != param->map_group_bwd->end(); ++it) {
        writeString(fid, it->second);
    }

    // voc_skill
    for (it = param->map_skill_bwd->begin(); it != param->map_skill_bwd->end(); ++it) {
        writeString(fid, it->second);
    }
    
    // voc_item
    for (it =  param->map_step_bwd->begin(); it != param->map_step_bwd->end(); ++it) {
        writeString(fid, it->second);
    }
    
    
    fclose(fid);
    return true;
}

// experimental
// for a skill, write all student sequences as a matrix
// Nstudents * Max attempts, 1 - correct, 2 - incorrect, 0 - empty
// space separated coumns
void InputUtil::writeInputMatrix(const char *filename, struct param* p, NCAT xndat, struct data** x_data) {
//    FILE *fid = fopen(filename,"w");
//    if(fid == NULL) {
//        fprintf(stderr,"Can't write output model file %s\n",filename);
//        exit(1);
//    }
    
    std::ofstream file;
    file.open(filename);
    
    NDAT nmax = 0;
    for(NCAT x=0; x<xndat; x++)
        if(nmax < x_data[x]->n)
            nmax = x_data[x]->n;
    
    for(NCAT x=0; x<xndat; x++) {
        std::stringstream ss;
        for(NDAT t=0; t<nmax; t++) {
            
            if(t<x_data[x]->n) {
                ss << ((t>0)?" ":"") << (int)(1+p->dat_obs[ x_data[x]->ix[t] ]);
            } else {
                ss << " " << 0;
            }
        }
        ss << "\n";
        file << ss.str();
    } // for all groups in skill
    
//    fclose(fid);
    file.close();
}
