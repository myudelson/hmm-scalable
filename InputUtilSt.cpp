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

#include "InputUtilSt.h"
#include "utilsSt.h"
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

void InputUtilSt::writeString(FILE *f, string str) {
    // Create char pointer from string.
    char* text = const_cast<char*>(str.c_str());
    // Find the length of the string.
    NDAT size = (NDAT)str.size();
    // Write the string's size to the file.
    fwrite(&size, sizeof(NDAT), 1, f);
    // Followed by the string itself.
    fwrite(text, 1, (size_t)size, f);
}

string InputUtilSt::readString(FILE *f) {
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

bool InputUtilSt::readTxt(const char *fn, struct task *task) {
    FILE *fid = fopen(fn,"r");
    if( fid == NULL) {
        fprintf(stderr,"Could not read input file (%s).\n",fn);
        return false;
    }
    
    int number_columns = 0;
    max_line_length = 1024;
    char *col;
    NDAT co = 0;
    
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
    if(task->multiskill==0)
        striped_dat_skill = new StripedArray<NCAT>();
    else {
//        task->dat_multiskill = new StripedArray< NCAT* >(true);
        striped_dat_skill_stacked = new StripedArray<NCAT>();
        striped_dat_skill_rcount  = new StripedArray<NCAT>();
        striped_dat_skill_rix     = new StripedArray<NCAT>();
    }
    StripedArray<NCAT> * striped_dat_item = new StripedArray<NCAT>();
    StripedArray<NPAR> *striped_dat_slice = NULL;
    if(task->sliced) {
        striped_dat_slice = new StripedArray<NPAR>();
    }
    task->map_group_fwd = new map<string,NCAT>();
    task->map_group_bwd = new map<NCAT,string>();
    task->map_skill_fwd = new map<string,NCAT>();
    task->map_skill_bwd = new map<NCAT,string>();
    task->map_step_fwd = new map<string,NCAT>();
    task->map_step_bwd = new map<NCAT,string>();
    string s_group, s_step, s_skill;
    NPAR slice = 0;
    map<string,NCAT>::iterator it;
    map<string,NCAT>::iterator it2;
    bool wrong_no_columns = false;
    task->N = 0;
    task->Nst = 0;
    NDAT Nst_alt = 0;
    task->N_null = 0;
    bool res = true;
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
            res=false;
            goto recycle;
        }
        striped_dat_obs->add(obs); // dat_obs[t] = (NPAR)obs;
        if( (obs >= 0) && ((task->nO-1) < obs) )
            task->nO = (NPAR)(obs + 1); // obs[t] + 1;

        // Group
        col = strtok(NULL,"\t\n\r");
        if(col == NULL) {
            wrong_no_columns = true;
            break;
        }
        number_columns++;
        s_group = string( col );
        it = task->map_group_fwd->find(s_group);
        if( it==task->map_group_fwd->end() ) { // not found
            if(task->map_group_fwd->size()==NCAT_MAX) {
                fprintf(stderr,"Number of unique groups exceeds allowed maximum of %d.\n",NCAT_MAX);
                res = false;
                goto recycle;
            }
            NCAT newg = (NCAT)task->map_group_fwd->size();
            striped_dat_group->add(newg); //[t] = task->map_group_fwd.size();
            task->map_group_fwd->insert(pair<string,NCAT>(s_group, newg));
            task->map_group_bwd->insert(pair<NCAT,string>(newg, s_group));
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
        it2 = task->map_step_fwd->find(s_step);
        if( it2==task->map_step_fwd->end() ) { // not found
            if(task->map_step_fwd->size()==NCAT_MAX) {
                fprintf(stderr,"Number of unique steps exceeds allowed maximum of %d.\n",NCAT_MAX);
                res=false;
                goto recycle;
            }
            NCAT news = (NCAT)task->map_step_fwd->size();
            striped_dat_item->add(news); //[t] = task->map_group_fwd.size();
            task->map_step_fwd->insert(pair<string,NCAT>(s_step, news));
            task->map_step_bwd->insert(pair<NCAT,string>(news, s_step));
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
        if(task->sliced) {
            col = strtok(NULL,"\t\n\r");
            if(col == NULL) {
                wrong_no_columns = true;
                break;
            }
            number_columns++;
            slice = (NPAR) atoi( col );
            if( slice<0 ) {
                fprintf(stderr,"Slice cannot be negative (line %d).\n",task->N+1);
                res=false;
                goto recycle;
            }
            if( (slice >= 0) && ((task->nZ-1) < slice) ) // update slice count
                task->nZ = (NPAR)(slice + 1);
            striped_dat_slice->add(slice);
        } // time
        
        // back to skill processing
        if( (s_skill.empty() || ( s_skill.size()==1 && (s_skill[0]=='.' || s_skill[0]==' ') ) ) ) { // null skill
            task->N_null++;
            task->Nst++;    // increase stackd count too
            if(task->multiskill == 0) {
                striped_dat_skill->add(-1); // [t] = -1;
            }
            else {
//                NCAT* a_skills = Malloc(NCAT, 2);
//                a_skills[0] = 1; // count
//                a_skills[1] = -1; // value
//                task->dat_multiskill->add(a_skills);
                // add whole multi-skill row
                striped_dat_skill_rcount->add(1);
                striped_dat_skill_rix->add(current_stacked++); // first add pointer, then increase
                striped_dat_skill_stacked->add(-1);
            }
        } // empty skill
        else { // non-empty skill
            // multiskill
            if(task->multiskill != 0) {
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
                    task->Nst++;    // increase line count
                    // adding vvvv
                    it = task->map_skill_fwd->find(s_kc);
                    if( it==task->map_skill_fwd->end() ) { // not found
                        if(task->map_skill_fwd->size()==NCAT_MAX) {
                            fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                            res = false;
                            goto recycle;
                        }
                        a_skills.insert(a_skills.end(), (NCAT)task->map_skill_fwd->size()); //dat_skill->add(task->map_skill_fwd->size());
                        task->map_skill_fwd->insert(pair<string,NCAT>(s_kc, (NCAT)task->map_skill_fwd->size()));
                        task->map_skill_bwd->insert(pair<NCAT,string>((NCAT)task->map_skill_bwd->size(),s_kc));
                        // add stacked skill
                        striped_dat_skill_stacked->add((NCAT)task->map_skill_fwd->size()-1); // -1 because after adding
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
//                task->dat_multiskill->add(b_skills);
                // add stacked count
                striped_dat_skill_rcount->add(skill_count);
                Nst_alt+=skill_count;
                // multi skill
            } else {
                // single skill
                it = task->map_skill_fwd->find(s_skill);
                if( it==task->map_skill_fwd->end() ) { // not found
                    if(task->map_skill_fwd->size()==NCAT_MAX) {
                        fprintf(stderr,"Number of unique skills exceeds allowed maximum of %d.\n",NCAT_MAX);
                        res = false;
                        goto recycle;
                    }
                    striped_dat_skill->add((NCAT)task->map_skill_fwd->size()); //[t] = task->map_skill_fwd.size();
                    task->map_skill_fwd->insert(pair<string,NCAT>(s_skill, (NCAT)task->map_skill_fwd->size()));
                    task->map_skill_bwd->insert(pair<NCAT,string>((NCAT)task->map_skill_bwd->size(),s_skill));
                }
                else
                    striped_dat_skill->add(it->second); //[t] = it->second;
            } // single skill
        } // non empty skill
        
        // count lines
        task->N++;    // increase line count
    }// reading loop
    if(wrong_no_columns) {
        fprintf(stderr,"Wrong number of columns in line %u. Expected %d, found %d\n",task->N+1,COLUMNS+task->sliced, number_columns);
        res = false;
        goto recycle;
    }
    task->nG = (NCAT)task->map_group_fwd->size();
    task->nK = (NCAT)task->map_skill_fwd->size();
    task->nI = (NCAT)task->map_step_fwd->size();
    
    // copy striped to lined
    task->dat_obs = striped_dat_obs->toArray();
    task->dat_obs_stacked = init1D<NPAR>(task->Nst);
    task->dat_predict = initToValue2D<NUMBER>(task->N, task->nO, 0);

    task->dat_group = striped_dat_group->toArray();
    if(task->multiskill==0) {
        task->dat_skill = striped_dat_skill->toArray();
    } else {
        task->dat_skill_stacked = striped_dat_skill_stacked->toArray();
        task->dat_skill_rcount  = striped_dat_skill_rcount->toArray();
        task->dat_skill_rix     = striped_dat_skill_rix->toArray();
    }
    task->dat_item = striped_dat_item->toArray();
    if(task->sliced) {
        task->dat_slice = striped_dat_slice->toArray();
    }
    // build dat_obs_stacked
    co = 0;
    if(task->multiskill!=0) {
        for(NDAT t=0;t<task->N;t++) {
            for(NDAT tt=0; tt<task->dat_skill_rcount[t]; tt++) {
                task->dat_obs_stacked[co] = task->dat_obs[t];
                co++;
            }
        }
    } else {
        memcpy( task->dat_obs_stacked, task->dat_obs , sizeof(NPAR)*(size_t)task->N );
    }

    //
    // create connectivities: skill-skill, student-skill
    //
    createConnectivities(task);

    recycle: delete striped_dat_obs;
    delete striped_dat_group;
    if(striped_dat_skill!=NULL) delete striped_dat_skill;
    if(striped_dat_skill_stacked!=NULL) delete striped_dat_skill_stacked;
    if(striped_dat_skill_rcount!=NULL) delete striped_dat_skill_rcount;
    if(striped_dat_skill_rix!=NULL) delete striped_dat_skill_rix;
    if(striped_dat_item!=NULL) delete striped_dat_item;
    if(striped_dat_slice!=NULL) delete striped_dat_slice;
    
    if(fid!=NULL) fclose(fid);
    free(line);
    return res;
}

void InputUtilSt::createConnectivities(struct task *a_task) {
    // init
    a_task->n_connectivities=2;
    a_task->n_connectivity_X = Calloc(NDAT, (size_t)a_task->n_connectivities);
    a_task->n_connectivity_X[0] = a_task->nK;
    a_task->n_connectivity_X[1] = a_task->nG;
    a_task->n_connectivity_Y = Calloc(NDAT, (size_t)a_task->n_connectivities);
    a_task->n_connectivity_Y[0] = a_task->nK;
    a_task->n_connectivity_Y[1] = a_task->nK;
    a_task->connectivities = Calloc(NPAR**, (size_t)a_task->n_connectivities);
    a_task->connectivities[0] = initToValue2D<NPAR>(a_task->n_connectivity_X[0], a_task->n_connectivity_Y[0], 0);
    a_task->connectivities[1] = initToValue2D<NPAR>(a_task->n_connectivity_X[1], a_task->n_connectivity_Y[1], 0);

    // create
    for(NDAT t=0; t<a_task->N; t++) { // vv for all non-stacked rows
        // skill-skill connectivity
        NCAT g = a_task->dat_group[t]; // -1, because in data they were 1-starting
        // grab skill array (if exists)
        NCAT *ar;
        NPAR n;
        getSkillsAtRow(a_task, t, &ar, &n);
        if(n>1) {
            for(NPAR i1=0;i1<(n-1);i1++) {
                for(NPAR i2=(NPAR)(i1+1);i2<n;i2++) {
                    a_task->connectivities[0][ ar[i1] ][ ar[i2] ] = 1;
                    a_task->connectivities[0][ ar[i2] ][ ar[i1] ] = 1;
                }
            }
        }
        // student-skill connectivity
        for(NPAR l=0; l<n; l++) {
            a_task->connectivities[1][g][l] = 1;
        }
    } // ^^ for all non-stacked rows
//    std::ofstream file; // DEBUG
//    file.open("connectivity.txt");  // DEBUG
//    for(NCAT i=0; i<a_task->n_connectivity_X[0]; i++) {  // DEBUG
//        std::stringstream ss;  // DEBUG
//        for(NDAT j=0; j<a_task->n_connectivity_Y[0]; j++) {  // DEBUG
//            ss << ((j>0)?" ":"") << std::to_string(a_task->connectivities[0][i][j]);  // DEBUG
//        }  // DEBUG
//        ss << "\n";  // DEBUG
//        file << ss.str();  // DEBUG
//    }  // DEBUG
//    file.close(); // DEBUG
    
    // turn connectivity[0] to reachability
    // vvv Warshall’s Algorithm $O(n^3)$ https://en.wikipedia.org/wiki/Reachability
    // https://cs.winona.edu/lin/cs440/ch08-2.pdf
    NDAT x1 = 0; // count 1's
    for(NCAT k=0; k<a_task->nK; k++) {
        x1 = 0;
        for(NCAT i=0; i<a_task->nK; i++) {
            for(NCAT j=0; j<a_task->nK; j++) {
                a_task->connectivities[0][i][j] = a_task->connectivities[0][i][j] || (a_task->connectivities[0][i][k] && a_task->connectivities[0][k][j]);
                x1+=a_task->connectivities[0][i][j];
            }
        }
    }// ^^^ Warshall’s Algorithm
//    std::ofstream file; // DEBUG
//    file.open("reachability.txt");  // DEBUG
//    for(NCAT i=0; i<a_task->n_connectivity_X[0]; i++) {  // DEBUG
//        std::stringstream ss;  // DEBUG
//        for(NDAT j=0; j<a_task->n_connectivity_Y[0]; j++) {  // DEBUG
//            ss << ((j>0)?" ":"") << std::to_string(a_task->connectivities[0][i][j]);  // DEBUG
//        }  // DEBUG
//        ss << "\n";  // DEBUG
//        file << ss.str();  // DEBUG
//    }  // DEBUG
//    file.close(); // DEBUG
}

bool InputUtilSt::readBin(const char *fn, struct task *task) {
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
    task->N = (NDAT)i;
    
    // Nst
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of rows (stacked) from %s\n",fn);
        return false;
    }
    task->Nst = (NDAT)i;
    
    // N_null
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of rows with null skills from %s\n",fn);
        return false;
    }
    task->N_null = (NDAT)i;
    
    // nO
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of observations from %s\n",fn);
        return false;
    }
    task->nO = (NPAR)i;
    
    // nG
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of groups (students) from %s\n",fn);
        return false;
    }
    task->nG = (NCAT)i;
    
    // nI
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of items from %s\n",fn);
        return false;
    }
    task->nI = (NCAT)i;
    
    // nK
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of skills from %s\n",fn);
        return false;
    }
    task->nK = (NCAT)i;
    
    // nZ
    nread = (NDAT)fread (&i, sizeof(NDAT), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading number of slices from %s\n",fn);
        return false;
    }
    task->nZ = (NPAR)i;
    if(task->nZ<1) {
        fprintf(stderr,"Number of slices should be at least 1\n");
        return true;
    }
    
    // multiskill
    nread = (NDAT)fread (&c, sizeof(char), (size_t)1, fid);
    if(nread != 1) {
        fprintf(stderr,"Error reading multiskill flag from %s\n",fn);
        return false;
    }
    task->multiskill = (NPAR)c;
    
    
    // dat_obs
//    StripedArray<NPAR> *striped_dat_obs = new StripedArray<NPAR>(fid, task->N);
//    task->dat_obs = striped_dat_obs->toArray();
    task->dat_obs = StripedArray<NPAR>::fromFileToArray(fid, task->N);
    task->dat_obs_stacked = init1D<NPAR>(task->Nst);
    task->dat_predict = initToValue2D<NUMBER>(task->N, task->nO, 0);
//    delete striped_dat_obs;

    
//    if(v==1) { // older NCAT of unsigned short
//        StripedArray<short> *dat_group_legacy = NULL;
//        StripedArray<short> *dat_skill_legacy = NULL;
//        StripedArray<short*> *dat_multiskill_legacy = NULL;
//        dat_group_legacy = new StripedArray<short>(fid, task->N);
//        for(t=0; t<task->N; t++)
//            task->dat_group[t] = (NCAT)dat_group_legacy->get(t);
//        if(task->multiskill == 0) {
//            dat_skill_legacy = new StripedArray<short>(fid, task->N);
//            for(t=0; t<task->N; t++)
//                task->dat_skill[t] = (NCAT)dat_skill_legacy->get(t);
//            delete dat_skill_legacy;
//        } else {
//            task->dat_multiskill = new StripedArray<NCAT*>(task->N,true);
//            nread = readMultiSkill(fid, param, v);
//            delete dat_multiskill_legacy;
//        }
//        delete dat_group_legacy;
//    } else {
        // dat_group
//        StripedArray<NCAT> *striped_dat_group = new StripedArray<NCAT>(fid, task->N);
//        task->dat_group = striped_dat_group->toArray();
//        delete striped_dat_group;
        task->dat_group = StripedArray<NCAT>::fromFileToArray(fid, task->N);
//    }
    
    // dat_item
//    StripedArray<NCAT> *striped_dat_item = new StripedArray<NCAT>(fid, task->N);
//    task->dat_item = striped_dat_item->toArray();
//    delete striped_dat_item;
    task->dat_item = StripedArray<NCAT>::fromFileToArray(fid, task->N);
    
    // dat_skill
    if(task->multiskill == 0) {
//        StripedArray<NCAT> *striped_dat_skill = new StripedArray<NCAT>(fid, task->N);
//        task->dat_skill = striped_dat_skill->toArray();
//        delete striped_dat_skill;
        task->dat_skill = StripedArray<NCAT>::fromFileToArray(fid, task->N);
    } else {
//            task->dat_multiskill = new StripedArray< NCAT* >(task->N,true);
//            nread = readMultiSkill(fid, param, v);
//        StripedArray<NCAT> *striped_dat_skill_stacked = new StripedArray<NCAT>(fid, task->Nst);
//        task->dat_skill_stacked = striped_dat_skill_stacked->toArray();
//        StripedArray<NCAT> *striped_dat_skill_rcount = new StripedArray<NCAT>(fid, task->N);
//        task->dat_skill_rcount = striped_dat_skill_rcount->toArray();
//        StripedArray<NCAT> *striped_dat_skill_rix = new StripedArray<NCAT>(fid, task->N);
//        task->dat_skill_rix = striped_dat_skill_rix->toArray();
//        delete striped_dat_skill_stacked;
//        delete striped_dat_skill_rcount;
//        delete striped_dat_skill_rix;
        task->dat_skill_stacked = StripedArray<NCAT>::fromFileToArray(fid, task->Nst);
        task->dat_skill_rcount = StripedArray<NCAT>::fromFileToArray(fid, task->N);
        task->dat_skill_rix = StripedArray<NCAT>::fromFileToArray(fid, task->N);
    }
    // dat_slices, only of nZ > 1
    if(task->nZ > 1) {
        NDAT szZ = (task->multiskill == 0)?task->N:task->Nst;
//        StripedArray<NPAR> *striped_dat_slice = new StripedArray<NPAR>(fid, szZ);
//        task->dat_slice = striped_dat_slice->toArray();
//        delete striped_dat_slice;
        task->dat_slice = StripedArray<NPAR>::fromFileToArray(fid, szZ);
    }
        
    string str;
    task->map_group_fwd = new map<string,NCAT>();
    task->map_group_bwd = new map<NCAT,string>();
    task->map_skill_fwd = new map<string,NCAT>();
    task->map_skill_bwd = new map<NCAT,string>();
    task->map_step_fwd = new map<string,NCAT>();
    task->map_step_bwd = new map<NCAT,string>();
    // voc_group
    for(NCAT g=0; g<task->nG; g++) {
        str = readString(fid);
        task->map_group_fwd->insert(pair<string,NCAT>(str, g));
        task->map_group_bwd->insert(pair<NCAT,string>(g, str));
    }
    
    // voc_skill
    for(NCAT k=0; k<task->nK; k++) {
        str = readString(fid);
        task->map_skill_fwd->insert(pair<string,NCAT>(str, k));
        task->map_skill_bwd->insert(pair<NCAT,string>(k, str));
    }
    
    // voc_item
    for(NCAT i=0; i<task->nI; i++) {
        str = readString(fid);
        task->map_step_fwd->insert(pair<string,NCAT>(str, i));
        task->map_step_bwd->insert(pair<NCAT,string>(i, str));
    }
    
    // build dat_obs_stacked
    NDAT co = 0;
    if(task->multiskill!=0) {
        for(NDAT t=0;t<task->N;t++) {
            for(NDAT tt=0; tt<task->dat_skill_rcount[t]; tt++) {
                task->dat_obs_stacked[c] = task->dat_obs[t];
                co++;
            }
        }
    } else {
        memcpy( task->dat_obs_stacked, task->dat_obs , sizeof(NPAR)*(size_t)task->N );
    }

    //
    // create connectivities: skill-skill, student-skill
    //
    createConnectivities(task);

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

bool InputUtilSt::toBin(struct task *task, const char *fn) {
    char c;
    NDAT i;
    FILE *fid = fopen(fn,"wb");

    // version
    c = bin_input_file_verstion;
    fwrite (&c , sizeof(char), 1, fid);
    
    // N
    i = task->N;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // Nst
    i = task->Nst;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // N_null
    i = task->N_null;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nO
    i = task->nO;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nG
    i = task->nG;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nI
    i = task->nI;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nK
    i = task->nK;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // nZ
    i = task->nZ;
    fwrite (&i , sizeof(NDAT), 1, fid);
    
    // multiskill
    c = task->multiskill;
    fwrite (&c , sizeof(char), 1, fid);
    
//    NDAT nwrit;
    // dat_obs
    /*nwrit = */StripedArray<NPAR>::arrayToBinFile(task->dat_obs, task->N, fid);//task->dat_obs->toBinFile(fid);
    // dat_group
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_group, task->N, fid);//->toBinFile(fid);
    // dat_item
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_item, task->N, fid);//->toBinFile(fid);
    // dat_skill
    if(task->multiskill == 0)
    /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_skill, task->N, fid);//->toBinFile(fid);
    else {
        //        /*nwrit = */writeMultiSkill(fid, param);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_skill_stacked, task->Nst, fid);//->toBinFile(fid);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_skill_rcount , task->N,        fid);//->toBinFile(fid);
        /*nwrit = */StripedArray<NCAT>::arrayToBinFile(task->dat_skill_rix    , task->N,        fid);//->toBinFile(fid);
    }
    // dat_slices, only of nZ > 1
    if(task->nZ > 1) {
        NDAT szZ = (task->multiskill == 0)?task->N:task->Nst;
        StripedArray<NPAR>::arrayToBinFile(task->dat_slice, szZ, fid);
    }

    map<NCAT,string>::iterator it;
    // voc_group
    for (it = task->map_group_bwd->begin(); it != task->map_group_bwd->end(); ++it) {
        writeString(fid, it->second);
    }

    // voc_skill
    for (it = task->map_skill_bwd->begin(); it != task->map_skill_bwd->end(); ++it) {
        writeString(fid, it->second);
    }
    
    // voc_item
    for (it =  task->map_step_bwd->begin(); it != task->map_step_bwd->end(); ++it) {
        writeString(fid, it->second);
    }
    
    
    fclose(fid);
    return true;
}

// experimental
// for a skill, write all student sequences as a matrix
// Nstudents * Max attempts, 1 - correct, 2 - incorrect, 0 - empty
// space separated coumns
//void InputUtilSt::writeInputMatrix(const char *filename, struct task* task, NCAT xndat, struct data** x_data) {
////    FILE *fid = fopen(filename,"w");
////    if(fid == NULL) {
////        fprintf(stderr,"Can't write output model file %s\n",filename);
////        exit(1);
////    }
//
//    std::ofstream file;
//    file.open(filename);
//
//    NDAT nmax = 0;
//    for(NCAT x=0; x<xndat; x++)
//        if(nmax < x_data[x]->n)
//            nmax = x_data[x]->n;
//
//    for(NCAT x=0; x<xndat; x++) {
//        std::stringstream ss;
//        for(NDAT t=0; t<nmax; t++) {
//
//            if(t<x_data[x]->n) {
//                ss << ((t>0)?" ":"") << (int)(1+p->dat_obs[ x_data[x]->ix[t] ]);
//            } else {
//                ss << " " << 0;
//            }
//        }
//        ss << "\n";
//        file << ss.str();
//    } // for all groups in skill
//
////    fclose(fid);
//    file.close();
//}
