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

//
//  InputUtil.h
//  This is a support class that helps read input text for HMM, convert it
//  to a more compact binary file, and read the binary file as well
//  HMM

#ifndef __HMM__InputUtil__
#define __HMM__InputUtil__

#include "utils.h"

//#define bin_input_file_verstion 1
//#define bin_input_file_verstion 2 // increase number of skills/students to a 4 byte integer
#define bin_input_file_verstion 3 // added Nstacked, changed how multi-skills are stored and added slices (single and multi-coded)

class InputUtil {
public:
    static bool readTxt(const char *fn, struct param * param); // read txt into param
    static bool readBin(const char *fn, struct param * param); // read bin into param
    static bool toBin(struct param * param, const char *fn);// writes data in param to bin file
    // experimental
    static void writeInputMatrix(const char *filename, struct param* p, NCAT xndat, struct data** x_data);
private:
    static void writeString(FILE *f, string str);
    static string readString(FILE *f);
//    static NDAT writeMultiSkill(FILE *f, struct param * param);
//    static NDAT  readMultiSkill(FILE *f, struct param * param, char version);
};
#endif /* defined(__HMM__InputUtil__) */
