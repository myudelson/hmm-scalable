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


// Parameter union is a helper class that "squishes" multiple probabilistic
// parameters together to produce a single probabilistic "superposition" value.
// We make distinction between standard parameters and "pooled" parameters.
// Pooled parameters denote phenomenae sampled from a population (like subjects).

#include "utils.h"

#ifndef _PARAMETERUNION_H
#define _PARAMETERUNION_H


class ParameterUnion {
public:
	// unions
	virtual NUMBER pair(NUMBER standard, NUMBER pooled) = 0; // a shortcut for uniting one pooled one standatd parameter
	virtual NUMBER unite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled) = 0;
	// derivatives of unions
	virtual NUMBER derivativePair(NUMBER standard, NUMBER pooled, NPAR which) = 0; // a shortcut for uniting one pooled one standatd parameter
	virtual NUMBER derivativeUnite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled, NPAR which, NCAT index) = 0;
};

#endif
