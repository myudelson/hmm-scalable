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


// This parameter union uses sum in the logistic space as a union function.
// Probabilities are converted into logits ( m=ln(p/(1-p)) ), summed ( sum_m )
// and converted back to probability using logistic function ( p = 1 / (1 + e^{sum_m}) )

#include "utils.h"
#include "PULogistic.h"


// a shortcut for uniting one pooled one standatd parameter
NUMBER PULogistic::pair(NUMBER standard, NUMBER pooled) {
	NUMBER S = safe01num(standard);
	NUMBER P = safe01num(pooled);
	return 1/( 1 + (1-S)*(1-P)/(S*P) );
}

NUMBER PULogistic::unite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled) {
	NUMBER NUM = 1;
	NUMBER DEN = 1;
	if(standard != NULL && nstandard > 0) {
		for(NCAT i=0; i<=nstandard; i++) {
			NUM *= safe01num(1-standard[i]);
			DEN *= safe01num(  standard[i]);
		}
	}
	if(pooled != NULL && npooled > 0) {
		for(NCAT i=0; i<=nstandard; i++) {
			NUM *= safe01num(1-pooled[i]);
			DEN *= safe01num(  pooled[i]);
		}
	}
	return 1/( 1 + NUM/DEN );
}


// a shortcut for uniting one pooled one standatd parameter
NUMBER PULogistic::derivativePair(NUMBER standard, NUMBER pooled, NPAR which) {
	NUMBER united = pair(standard, pooled);
	// which 0-standard, 1-pooled
	NUMBER deriv_param = ((which==0)?
						  1 / safe0num( standard * (1-standard) ) :
						  1 / safe0num(   pooled * (1-  pooled) )
						  );
	return united * (1-united) * deriv_param;
}


NUMBER PULogistic::derivativeUnite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled, NPAR which, NCAT index) {
	NUMBER united = unite(standard, nstandard, pooled, npooled);
	// which 0-standard, 1-pooled
	NUMBER deriv_param = ((which==0)?
						  1 / safe0num( standard[index] * (1-standard[index]) ) :
						  1 / safe0num(   pooled[index] * (1-  pooled[index]) )
						  );
	return united * (1-united) * deriv_param;
}


NUMBER PULogistic::derivativeUnite(NUMBER united, NUMBER standard, NUMBER pooled, NPAR which) {
	// which 0-standard, 1-pooled
	NUMBER deriv_param = ((which==0)?
						  1 / safe0num( standard * (1-standard) ) :
						  1 / safe0num(   pooled * (1-  pooled) )
						  );
	return united * (1-united) * deriv_param;
}

