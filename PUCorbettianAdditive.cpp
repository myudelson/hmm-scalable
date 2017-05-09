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


// This parameter union uses the approach described in Corbett & Adnerson 1996
// Bayesian Knowledge Tracing paper. Instead of multiplicative introduction of
// pooled effects, it uses additition.
// Probabilities for standard parameters are converted into odds. Odds and raw
// pooled parameters are added. The sum is then covnerted into probability space.

#include <math.h>
#include "utils.h"
#include "PUCorbettianAdditive.h"


// a shortcut for uniting one pooled one standatd parameter
NUMBER PUCorbettianAdditive::pair(NUMBER standard, NUMBER pooled) {
	NUMBER sum_odds =  (standard / safe0num(1-standard)) + pooled;
	return sum_odds / (1+sum_odds);
}

NUMBER PUCorbettianAdditive::sumOdds(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled) {
	NUMBER sum_odds = 0;
	if(standard != NULL && nstandard > 0) {
		for(NCAT i=0; i<=nstandard; i++) {
			sum_odds += standard[i] / safe0num(1-standard[i]);
		}
	}
	if(pooled != NULL && npooled > 0) {
		for(NCAT i=0; i<=nstandard; i++) {
			sum_odds += pooled[i];
		}
	}
	return sum_odds;
}

NUMBER PUCorbettianAdditive::unite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled) {
	// W_{ x_i ,g_j}^{ add }=O_{ x_i ,g_j}^{ add }/( 1+O_{ x_i ,g_j}^{ add })
	NUMBER sum_odds = sumOdds(standard, nstandard, pooled, npooled);
	return sum_odds / (1+sum_odds);
}


// a shortcut for uniting one pooled one standatd parameter
NUMBER PUCorbettianAdditive::derivativePair(NUMBER standard, NUMBER pooled, NPAR which) {
	NUMBER paired = pair(standard, pooled);
	// which 0-standard, 1-pooled
	return 1 / ( pow(1+paired,2)) * ((which==0)?( 1/(pow( safe0num(1-standard) ,2)) ):1);
}


NUMBER PUCorbettianAdditive::derivativeUnite(NUMBER* standard, NCAT nstandard, NUMBER* pooled, NCAT npooled, NPAR which, NCAT index) {
	NUMBER sum_odds = sumOdds(standard, nstandard, pooled, npooled);
	// which 0-standard
	// \frac{1}{ (1+O_{ x_{ i },g_{ j } }^{ add } )^2 } \frac { 1 }{ (1-x_{ i })^2 }
	// which 1-pooled
	// \frac{1}{ (1+O_{ x_{ i },g_{ j } }^{ add } )2 }
	return 1 / ( pow(1+sum_odds,2) ) * ((which==0) ? ( 1/(pow( safe0num(1-standard[index]) ,2)) ) : 1 );
}

NUMBER PUCorbettianAdditive::derivativeUnite(NUMBER united, NUMBER standard, NUMBER pooled, NPAR which) {
	NUMBER sum_odds = united / safe0num(1 - united);
	// which 0-standard, 1-pooled
	return 1 / ( pow(1+sum_odds,2) ) * ((which==0) ? ( 1/(pow( safe0num(1-standard) ,2)) ) : 1 );
}

