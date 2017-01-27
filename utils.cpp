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

#include "utils.h"
using namespace std;

// project of others
int compareNumber (const void * a, const void * b) {
	NUMBER cmp = ( *(NUMBER*)a - *(NUMBER*)b );
	return -1*(cmp<0) + 0 + (cmp>0)*1;
}

int compareNumberRev (const void * a, const void * b) {
    //	return ( *(NUMBER*)b - *(NUMBER*)a );
	NUMBER cmp = ( *(NUMBER*)b - *(NUMBER*)a );
	return -1*(cmp<0) + 0 + (cmp>0)*1;
}

void qsortNumber(NUMBER* ar, NPAR size) {
	qsort (ar, (size_t)size, sizeof(NUMBER), compareNumber);
}

void qsortNumberRev(NUMBER* ar, NPAR size) {
	qsort (ar, (size_t)size, sizeof(NUMBER), compareNumberRev);
}

int compareNcat (const void * a, const void * b) {
	int cmp = ( *(NCAT*)a - *(NCAT*)b );
	return -1*(cmp<0) + 0 + (cmp>0)*1;
}


void qsortNcat(NCAT* ar, NPAR size) {
	qsort (ar, (size_t)size, sizeof(NCAT), compareNcat);
}

// refer to http://arxiv.org/abs/1101.6081 for source
void projsimplex(NUMBER* y, NPAR size) {
    
	bool bget = false;
	NUMBER *s = init1D<NUMBER>((NDAT)size);
	cpy1D<NUMBER>(y, s, (NDAT)size);
	
	qsortNumberRev(s, size);
	NUMBER tmpsum = 0, tmax = 0;
    
	for(NPAR i=0; i<(size-1); i++) {
		tmpsum = tmpsum + s[i];
		tmax = (tmpsum - 1)/(i+1);
		if(tmax >= s[i+1]) {
			bget = true;
			break;
		}
	}
	if(!bget) tmax = (tmpsum + s[size-1] -1)/size;
	free(s);
    
	for(NPAR i=0; i<size; i++)
		y[i] = ((y[i]-tmax)<0)?0:(y[i]-tmax);
}


// project of my own
bool issimplex(NUMBER* ar, NPAR size) {
	NUMBER sum = 1;
	for(NPAR i=0; i<size; i++) {
		sum -= ar[i];
		if( ar[i] < 0 || ar[i] > 1)
			return false;
	}
//    for(NPAR i=0; i<size && fabs(sum)<SAFETY && fabs(sum)>0; i++) {
//        if( (ar[i]+sum) >=0 ) {
//            ar[i]+=sum;
//            break;
//        }
//        if( (ar[i]+sum) <=1 ) {
//            ar[i]+=sum;
//            break;
//        }
//    }
	return fabs(sum)<SAFETY;
}

bool issimplexbounded(NUMBER* ar, NUMBER *lb, NUMBER *ub, NPAR size) {
	NUMBER sum = 1;
	for(NPAR i=0; i<size; i++) {
		sum -= ar[i];
		if( ar[i] < lb[i] || ar[i] > ub[i])
			return false;
	}
//    for(NPAR i=0; i<size && fabs(sum)<SAFETY && fabs(sum)>0; i++) {
//        if( (ar[i]+sum) >=lb[i] ) {
//            ar[i]+=sum;
//            break;
//        }
//        if( (ar[i]+sum) <=ub[i] ) {
//            ar[i]+=sum;
//            break;
//        }
//    }
	return fabs(sum)<SAFETY;
}


void projectsimplexOld(NUMBER* ar, NPAR size) {
	NPAR i, num_at_hi, num_at_lo; // number of elements at lower,upper boundary
	NPAR *at_hi = Calloc(NPAR, (size_t)size);
	NPAR *at_lo = Calloc(NPAR, (size_t)size);
	NUMBER err, lambda;
	while( !issimplex(ar, size)) {
        lambda = 0;
		num_at_hi = 0;
		num_at_lo = 0;
		err = -1; // so that of the sum is over 1, the error is 1-sum
		// threshold
		for(i=0; i<size; i++) {
			at_lo[i] = (ar[i]<0)?1:0;
			ar[i] = (at_lo[i]==1)?0:ar[i];
			num_at_lo = (NPAR)(num_at_lo + at_lo[i]);
            
			at_hi[i] = (ar[i]>1)?1:0;
			ar[i] = (at_hi[i]==1)?1:ar[i];
			num_at_hi = (NPAR)(num_at_hi + at_hi[i]);
			
			err += ar[i];
		}
		if (size > (num_at_hi + num_at_lo) )
			lambda = err / (size - (num_at_hi + num_at_lo));
		for(i=0; i<size; i++)
			ar[i] -= (at_lo[i]==0 && at_hi[i]==0)?lambda:0;
		
	} // until satisfied
    // force last to be 1 - sum of all but last
    err = 1;
    for(i=0; i<(size-1); i++) err -= ar[i];
    ar[size-1] = err;
    err = 1;
    for(i=1; i<(size-0); i++) err -= ar[i];
    ar[0] = err;

	free(at_hi);
	free(at_lo);
}

void projectsimplex(NUMBER* ar, NPAR size) {
    NPAR i, num_at_hi, num_at_lo; // number of elements at lower,upper boundary
    NPAR *at_hi = Calloc(NPAR, (size_t)size);
    NPAR *at_lo = Calloc(NPAR, (size_t)size);
    NUMBER err, lambda;
    NUMBER* ar_copy = Calloc(NUMBER, (size_t)size);
    memcpy(ar_copy, ar, (size_t)size);
    int iter = 0;
    while( !issimplex(ar, size) ) {
        lambda = 0;
        num_at_hi = 0;
        num_at_lo = 0;
        err = -1;
        // threshold
        for(i=0; i<size; i++) {
            at_lo[i] = (ar[i]<=/*lb[i]*/0)?1:0;
            ar[i] = (at_lo[i]==1)?/*lb[i]*/0:ar[i];
            num_at_lo = (NPAR)( num_at_lo + at_lo[i] );
            
            at_hi[i] = (ar[i]>=/*ub[i]*/1)?1:0;
            ar[i] = (at_hi[i]==1)?/*ub[i]*/1:ar[i];
            num_at_hi = (NPAR)( num_at_hi + at_hi[i] );
            
            err += ar[i];
        }
        //		if (size > (num_at_hi + num_at_lo) )
        //			lambda = err / (size - (num_at_hi + num_at_lo));
        if ( err > 0 && size > num_at_lo )
            lambda = err / (size - num_at_lo);
        else if ( err < 0 && size > num_at_hi )
            lambda = err / (size - num_at_hi);
        
        int will_suffer_from_lambda = 0; // those values that, if lessened by lambda, will be below 0 or over 1
        NUMBER err2 = 0.0, lambda2 = 0.0;
        for(i=0; i<size; i++) {
            if( (at_lo[i]==0 && at_hi[i]==0) ) {
                if( ar[i] < lambda ) {
                    will_suffer_from_lambda++;
                    err2 += ar[i];
                }
            }
        }
        
        if( will_suffer_from_lambda == 0 ) {
            for(i=0; i<size; i++) {
                ar[i] -= (at_lo[i]==0 && err>0)?lambda:0;
                ar[i] -= (at_hi[i]==0 && err<0)?lambda:0;
            }
        } else {
            lambda2 = (err-err2) / ( size - (num_at_hi + num_at_lo) - will_suffer_from_lambda );
            for(i=0; i<size; i++) {
                if( at_lo[i]==0 && at_hi[i]==0 ) {
                    ar[i] = (ar[i]<lambda)?0:(ar[i]-lambda2);
                }
            }
        }
        iter++;
        if(iter==100) {
            fprintf(stderr,"WARNING! Stuck in projectsimplex().\n");
            exit(1);
        }
    } // until satisfied
    // force last to be 1 - sum of all but last -- this code, actually breaks things
    //    err = 0;
    //    for(i=0; i<(size); i++) {
    //        err += ar[i];
    //    }
    //    err = 1;
    //    for(i=0; i<(size-1); i++) {
    //        err -= ar[i];
    //    }
    //    ar[size-1] = err;
    //    err = 1;
    //    for(i=1; i<(size-0); i++) {
    //        err -= ar[i];
    //    }
    //    ar[0] = err;
    
    NUMBER sum = 0.0;
    for(i=0; i<size; i++) {
        sum += ar[i];
        if(ar[i]<0 || ar[i] >1)
            fprintf(stderr, "ERROR! projected value is not within [0, 1] range!\n");
    }
//    if( fabs(sum-1)>SAFETY)
//        fprintf(stderr, "ERROR! projected simplex does not sum to 1!\n");
	
    free(at_hi);
    free(at_lo);
}

void projectsimplexbounded(NUMBER* ar, NUMBER *lb, NUMBER *ub, NPAR size) {
	NPAR i, num_at_hi, num_at_lo; // number of elements at lower,upper boundary
	NPAR *at_hi = Calloc(NPAR, (size_t)size);
	NPAR *at_lo = Calloc(NPAR, (size_t)size);
	NUMBER err, lambda, v;
    NUMBER* ar_copy = Calloc(NUMBER, (size_t)size);
    memcpy(ar_copy, ar, (size_t)size);
    int iter = 0;
    for(i=0; i<size; i++)
        if(ar[i]!=ar[i]) {
            fprintf(stderr,"WARNING! NaN detected!\n");
        }
    bool doexit = false;
	while( !issimplexbounded(ar, lb, ub, size) && !doexit ) {
        lambda = 0;
		num_at_hi = 0;
		num_at_lo = 0;
		err = -1;
		// threshold
		for(i=0; i<size; i++) {
			at_lo[i] = (ar[i]<=lb[i])?1:0;
            
            v = (at_lo[i]==1)?lb[i]:ar[i];
            if(v!=v) {
                fprintf(stderr,"WARNING! NaN to be set!\n");
            }
			
            ar[i] = (at_lo[i]==1)?lb[i]:ar[i];
			num_at_lo = (NPAR)( num_at_lo + at_lo[i] );
			
			at_hi[i] = (ar[i]>=ub[i])?1:0;
            
            v = (at_hi[i]==1)?ub[i]:ar[i];
            if(v!=v) {
                fprintf(stderr,"WARNING! NaN to be set!\n");
            }
            
			ar[i] = (at_hi[i]==1)?ub[i]:ar[i];
			num_at_hi = (NPAR)( num_at_hi + at_hi[i] );
			
			err += ar[i];
		}
        if ( err > 0 && size > num_at_lo )
            lambda = err / (size - num_at_lo);
        else if ( err < 0 && size > num_at_hi )
            lambda = err / (size - num_at_hi);

        int will_suffer_from_lambda = 0; // those values that, if lessened by lambda, will be below 0 or over 1
        NUMBER err2 = 0.0, lambda2 = 0.0;
        for(i=0; i<size; i++) {
            if( (at_lo[i]==0 && at_hi[i]==0) ) {
                if( ar[i] < lambda ) {
                    will_suffer_from_lambda++;
                    err2 += ar[i];
                }
            }
        }
        
        if( will_suffer_from_lambda == 0 ) {
            for(i=0; i<size; i++) {

                v = ar[i] - ((at_lo[i]==0 && err>0)?lambda:0);
                if(v!=v) {
                    fprintf(stderr,"WARNING! NaN to be set!\n");
                }

                ar[i] -= (at_lo[i]==0 && err>0)?lambda:0;

                v = ar[i] - ((at_hi[i]==0 && err<0)?lambda:0);
                if(v!=v) {
                    fprintf(stderr,"WARNING! NaN to be set!\n");
                }

                ar[i] -= (at_hi[i]==0 && err<0)?lambda:0;
            }
        } else {
            lambda2 = (err-err2) / ( size - (num_at_hi + num_at_lo) - will_suffer_from_lambda );
            for(i=0; i<size; i++) {
                if( at_lo[i]==0 && at_hi[i]==0 ) {

                    v = (ar[i]<lambda)?0:(ar[i]-lambda2);
                    if(v!=v) {
                        fprintf(stderr,"WARNING! NaN to be set!\n");
                    }
                    
                    ar[i] = (ar[i]<lambda)?0:(ar[i]-lambda2);
                }
            }
        }
        iter++;
        if(iter==100) {
            fprintf(stderr,"WARNING! Stuck in projectsimplexbounded().\n");
//            doexit = true;
            exit(1);
        }
	} // until satisfied
    // force last to be 1 - sum of all but last -- this code, actually breaks things
//    err = 0;
//    for(i=0; i<(size); i++) {
//        err += ar[i];
//    }
//    err = 1;
//    for(i=0; i<(size-1); i++) {
//        err -= ar[i];
//    }
//    ar[size-1] = err;
//    err = 1;
//    for(i=1; i<(size-0); i++) {
//        err -= ar[i];
//    }
//    ar[0] = err;
    
    NUMBER sum = 0.0;
    for(i=0; i<size; i++) {
        sum += ar[i];
        if(ar[i]<0 || ar[i] >1)
            fprintf(stderr, "ERROR! projected value is not within [0, 1] range!\n");
    }
//    if( fabs(sum-1)>SAFETY) {
//        fprintf(stderr, "ERROR! projected simplex does not sum to 1!\n");
//    }
	
	free(at_hi);
	free(at_lo);
}


// protect against hard 0 or hard 1
NUMBER safe01num(NUMBER val) {
    //    val = (val<0)?0:((val>1)?1:val); // squeeze into [0,1]
    //	return val + SAFETY*(val==0) - SAFETY*(val==1); // then futher in
//	return (val<=0)? SAFETY : ( (val>=1)? (1-SAFETY) : val );
	return (val<=SAFETY)? SAFETY : ( (val>=(1-SAFETY))? (1-SAFETY) : val );
}

// protect against hard 0
NUMBER safe0num(NUMBER val) {
    //    return (fabs(val)<SAFETY)?(SAFETY*(val>=0) + SAFETY*(val<0)*(-1)):val;
//    NUMBER a_sign = (val<0)?-1:1;
//    return (fabs(val)<SAFETY)?a_sign*SAFETY:val;
	return (val<=SAFETY)? SAFETY : val;
}

NUMBER itself(NUMBER val) {
	return val;
}

NUMBER deprecated_fsafelog(NUMBER val) {
	return safelog(val + (val<=0)*SAFETY);
}

NUMBER safelog(NUMBER val) {
	return log(val + (val<=0)*SAFETY);
}

NUMBER sigmoid(NUMBER val) {
    return 1 / (1 + exp(-val));
}

NUMBER logit(NUMBER val) {
    NUMBER prob = (val<=0)? SAFETY : ( (val>=1)? (1-SAFETY) : val );
    return log( prob / (1-prob) );
    //    return fsafelog(val/safe0num(1-val));
}

// computes sigmoid( logit(p) + logit(q) )
NUMBER pairing(NUMBER p, NUMBER q) {
    NUMBER P = safe01num(p);
    NUMBER Q = safe01num(q);
    return 1/( 1 + (1-P)*(1-Q)/(P*Q) );
}

// computes sigmoid( logit(p1) + logit(p2) + ... )
//NUMBER squishing(NUMBER* p, NCAT n) {
//	NUMBER* P = Calloc(NUMBER, (size_t)n);
//	for(NCAT i=0; i<n; i++) P[i] = safe01num(p[i]);
//	NUMBER num = 1, den = 1;
//	for(NCAT i=0; i<n; i++) {
//		num *= (1-P[i]);
//		den *= P[i];
//	}
//	free(P);
//	return 1/( 1 + num/safe0num(den) );
//}

// max value of n
NUMBER maxn(NUMBER *ar, NDAT n) {
	NUMBER mx = ar[0];
	for(NDAT i=1; i<n; i++)
		if(ar[i]>mx) mx = ar[i];
	return mx;
}

//#define logit(y)
//#define logit(y)
//NUMBER logit(NUMBER val) {
//    //	return fsafelog( val / safe0num(1-val) );
//    return fastsafelog((val+(val==0)*SAFETY)/(1-val+(val>=1)*SAFETY));
//}

NUMBER sgn(NUMBER val) {
	return (0 < val) - (val < 0);
}

void add1DNumbersWeighted(NUMBER* sourse, NUMBER* target, NPAR size, NUMBER weight) {
	for(NPAR i=0; i<size; i++)
		target[i] = target[i] + sourse[i]*weight;
}

void add2DNumbersWeighted(NUMBER** sourse, NUMBER** target, NPAR size1, NPAR size2, NUMBER weight) {
    for(NPAR i=0; i<size1; i++)
        for(NPAR j=0; j<size2; j++)
            target[i][j] = target[i][j] + sourse[i][j]*weight;
}

void add3DNumbersWeighted(NUMBER*** sourse, NUMBER*** target, NPAR size1, NPAR size2, NPAR size3, NUMBER weight) {
    for(NPAR i=0; i<size1; i++)
        for(NPAR j=0; j<size2; j++)
            for(NPAR l=0; l<size3; l++)
                target[i][j][l] = target[i][j][l] + sourse[i][j][l]*weight;
}

bool isPasses(NUMBER* ar, NPAR size) {
	NUMBER sum = 0;
	for(NPAR i=0; i<size; i++) {
		if( ar[i]<0 || ar[i]>1)
			return false;
		sum += ar[i];
	}
	return sum==1;
}

bool isPassesLim(NUMBER* ar, NPAR size, NUMBER *lb, NUMBER* ub) {
	NUMBER sum = 0;
	for(NPAR i=0; i<size; i++) {
		if( (lb[i]-ar[i])>SAFETY || (ar[i]-ub[i])>SAFETY )
			return false;
		sum += ar[i];
	}
	return fabs(sum-1)<SAFETY;
}

// scale by smallest factor of 10 (max scaling by default
NUMBER doLog10Scale1D(NUMBER *ar, NPAR size) {
	NPAR i;
	NUMBER min_10_scale = 1000, max_10_scale = 0, candidate;
	for(i=0; i<size; i++) {
		if( fabs(ar[i]) < SAFETY ) // 0 gradient
			continue;
		candidate = floor( log10( fabs(ar[i]) ) );
		if(candidate < min_10_scale)
			min_10_scale = candidate;
		candidate = ceil( log10( fabs(ar[i]) ) );
		if(candidate > max_10_scale)
			max_10_scale = candidate;
	}
	min_10_scale++;
	max_10_scale++;
	if(max_10_scale > 0)
		for(i=0; i<size; i++)
			ar[i] = ar[i] / pow(10, max_10_scale);
    //	if(min_10_scale<1000)
    //		for(i=0; i<size; i++)
    //			ar[i] = ar[i] / pow(10, min_10_scale);
	return pow(10, max_10_scale);
}

// scale by smallest factor of 10  (max scaling by default
NUMBER doLog10Scale2D(NUMBER **ar, NPAR size1, NPAR size2) {
	NPAR i,j;
	NUMBER min_10_scale = 1000, max_10_scale = 0, candidate;
	for(i=0; i<size1; i++)
		for(j=0; j<size2; j++) {
			if( fabs(ar[i][j]) < SAFETY ) // 0 gradient
				continue;
			candidate = floor( log10( fabs(ar[i][j]) ) );
			if(candidate < min_10_scale)
				min_10_scale = candidate;
			candidate = ceil( log10( fabs(ar[i][j]) ) );
			if(candidate > max_10_scale)
				max_10_scale = candidate;
		}
	min_10_scale++;
	max_10_scale++;
	if(max_10_scale >0 )
		for(i=0; i<size1; i++)
			for(j=0; j<size2; j++)
				ar[i][j] = ar[i][j] / pow(10, max_10_scale);
    //	if(min_10_scale<1000)
    //		for(i=0; i<size1; i++)
    //			for(j=0; j<size2; j++)
    //				ar[i][j] = ar[i][j] / pow(10, min_10_scale);
	return pow(10, max_10_scale);
}


// Gentle - as per max distance to go toward extreme value of 0 or 1
NUMBER doLog10Scale1DGentle(NUMBER *grad, NUMBER *par, NPAR size) {
	NPAR i;
	NUMBER max_10_scale = 0, candidate, min_delta = 1, max_grad = 0, scale;
	for(i=0; i<size; i++) {
		if( fabs(grad[i]) < SAFETY ) // 0 gradient
			continue;
        // scale
        if(max_grad < fabs(grad[i]))
            max_grad = fabs(grad[i]);
		candidate = ceil( log10( fabs(grad[i]) ) );
		if(candidate > max_10_scale)
			max_10_scale = candidate;
        // delta: if grad<0 - distance to 1, if grad>0 - distance to 0
        candidate = (grad[i]<0)*(1-par[i]) + (grad[i]>0)*(par[i]) +
        ( (fabs(par[i])< SAFETY) || (fabs(1-par[i])< SAFETY) ); // these terms are there to avoid already extreme 0, 1 values;
        if( candidate < min_delta)
            min_delta = candidate;
	}
    max_grad = max_grad / pow(10, max_10_scale);
    
    scale = (max_grad!=0) ? ( 0.95 * min_delta / max_grad ) / pow(10, max_10_scale) : 1;
    if(scale==std::numeric_limits<double>::infinity() || scale==(-1*std::numeric_limits<double>::infinity()) ) {
        fprintf(stderr,"WARNING! scale is Infinite\n");
    }
    
    NUMBER v;
    if(max_10_scale > 0) {
        for(i=0; i<size; i++) {
            
            v = ( 0.95 * min_delta / max_grad) * grad[i] / pow(10, max_10_scale);
            if(v!=v) {
                fprintf(stderr,"WARNING! NaN to be set!\n");
            }
            
            grad[i] = ( 0.95 * min_delta / max_grad) * grad[i] / pow(10, max_10_scale);
        }
    }
    return scale; //( 0.95 * min_delta / max_grad ) / pow(10, max_10_scale);
}

// Gentle - as per max distance to go toward extreme value of 0 or 1
NUMBER doLog10Scale2DGentle(NUMBER **grad, NUMBER **par, NPAR size1, NPAR size2) {
	NPAR i,j;
	NUMBER max_10_scale = 0, candidate, min_delta = 1, max_grad = 0;
	for(i=0; i<size1; i++)
		for(j=0; j<size2; j++) {
			if( fabs(grad[i][j]) < SAFETY ) // 0 gradient
				continue;
            // scale
            if(max_grad < fabs(grad[i][j]))
                max_grad = fabs(grad[i][j]);
			candidate = ceil( log10( fabs(grad[i][j]) ) );
			if(candidate > max_10_scale)
				max_10_scale = candidate;
            // delta: if grad<0 - distance to 1, if grad>0 - distance to 0
            candidate = (grad[i][j]<0)*(1-par[i][j]) + (grad[i][j]>0)*(par[i][j]) +
            ( (fabs(par[i][j])< SAFETY) || (fabs(1-par[i][j])< SAFETY) ); // these terms are there to avoid already extreme 0, 1 values
            if( candidate < min_delta)
                min_delta = candidate;
		}
    max_grad = max_grad / pow(10, max_10_scale);
    if(max_10_scale >0 )
		for(i=0; i<size1; i++)
			for(j=0; j<size2; j++)
				grad[i][j] = ( 0.95 * min_delta / max_grad) * grad[i][j] / pow(10, max_10_scale);
	return ( 0.95 * min_delta / max_grad ) / pow(10, max_10_scale);
}

// Gentle - as per max distance to go toward extreme value of 0 or 1
NUMBER doLog10Scale3DGentle(NUMBER ***grad, NUMBER ***par, NPAR size1, NPAR size2, NPAR size3) {
    NPAR i,j,k;
    NUMBER max_10_scale = 0, candidate, min_delta = 1, max_grad = 0;
    for(i=0; i<size1; i++)
        for(j=0; j<size2; j++)
            for(k=0; k<size3; k++) {
                if( fabs(grad[i][j][k]) < SAFETY ) // 0 gradient
                    continue;
                // scale
                if(max_grad < fabs(grad[i][j][k]))
                    max_grad = fabs(grad[i][j][k]);
                candidate = ceil( log10( fabs(grad[i][j][k]) ) );
                if(candidate > max_10_scale)
                    max_10_scale = candidate;
                // delta: if grad<0 - distance to 1, if grad>0 - distance to 0
                candidate = (grad[i][j][k]<0)*(1-par[i][j][k]) + (grad[i][j][k]>0)*(par[i][j][k]) +
                ( (fabs(par[i][j][k])< SAFETY) || (fabs(1-par[i][j][k])< SAFETY) ); // these terms are there to avoid already extreme 0, 1 values
                if( candidate < min_delta)
                    min_delta = candidate;
            }
    max_grad = max_grad / pow(10, max_10_scale);
    if(max_10_scale >0 )
        for(i=0; i<size1; i++)
            for(j=0; j<size2; j++)
                for(k=0; k<size3; k++)
                    grad[i][j][k] = ( 0.95 * min_delta / max_grad) * grad[i][j][k] / pow(10, max_10_scale);
    return ( 0.95 * min_delta / max_grad ) / pow(10, max_10_scale);
}


// for skill of group
void zeroLabels(NCAT xdat, struct data** x_data) { // set counts in data sequences to zero
	NCAT x;
	for(x=0; x<xdat; x++)
		x_data[x][0].cnt = 0;
}

//
// The heavy end - common functionality
//

void set_param_defaults(struct param *param) {
    param->item_complexity = NULL;
	// configurable - set
    param->tol                   = 0.01;
    param->tol_mode              = 'p';
	param->scaled                = 0;
    param->do_not_check_constraints = 0;
	param->sliced                = 0;
    param->duplicate_console     = 0;
	param->maxiter               = 200;
	param->quiet                 = 0;
	param->single_skill          = 0;
	param->structure			 = 1; // default is by skill
	param->solver				 = 2; // default is Gradient Descent
	param->solver_setting		 = -1; // -1 - not set
    param->metrics               = 0;
    param->metrics_target_obs    = 0;
    param->predictions           = 0;
    param->update_known          = 'r';
    param->update_unknown        = 't';
    param->binaryinput           = 0;
	param->Cw                     = Calloc(NUMBER, (size_t)1);
    param->Cw[0]                  = 0;
    param->Ccenters              = NULL;
    param->initfile[0]           = 0; // 1st bit is 0 - no file
	param->init_params			 = Calloc(NUMBER, (size_t)5);
	param->init_params[0] = 0.5; // PI[0]
	param->init_params[1] = 1.0; // p(not forget)
	param->init_params[2] = 0.4; // p(learn)
	param->init_params[3] = 0.8; // p(not slip)
	param->init_params[4] = 0.2; // p(guess)
	param->param_lo				= Calloc(NUMBER, (size_t)10);
	param->param_lo[0] = 0; param->param_lo[1] = 0; param->param_lo[2] = 1; param->param_lo[3] = 0; param->param_lo[4] = 0;
	param->param_lo[5] = 0; param->param_lo[6] = 0; param->param_lo[7] = 0; param->param_lo[8] = 0; param->param_lo[9] = 0;
	param->param_hi				= Calloc(NUMBER, (size_t)10);
	param->param_hi[0] = 1.0; param->param_hi[1] = 1.0; param->param_hi[2] = 1.0; param->param_hi[3] = 0.0; param->param_hi[4] = 1.0;
	param->param_hi[5] = 1.0; param->param_hi[6] = 1.0; param->param_hi[7] = 0.3; param->param_hi[8] = 0.3; param->param_hi[9] = 1.0;
	param->cv_folds = 0;
	param->cv_strat = 'g'; // default group(student)-stratified
    param->cv_target_obs = 0; // 1st state to validate agains by default, cv_folds enables cross-validation
    param->cv_folds_file[0] = 0; // empty folds file
    param->cv_inout_flag = 'o'; // default rule, we're writing folds out
    param->multiskill = 0; // single skill per ovservation by default
    param->parallel = 0; // parallelization flag, no parallelization (0) by default
    // parse running settings
    param->init_reset = false; // init parameters specified
    param->lo_lims_specd = false; // parameter limits s`pecified
    param->hi_lims_specd = false; // parameter limits s`pecified
    param->stat_specd_gt2 = false; // number of states specified to be >2
    // vocabilaries
    param->map_group_fwd = NULL;
    param->map_group_bwd = NULL;
    param->map_step_fwd = NULL;
    param->map_step_bwd = NULL;
    param->map_skill_fwd = NULL;
    param->map_skill_bwd = NULL;
	// derived from data - set to 0
    param->N  = 0; //will be dynamically set in read_data_...()
    param->Nstacked  = 0; //will be dynamically set in read_data_...()
	param->nS = 2;
	param->nO = 0;
	param->nG = 0;
	param->nK = 0;
	param->nI = 0;
    param->nZ = 1;
	// data
    param->all_data = NULL;
    param->nSeq = 0;
	param->k_numg = NULL;
	param->k_data = NULL;
	param->k_g_data = NULL;
	param->g_numk = NULL;
	param->g_data = NULL;
	param->g_k_data = NULL;
    param->N_null = 0;
    param->n_null_skill_group = 0;
    param->null_skills = NULL;
	// fitting specific - Armijo rule
	param->ArmijoC1            = 1e-4;
	param->ArmijoC2            = 0.9; // since we're not newton method
	param->ArmijoReduceFactor  = 2;//1/0.9;//
	param->ArmijoSeed          = 1.0; //1; - since we use smooth stepping 1 is the only thing we need
    param->ArmijoMinStep       = 0.001; //  0.000001~20steps, 0.001~10steps
    // coord descend
    param->first_iteration_qualify = 0;
    param->iterations_to_qualify   = 2;
    // block fitting of some parameters
    param->block_fitting_type = 0; // no bocking of fitting - TODO, enable diff block types
    param->block_fitting[0] = 0; // no bocking fitting for PI
    param->block_fitting[1] = 0; // no bocking fitting for A
    param->block_fitting[2] = 0; // no bocking fitting for B
    param->per_kc_rmse_acc = false;
    // data
    param->dat_obs = NULL;
    param->dat_group = NULL;
    param->dat_skill = NULL;
    param->dat_skill_stacked = NULL;
    param->dat_skill_rcount = NULL;
    param->dat_skill_rix = NULL;
    param->dat_item = NULL;
//    param->dat_multiskill = NULL;
    param->dat_slice = NULL;
	param->tag1 = 0;
}

void destroy_input_data(struct param *param) {
    if(param->item_complexity) free(param->item_complexity);
	if(param->init_params != NULL) free(param->init_params);
	if(param->param_lo != NULL) free(param->param_lo);
	if(param->param_hi != NULL) free(param->param_hi);
	
    // data - checks if pointers to data are null anyway (whether we delete linear columns of data or not)
	if(param->dat_obs != NULL) free( param->dat_obs );
	if(param->dat_group != NULL) free( param->dat_group );
	if(param->dat_item != NULL) free( param->dat_item );
	if(param->dat_skill != NULL) free( param->dat_skill );
//	if(param->dat_multiskill != NULL) delete param->dat_multiskill;
    if(param->dat_skill_stacked != NULL) free( param->dat_skill_stacked );
    if(param->dat_skill_rcount != NULL) free( param->dat_skill_rcount );
    if(param->dat_skill_rix != NULL) free( param->dat_skill_rix );
	if(param->dat_slice != NULL) free( param->dat_slice );
    if( param->Cslices>0 ) {
        free( param->Cw );
        free( param->Ccenters );
    }
    
    // not null skills
    for(NDAT kg=0;kg<param->nSeq; kg++) {
		free(param->all_data[kg].ix); // was obs;
		if( param->all_data[kg].ix_stacked != NULL ) free(param->all_data[kg].ix_stacked); // was obs;
//        if(param->sliced) // handled via one global array and ix indexing
//            free(param->all_data[kg].time);
    }
    if(param->all_data != NULL) free(param->all_data); // ndat of them
    if(param->k_data != NULL)   free(param->k_data); // ndat of them (reordered by k)
    if(param->g_data != NULL)   free(param->g_data); // ndat of them (reordered by g)
    if(param->k_g_data != NULL) free(param->k_g_data); // nK of them
    if(param->g_k_data != NULL) free(param->g_k_data); // nG of them
    
	if(param->k_numg != NULL)   free(param->k_numg);
	if(param->g_numk != NULL)   free(param->g_numk);
    // null skills
    for(NCAT g=0;g<param->n_null_skill_group; g++)
        free(param->null_skills[g].ix); // was obs
    if(param->null_skills != NULL) free(param->null_skills);
    // vocabularies
    delete param->map_group_fwd;
    delete param->map_group_bwd;
    delete param->map_step_fwd;
    delete param->map_step_bwd;
    delete param->map_skill_fwd;
    delete param->map_skill_bwd;
}


//
// read/write solver info to a file
//
void writeSolverInfo(FILE *fid, struct param *param) {
    // solver id
    if( param->solver_setting>0 ) {
        fprintf(fid,"SolverId\t%d.%d.%d\n",param->structure,param->solver,param->solver_setting);
    } else {
        fprintf(fid,"SolverId\t%d.%d\n",param->structure,param->solver);
    }
	// nK
    fprintf(fid,"nK\t%d\n",param->nK);
    // nG
    fprintf(fid,"nG\t%d\n",param->nG);
    // nS
    fprintf(fid,"nS\t%d\n",param->nS);
    // nO
    fprintf(fid,"nO\t%d\n",param->nO);
    // nZ
    fprintf(fid,"nZ\t%d\n",param->nZ);
}

void readSolverInfo(FILE *fid, struct param *param, NDAT *line_no) {
    string s;
    int c, i1, i2;
    // SolverId
    fscanf(fid,"SolverId\t%i.%i", &i1, &i2);
    param->structure = (NPAR) i1;
    param->solver    = (NPAR) i2;
    fscanf(fid,"SolverId\t%hhu.%hhu", &param->structure, &param->solver);
    c = fscanf(fid,".%hhu\n", &param->solver_setting);
    if( c<1 ) {
        fscanf(fid,"\n");
        param->solver_setting = -1;
    }
    (*line_no)++;
	// nK
    fscanf(fid,"nK\t%i\n",&param->nK);
    (*line_no)++;
    // nG
    fscanf(fid,"nG\t%i\n",&param->nG);
    (*line_no)++;
    // nS
    fscanf(fid,"nS\t%hhu\n",&param->nS);
    (*line_no)++;
    // nO
    fscanf(fid,"nO\t%hhu\n",&param->nO);
    (*line_no)++;
    // nZ
    fscanf(fid,"nZ\t%hhu\n",&param->nZ);
    (*line_no)++;
}

//
// Handling blocking labels
//
void zeroLabels(struct param* param) { // set counts in data sequences to zero
	NCAT k,g;
	for(k=0; k<param->nK; k++)
		for(g=0; g<param->k_numg[k]; g++)
			param->k_g_data[k][g]->cnt = 0;
}

//
// clear up all forward/backward/etc variables for a skill-slice
//
void RecycleFitData(NCAT xndat, struct data** x_data, struct param *param) {
	NCAT x;
	NDAT t;
	for(x=0; x<xndat; x++) {
        //        if( x_data[x][0].cnt != 0)
        //            continue;
		if( x_data[x][0].alpha != NULL ) {
			free2D<NUMBER>(x_data[x][0].alpha, x_data[x][0].n);  // only free data here
			x_data[x][0].alpha = NULL;
		}
		if( x_data[x][0].c != NULL ) {
			free(x_data[x][0].c);  // only free data here
			x_data[x][0].c = NULL;
		}
		if( x_data[x][0].beta != NULL ) {
			free2D<NUMBER>(x_data[x][0].beta,  x_data[x][0].n); // only free data here
			x_data[x][0].beta = NULL;
		}
		if( x_data[x][0].gamma != NULL ) {
			free2D<NUMBER>(x_data[x][0].gamma, x_data[x][0].n); // only free data here
			x_data[x][0].gamma = NULL;
		}
		if( x_data[x][0].xi != NULL ) {
			for(t=0;t<x_data[x][0].n; t++)
				free2D<NUMBER>(x_data[x][0].xi[t],  (NDAT)param->nS); // only free data here
			x_data[x][0].xi = NULL;
		}
	}
}

//
// working with time
//

// limits are the borders of time bins, there are nlimits+1 bins total,  bins 0:nlimits
NPAR sec_to_linear_interval(int time, int *limits, NPAR nlimits){
    for(NPAR i=0; i<nlimits; i++)
        if( time < limits[i] )
            return i;
    return nlimits;
}

// limits are the borders of time bins, there are nlimits+1 bins total
int time_lim_20HDWM[5] = {20*60, 60*60, 24*60*60, 7*24*60*60, 30*24*60*60}; // 20min, hour, day, week, month

// 9 categories: <2m, <20m, <1h, same day, next day, same week, next week, <30d, >=30d
NPAR sec_to_9cat(int time1, int time2, int *limits, NPAR nlimits) {
    int diff = time2 - time1;
    
    if(diff <0 ) {
        fprintf(stderr,"ERROR! time 1 should be smaller than time 2\n");
        return 0;
    } else if ( diff < (2*60) ) { // 20min
        return 0;
    } else if ( diff < (20*60) ) { // 20min
        return 1;
    } else if (diff < (60*60) ) { // 1h
        return 2;
    } else if(diff < (14*24*60*60)) { // detect date structures
        time_t t1 = (time_t)time1;
        time_t t2 = (time_t)time2;
        struct tm * ttm1 = localtime (&t1);
        struct tm * ttm2 = localtime (&t2);
        if( diff < (24*60*60) && ttm1->tm_mday==ttm2->tm_mday ) { // same day
            return 3;
        } else if(    diff < (2*24*60*60) &&                      // next day
                  ( (ttm2->tm_wday == (ttm1->tm_wday+1)) ||
                   (ttm1->tm_wday==6 && ttm2->tm_wday==0)
                   )
                  ) {
            return 4;
        } else if( diff < (7*24*60*60) && (ttm1->tm_wday<ttm2->tm_wday) ) { // this week
            return 5;
        } else { // next week
            return 6;
        }
    } else if(diff < (30*24*60*60)) { // less than 30 days
        return 7;
    } else if(diff >= (30*24*60*60)) { // 30 days or more
        return 8;
    }
    return 0;
}


//// write time intervals to file
//void write_time_interval_data(param* param, const char *file_name) {
//    if(param->time != 1) {
//        fprintf(stderr,"ERROR! Time data has not been read.\n");
//        return;
//    }
//    std::map<NCAT,std::string>::iterator it_g;
//    std::map<NCAT,std::string>::iterator it_k;
//    string group, skill;
//    data *dt;
//    // open file
//    FILE *fid = fopen(file_name,"w");
//    fprintf(fid,"Group\tKC\ttime1\ttime2\ttimediff\ttime_lim_20HDWM\tOutcome\n");
//    // for all groups
//    for(NCAT g=0; g<param->nG; g++) {
//        // for all KCs
//        it_g = param->map_group_bwd->find(g);
//        for(NCAT k=0; k<param->g_numk[g]; k++) {
//            it_k = param->map_skill_bwd->find(k);
//            dt = param->g_k_data[g][k];
//            // for times from 2 to N
//            for(NDAT t=1; t<dt->n; t++) {
//                //                NPAR code = sec_to_linear_interval(dt->time[t]-dt->time[t-1], time_lim_20HDWM, sizeof(time_lim_20HDWM)/sizeof(int));
//                NPAR code = sec_to_9cat(dt->time[t-1], dt->time[t], time_lim_20HDWM, sizeof(time_lim_20HDWM)/sizeof(int));
//                fprintf(fid,"%s\t%s\t%d\t%d\t%d\t%d\t%d\n", it_g->second.c_str(), it_k->second.c_str(), dt->time[t-1], dt->time[t], (dt->time[t]-dt->time[t-1]), code, 1-param->dat_obs[ dt->ix[t] ]/*->get( dt->ix[t] )*/ );
//            }// for times from 2 to N
//        }// for all KCs
//    }// for all groups
//    // close file
//    fclose(fid);
//}
//

// penalties

// // uniform -- no longer up to date
//NUMBER L2penalty(param* param, NUMBER w) {
//    NUMBER penalty_offset = 0.5;
//    return (param->C != 0)? 0.5*param->C*fabs((w-penalty_offset)) : 0;
//}

// pre-specified
NUMBER L2penalty(NUMBER C, NUMBER w, NUMBER Ccenter) {
    return 0.5*C*fabs((w-Ccenter));
}

// for fitting larger portions first
int compareSortBitInv (const void * a, const void * b) {
    return -( ((sortbit*)a)->ndat - ((sortbit*)b)->ndat );
}


// random NUMBER in range
NUMBER NRand(NUMBER NMin, NUMBER NMax) {
    NUMBER N = (NUMBER)rand() / RAND_MAX;
    return NMin + N * (NMax - NMin);
}
