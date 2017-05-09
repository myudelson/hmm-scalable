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

/*
 * helper auto-expanded array without re-allocating the memory, by allocating 
 * chunks - stripes
 */

//#include "utils.h"
#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <math.h>


#ifndef STRIPEDARRAY_H
#define STRIPEDARRAY_H

#define NDAT_MAX INT_MAX
typedef signed int NDAT;  // number of data rows, now 4 bill max

template <typename T>
class StripedArray {
public:
	StripedArray();
	StripedArray(NDAT _size/*, bool a_complex*/); // predefined size
	StripedArray(FILE *f, NDAT N);
//	StripedArray(bool a_complex);
	~StripedArray();
	NDAT getSize();
	void add(T value);
	T& operator [] (NDAT idx);
	T get(NDAT idx);
	void set(NDAT idx, T value);
	void clear();
    NDAT toBinFile(FILE* f);
    static NDAT arrayToBinFile(T* ar, NDAT size, FILE* f);
    T* toArray();
private:
	NDAT size; // linear
	NDAT stripe_size;
	NDAT nstripes;
	NDAT size_last_stripe;
//    bool complex; // element stored is a pointer to an array (should be deleted further)
	T** stripes;
	void addStripe();
};


template <typename T>
StripedArray<T>::StripedArray() {
	this->size = 0;
	this->stripe_size = 20000;
	this->nstripes = 0;
	this->size_last_stripe = 0;
	this->stripes = NULL;
//    this->complex = false;
}

// predefined size
template <typename T>
StripedArray<T>::StripedArray(NDAT _size/*, bool a_complex*/) {
	stripe_size = 20000;
//    complex = a_complex;
	size = _size;
	nstripes = (NDAT)ceil((double)size/stripe_size);
	size_last_stripe = (NDAT)fmod(size, (size_t)stripe_size);
	stripes = (T **) calloc((size_t)nstripes, sizeof(T*));
    NDAT num;
    for(NDAT i=0; i<nstripes; i++) {
        num = (i<(nstripes-1)?stripe_size:size_last_stripe);
        stripes[i] = (T*)calloc((size_t)num, sizeof(T)); // alloc data
    }
}

// from file
template <typename T>
StripedArray<T>::StripedArray(FILE *f, NDAT N) {
	stripe_size = 20000;
//    complex = false;
	size = N;
	nstripes = (NDAT)ceil((double)N/stripe_size);
	size_last_stripe = (NDAT)fmod(N, stripe_size);
	stripes = (T **) calloc((size_t)nstripes, sizeof(T*));
    NDAT num;
    NDAT nread;
    for(NDAT i=0; i<nstripes; i++) {
        num = (i<(nstripes-1)?stripe_size:size_last_stripe);
        stripes[i] = (T*)calloc((size_t)num, sizeof(T)); // alloc data
        nread = (NDAT)fread (stripes[i], sizeof(T), (size_t)num, f);
        if(nread != num) {
            fprintf(stderr,"Error reading data from file\n");
            return;
        }
    }
}

//template <typename T>
//StripedArray<T>::StripedArray(bool a_complex) {
//	size = 0;
//	stripe_size = 20000;
//	nstripes = 0;
//	size_last_stripe = 0;
//	stripes = NULL;
//    complex = a_complex;
//}

template <typename T>
StripedArray<T>::~StripedArray() {
	for(NDAT i=0; i<this->nstripes;i++) {
//        if(this->complex)
//            for(NDAT j=0; j<( (i<(this->nstripes-1))?this->stripe_size:this->size_last_stripe );j++)
//                free(  (void *)(this->stripes[i][j]) );
		free(this->stripes[i]);
    }
	free(this->stripes);
//	for(NDAT i=0; i<this->nstripes;i++)
//		delete [] this->stripes[i];
//	delete [] this->stripes;
}

template <typename T>
void StripedArray<T>::add(T value) {
	if(this->size == (NDAT_MAX-1)) {
		fprintf(stderr, "Error! Maximum array size reached.\n");
		return;
	}
	
	// add stripe if necessary
	if( this->nstripes==0 || this->size_last_stripe==this->stripe_size)
		this->addStripe();
	// place data
	this->stripes[this->nstripes-1][this->size_last_stripe++] = value;
	this->size++;
}

template <typename T>
void StripedArray<T>::clear() {
	for(NDAT i=0; i<this->nstripes;i++)
		free(this->stripes[i]);
	free(this->stripes);
	this->size = 0;
	this->stripe_size = 20000;
	this->nstripes = 0;
	this->size_last_stripe = 0;
	this->stripes = NULL;
}

template <typename T>
T& StripedArray<T>::operator [](NDAT idx) {
	if(idx>(this->size-1)) {
		fprintf(stderr, "   %d exceeds array size %d.\n",idx,this->size);
		return NULL;
	}
	NDAT idx_stripe = idx / this->stripe_size;
	NDAT idx_within = idx % this->stripe_size;
//	if(idx_within>(this->size_last_stripe-1)) {
//		fprintf(stderr, "Exception! Element index in last stripe %d exceeds sltipe size %d.\n",idx_within,this->size_last_stripe);
//		return NULL;
//	}
	return this->stripes[idx_stripe][idx_within];
}

template <typename T>
T StripedArray<T>::get(NDAT idx) {
	if(idx>(this->size-1)) {
		fprintf(stderr, "Exception! Element index %u exceeds array size %u\n",idx,this->size);
		return (T)0;
	}
	NDAT idx_stripe = idx / this->stripe_size;
	NDAT idx_within = idx % this->stripe_size;
//	if(idx_within>(this->size_last_stripe-1)) {
//		fprintf(stderr, "Exception! Element index in last stripe %d exceeds sltipe size %d\n",idx_within,this->size_last_stripe);
//		return NULL;
//	}
	return this->stripes[idx_stripe][idx_within];
}

template <typename T>
void StripedArray<T>::set(NDAT idx, T value) {
	if(idx>(this->size-1)) {
		fprintf(stderr, "Exception! Element index %u exceeds array size %u\n",idx,this->size);
		return;
	}
	NDAT idx_stripe = idx / this->stripe_size;
	NDAT idx_within = idx % this->stripe_size;
    //	if(idx_within>(this->size_last_stripe-1)) {
    //		fprintf(stderr, "Exception! Element index in last stripe %d exceeds sltipe size %d\n",idx_within,this->size_last_stripe);
    //		return NULL;
    //	}
	this->stripes[idx_stripe][idx_within] = value;
}


template <typename T>
NDAT StripedArray<T>::getSize() {
	return this->size;
}

template <typename T>
void StripedArray<T>::addStripe() {
	if(this->size_last_stripe<this->stripe_size && this->nstripes>0) {
		fprintf(stderr, "Exception! You can only add stripes when none exits or previous is full.\n");
		return;
	}
	
	this->stripes = (T **) realloc(this->stripes,(size_t)(++this->nstripes)*sizeof(T*)); // add pointer
//	T** new_stripes = Calloc(T*, (this->nstripes));
//	memcpy(new_stripes, this->stripes, sizeof(T*)*this->nstripes );
//	free(this->stripes);
//	this->stripes = Malloc(T*,this->nstripes+1);
//	memcpy(this->stripes, new_stripes, sizeof(T*)*this->nstripes );
//	this->nstripes++;
//	free(new_stripes);
	
//	this->stripes[this->nstripes-1] = Calloc(T, (size_t)this->stripe_size); // alloc data
	this->stripes[this->nstripes-1] = (T*)calloc((size_t)this->stripe_size, sizeof(T)); // alloc data
	this->size_last_stripe = 0; // reset counter
}

template <typename T>
NDAT StripedArray<T>::toBinFile(FILE *f) {
    NDAT num;
    NDAT nwrit;
    NDAT all_nwrit = 0;
    for(NDAT i=0; i<nstripes; i++) {
        num = (NDAT)(i<(nstripes-1))?stripe_size:size_last_stripe;
        nwrit = (NDAT)fwrite (stripes[i] , sizeof(T), (size_t)num, f);
        all_nwrit += nwrit;
        if(num != nwrit) {
            fprintf(stderr, "Error writing. Attempted to write %u but %u were written.\n",num, nwrit);
            return 0;
        }
    }
    return all_nwrit;
}

template <typename T>
NDAT StripedArray<T>::arrayToBinFile(T* ar, NDAT size, FILE* f) {
    return (NDAT)fwrite (ar, sizeof(T), (size_t)size, f);
}

template <typename T>
T* StripedArray<T>::toArray() {
    T *result = (T*)malloc( (size_t)this->size*sizeof(T));
    for(NDAT i=0; i<this->nstripes; i++) {
        T* dest = &result[ i*this->stripe_size ];
        size_t sz = sizeof(T)*( (i<(this->nstripes-1))?this->stripe_size:this->size_last_stripe );
        memcpy( dest, this->stripes[i], sz );
    }
    return result;
}

#endif
