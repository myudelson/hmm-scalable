//
//  SparseArray2D.h
//  HMM
//
//  Created by Yudelson, Michael on 9/23/13.
//
//

#ifndef __HMM__SparseArray2D__
#define __HMM__SparseArray2D__

#include <stdlib.h>
#include "utils.h"

template <typename T>
class SparseArray2D {
public:
	SparseArray2D();
	~SparseArray2D();
    static NDAT binsearch(register const T *key, const void *base0, NDAT nmemb, NDAT *lim);
    void set(NDAT r, NDAT c, T value);
    
private:
    T *elements; // all elements
    NDAT size; // all elements
    NDAT n_rows;
    NDAT *row_p; // pointer to rows in index terms
    NDAT *row_ix; // row coordinate index
    NDAT *col_ix; // colum coordinate indexes
    NDAT row_size(NDAT r); // for row coordinate
    NDAT row_size_ix(NDAT r_ix); // for row index
    void add_row(NDAT r, NDAT c, T v); // add row and column and value in one pack
};

template <typename T>
SparseArray2D<T>::SparseArray2D() {
    this->elements = NULL;
    this->size     = 0;
    this->n_rows   = 0;
    this->row_p    = NULL;
    this->row_ix   = NULL;
    this->col_ix   = NULL;
}

template <typename T>
SparseArray2D<T>::~SparseArray2D() {
    free(this->elements);
    free(this->row_p);
    free(this->row_ix);
    free(this->col_ix);
}

template<typename T>
NDAT SparseArray2D<T>::binsearch(
               register const T *key,
               const void *base0,
               NDAT nmemb, NDAT *lim) {
    
	register const char *base = (const char*)base0;
	register const char *base00 = base;
//	register NDAT lim;
    register NDAT size = sizeof(T);
	register int cmp;
	register const void *p;
    
	for (*lim = nmemb; *lim != 0; *lim >>= 1) {
		p = base + (*lim >> 1) * size;
		cmp = *key -(*((T*)p));// (*compar)(key, p);
		if (cmp == 0) {
			return ((char*)p-base00)/size;
        }
		if (cmp > 0) {	/* key > p: move right */
			base = (char *)p + size;
			(*lim)++;
		}
        /* else move left */
	}
	return -1;
}

template <typename T>
NDAT SparseArray2D<T>::row_size(NDAT r) {
    NDAT lim;
    NDAT ix =  binsearch(&r, this->row_ix, this->n_row, sizeof(NDAT), &lim);
    if( ((int)ix)==-1 )
        return 0;
    else {
        if( this->nrows==1)
            return this->size;
        else if( ix==(this->n_rows-1) )
            return this->size - this->row_ix[ix-1];
        else
            return this->row_ix[ix] - this->row_ix[ix-1];
    }
}

template <typename T>
NDAT SparseArray2D<T>::row_size_ix(NDAT r_ix) {
    if( ((int)r_ix)==-1 )
        return 0;
    else {
        if( this->nrows==1)
            return this->size;
        else if( r_ix==(this->n_rows-1) )
            return this->size - this->row_ix[r_ix-1];
        else
            return this->row_ix[r_ix] - this->row_ix[r_ix-1];
    }
}


template <typename T>
void SparseArray2D<T>::set(NDAT r, NDAT c, T value) {
    // row exist?
    NDAT r_ix = binsearch(&r, this->row_ix, this->n_row, sizeof(NDAT));
    // add row
    if( (int)r_ix == -1) {
        
    }
    
    // column exist?
    // add column
    // set value
}

// add row and column and value in one pack
template <typename T>
void SparseArray2D<T>::add_row(NDAT r, NDAT c, T v) {
    
}

#endif /* defined(__HMM__SparseArray2D__) */
