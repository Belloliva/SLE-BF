//
//  chebyshev.h
//  Lattice enumeration via Forward-Backward programming
//
//  Created by moulay abdella chkifa on 9/22/22.
//  Copyright © 2022 moulay abdella chkifa. All rights reserved.
//

#ifndef chebyshev_h
#define chebyshev_h

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#include <math.h>
#include "alias.h"

/*
 *****************************************************************
 Permuted Tchebechev abscissas
 First we define the permutations of rows and columns as
 described in our HMJ paper
 *****************************************************************
 */

size_t rowPermutation(const size_t k){
    if (k==0){return 0;}
    else{
        size_t n = floor( log2(k));
        size_t N = pow(2,n);
        return 2*N -1-rowPermutation(k - N);
    }
}

size_t colPermutation(const size_t n, const size_t k){
    if (n==0 && k==0){return 0;}
    else{
        size_t N = pow(2,n-1);
        if (k < N)
            return 2*colPermutation(n-1, k);
        else
            return 2*colPermutation(n-1, k-N)+1;
    }
}

/*
 ******************************************************************************
 Permuted tchebechev abscissas, these are diagonal matrices Dn
 we use on the HMJ paper
 ******************************************************************************
 */

double Tcheb_root (size_t n, size_t i){
    return 2*cos((i+0.5)*M_PI/pow(2,n));
}

vec_real Tcheb_roots ( size_t n){
    size_t N = pow(2,n);
    vec_real wn (N);
    for (size_t i = 0; i <N; i++)
        wn[i] = Tcheb_root (n+1,rowPermutation(i));
    return wn;
}

vec_real Tcheb_roots_inv (size_t n){
    vec_real wn = Tcheb_roots (n);
    for (size_t i = 0; i < wn.size(); i++)
        wn[i] = 1./wn[i];
    return wn;
}

/*
 ****************************************
 Global variables
 ****************************************
 */
static std::vector<vec_real> vector_Dn;
static std::vector<vec_real> vector_inverse_Dn;
/*
 ******************************************
 and how they get initialized in main.cpp
 ******************************************
 */
std::vector<vec_real> initialize_tchebechev(size_t nmax){
    std::vector<vec_real> Tcheb_roots_table(nmax);
    for (size_t i = 0; i < nmax; i++)
        Tcheb_roots_table[i] = Tcheb_roots (i);
    return Tcheb_roots_table;
}


std::vector<vec_real> initialize_hadamard(size_t nmax){
    std::vector<vec_real> units_table(nmax);
    for (size_t i = 0; i < nmax; i++){
        vec_real wn ((size_t)pow(2,i),1);
        units_table[i] =wn;
    }
    return units_table;
}

#endif /* chebyshev_h */
