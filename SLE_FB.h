//
//  SLE_FB.h
//  Lattice enumeration via Forward-Backward programming
//
//  Created by moulay abdella chkifa on 9/23/22.
//  Copyright © 2022 moulay abdella chkifa. All rights reserved.
//

#ifndef SLE_FB_h
#define SLE_FB_h

#include "alias.h"
#include "chebyshev.h"

/*
 The FB_array is the main datastructure in the HMJ paper,
 It is an array consisting in (n+1) columns of size 2^n which
 are used in order to propagate Forward/Backward computation
 (that why FB) of the algorithm logic
 */

struct FB_array {
    
    static size_t n;
    std::vector<vec_real> Matrix;
    FB_array(size_t n_){
        n=n_;
        Matrix =  std::vector<vec_real>(n_+1) ;
        for (size_t i=0; i<n_+1; i++)
            Matrix[i] = vec_real(pow(2,n_));
    }
    
    FB_array(size_t n_,size_t m_){
        Matrix =  std::vector<vec_real>(n_) ;
        for (size_t i=0; i<n_; i++)
            Matrix[i] = vec_real(m_);
    }
    
    vec_real& operator[](size_t i) { return Matrix[i]; }
    const vec_real& operator[](size_t i) const { return Matrix[i]; }
    
    const vec_real& front() const { return Matrix[0]; }
    const vec_real& back() const { return Matrix[n]; }
};
size_t FB_array::n;

/*
 ****************************************************************************
 Updating functions acting on vector slices as described in HMJ paper
 UpdateForward reflect the application of functions tau_m to X
 UpdateBackward_1 reflect the application of function rho_m to B and C
 UpdateBackward_2 reflect the application of functions phi_m and psi_m
 to B, C and X
 ****************************************************************************
 */

/*
 Forevery integer i, vector_Dn[i] is a vector in dimension
 2^i. They are used in the recurrence definition of matrices
 A1,A2,etc
 */

inline void UpdateForward (size_t m, size_t p, FB_array& curr_X){
    size_t M     = pow(2,m-1);
    size_t pos1  = pow(2,m)*p;
    size_t pos2  = pow(2,m)*p + M ;
    
    
    // input iterators
    vec_real::iterator it_in_1=curr_X[m-1].begin()+pos1;
    vec_real::iterator it_in_2=curr_X[m-1].begin()+pos2;
    vec_real::iterator itD=vector_Dn[m-1].begin();
    
    // output iterators
    vec_real::iterator it1=curr_X[m].begin()+pos1;
    vec_real::iterator it2=curr_X[m].begin()+pos2;


    for (int i=0; i<M; i++) {
        *it1=(*it_in_1) + (*it_in_2) * (*itD);
        *it2=(*it_in_1) - (*it_in_2) * (*itD);
        it1++;it2++;it_in_1++;it_in_2++;itD++;
    }
    

    /*
     // if you use slicing and a science library
     // such as numpy in python, you would simply use
     // something like
    slice s1(pos1,M);
    slice s2(pos2,M);
    curr_X[m][s1] = curr_X[m-1][s1] + curr_X[m-1][s2] * vector_Dn[m-1];
    curr_X[m][s2] = curr_X[m-1][s1] - curr_X[m-1][s2] * vector_Dn[m-1];
    */
}

inline void UpdateBackward_1 (size_t m, size_t p, FB_array& curr_B, FB_array& curr_C){
    size_t M    = pow(2,m-1);
    size_t pos1 = pow(2,m)*p;
    size_t pos2 = pow(2,m)*p + M ;
    
    
    // input iterators
    vec_real::iterator itB1=curr_B[m].begin()+pos1;
    vec_real::iterator itB2=curr_B[m].begin()+pos2;
    
    vec_real::iterator itC1=curr_C[m].begin()+pos1;
    vec_real::iterator itC2=curr_C[m].begin()+pos2;
    
    // output iterators
    vec_real::iterator itB=curr_B[m-1].begin()+pos1;
    vec_real::iterator itC=curr_C[m-1].begin()+pos1;
    
    for (int i=0; i<M; i++) {
        *itB=0.5*((*itB1)+(*itB2));
        *itC=0.5*((*itC1)+(*itC2));
        itB++;itB1++;itB2++;
        itC++;itC1++;itC2++;
    }
    
    /*
     // if you use slicing and a science library
     // such as numpy in python, you would simply use
     // something like
    slice s1(pos1,M);
    slice s2(pos2,M);
    
    curr_B[m-1][s1] = 0.5*(curr_B[m][s1]+ curr_B[m][s2]);
    curr_C[m-1][s1] = 0.5*(curr_C[m][s1]+ curr_C[m][s2]);
     */
}

inline void UpdateBackward_2(size_t m, size_t p,
                             FB_array& curr_X,
                             FB_array& curr_B,
                             FB_array& curr_C){
    
    size_t M     = pow(2,m-1);
    size_t pos1  = pow(2,m)*p;
    size_t pos2  = pow(2,m)*p + M ;
    
    
    // input iterators
    vec_real::iterator itD=vector_Dn[m-1].begin();
    vec_real::iterator itX=curr_X[m-1].begin()+pos1;
    
    vec_real::iterator itB1=curr_B[m].begin()+pos1;
    vec_real::iterator itB2=curr_B[m].begin()+pos2;
    
    vec_real::iterator itC1=curr_C[m].begin()+pos1;
    vec_real::iterator itC2=curr_C[m].begin()+pos2;
    
    // output iterators
    vec_real::iterator itB=curr_B[m-1].begin()+pos2;
    vec_real::iterator itC=curr_C[m-1].begin()+pos2;

    for (int i=0; i<M; i++) {
        *itB=std::max((*itB1)-(*itX),(*itX)-(*itC2))/(*itD);
        *itC=std::min((*itX)-(*itB2),(*itC1)-(*itX))/(*itD);
        itB++;itB1++;itB2++;
        itC++;itC1++;itC2++;
        itX++;itD++;
    }
    
    /*
     // if you use slicing and a science library
     // such as numpy in python, you would simply use
     // something like
    slice s1(pos1,M);
    slice s2(pos2,M);
    curr_B[m-1][s2] = max_expr(curr_B[m][s1] - curr_X[m-1][s1],
                               curr_X[m-1][s1] - curr_C[m][s2])/vector_Dn[m-1];
    curr_C[m-1][s2] = min_expr(curr_X[m-1][s1] - curr_B[m][s2],
                               curr_C[m][s1] - curr_X[m-1][s1])/vector_Dn[m-1];
    */
}


/*
 *********************************************************************************
 Main function is a recursive function that encode a d-fold loop where at
 each loop, the FB_array X,B and C are modified according to the
 Forward-Backward logic described in the HMJ paper
 *********************************************************************************
 */

/*
 *********************************************************************************
 First identify p and r which are used as index of columns (p leading columns)
 and rows (2^r rows) and which are modified in the Forward Backward traversal
 of  X, B and C
 *********************************************************************************
 */

std::pair<size_t, size_t> decompose_rp (size_t i){
    if (i==0)
        return std::make_pair(0, i);
    size_t ri=0;
    while (i%2==0){
        ri+=1;
        i = i/2;
    }
    return  std::make_pair(ri, i);
}

/*
 ***************************************************************************
 Update functions for X, B and C
 **************************************************************************
 */
inline void UpdateX (size_t i, FB_array& curr_X){
    size_t r = decompose_rp(i).first;
    size_t curr_p;
    for (size_t m=1; m<=r; m++) {
        curr_p = i/pow(2, m) - 1;
        UpdateForward(m, curr_p, curr_X);
    }
}

inline void UpdateBC (size_t i,
                      FB_array& curr_X,
                      FB_array& curr_B,
                      FB_array& curr_C){
    if(i==0){return;}
    std::pair<size_t,size_t> rp = decompose_rp(i); // r and p=2q+1
    size_t r = rp.first, p= rp.second;
    size_t curr_m  = r+1; // r+1
    size_t curr_q  = (p-1)/2;// q
    UpdateBackward_2 (curr_m, curr_q, curr_X, curr_B, curr_C);
    
    for (size_t m=r; m>=1; m--) {
        curr_m = m;
        curr_q = i/pow(2,m);
        UpdateBackward_1(curr_m, curr_q, curr_B, curr_C);
    }
}

static size_t count_enumerated= 0 ;
static size_t count_dead_ends = 0;
static double integrale= 0. ;

/*
 ********************************************************************
 The main recursive function in the enumeration
 The action in the inner-most loop is either to
 store the node or to evaluate an integrand
 ********************************************************************
 */
#define _store_node
// LC used for latticeContainer, either a vector, a list, or set
template <typename LC>
void enumerate(size_t i,
               FB_array& curr_X,
               FB_array& curr_B,
               FB_array& curr_C,
               LC& curr_Pn){
    size_t d = curr_X[0].size();
    if (i == d ){
        count_enumerated+=1;
#ifdef _store_node
        //UpdateX(i,curr_X);
        //insert_node(curr_Pn, curr_X.front(),curr_X.back());
#else
        //insert_node(curr_Pn, vec_i(0),vec_d(0));
#endif
        // if (is_units (curr_X.front()))
        //print_array (curr_X.front());
        //count_succes+=1;
        //   integrale += pow(Nm(curr_X.back()),1./8)/200;
    }
    else{
        UpdateX (i,curr_X);
        UpdateBC(i, curr_X, curr_B, curr_C);
        
        int xi_min =  ceil(curr_B[0][i]);
        int xi_max = floor(curr_C[0][i]);
        for (int xi=xi_min; xi<=xi_max; xi++) {
            curr_X[0][i]=xi;
            enumerate(i+1, curr_X, curr_B, curr_C, curr_Pn);
        }
        if (xi_min == xi_max+1 ){count_dead_ends+=1;}
    }
}

/*
 *****************************************************************************
 Sequential Enumeration Algorithm using forward backward
 logic is called SLE_FB
 Enumeration algorithms : d-fold loop
 ****************************************************************************/
void resetAlgorithmParams(size_t n,
                          const vec_real& b,
                          const vec_real& c,
                          FB_array& curr_B,
                          FB_array& curr_C){
    
    std::copy(b.cbegin(),b.cend(),curr_B[n].begin());
    std::copy(c.cbegin(),c.cend(),curr_C[n].begin());
    for (size_t m=n; m>=1; m--)
        UpdateBackward_1(m, 0, curr_B, curr_C);
}

/*
 ***************************************************************
 LC used for latticeContainer, either a vector ,a list, or a set
 enumerate x s.t. b <= A_n x <= c or particular case
 *************************************************************
 */

template <typename LC>
LC SLE_FB(size_t n, const vec_real& b, const vec_real& c){
    FB_array curr_X(n), curr_B(n), curr_C(n);
    resetAlgorithmParams(n, b, c, curr_B, curr_C);
    
    LC curr_Pn;
    enumerate(0, curr_X, curr_B, curr_C, curr_Pn);
    return curr_Pn;
}

/*
 ***************************************************
 Sequential algrithm, other convenient calls  //
 ***************************************************
 */
template <typename LC>
LC SLE_FB(size_t n, const std::pair<vec_real,vec_real>& bc)
{return SLE_FB<LC>(n, bc.first, bc.second);}

template <typename LC>
LC SLE_FB(size_t n, const vec_real& v){
    vec_real w(v);
    for_each(w.begin(), w.end(), [](double x) {x*=-1;});
    return SLE_FB<LC>(n, w, v);}

template <typename LC>
LC SLE_FB(size_t n, double alpha)
{return SLE_FB<LC>(n, vec_real(pow(2,n),-alpha),vec_real(pow(2,n),alpha));}


#endif /* SLE_BF_h */
