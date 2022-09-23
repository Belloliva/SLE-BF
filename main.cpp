//
//  main.cpp
//  Lattice enumeration via Forward-Backward programming
//
//  Created by moulay abdella chkifa on 9/22/22.
//  Copyright © 2022 moulay abdella chkifa. All rights reserved.
//

#include <cmath>
#include <iostream>
#include <vector>
#include <numeric>
#include "SLE_FB.h"
#include "timer.h"

using namespace std;

int enumerate_hadamard(size_t n,size_t scaling_factor){
    size_t N=scaling_factor;
    size_t d= pow(2,n);
    double rho = 0.5*pow(N,1./d)*sqrt(d);
    std::vector<int> Pn = SLE_FB<std::vector<int>>(n,rho);
    // function for for sequential lattice enumeration using
    // forward backward computation as explained in my HMJ paper
    return 0;
}

int enumerate_tchebechev(size_t n,size_t scaling_factor){
    size_t N=scaling_factor;
    size_t d= pow(2,n);
    double rho = pow(N/sqrt(2),1./d)*sqrt(0.5*d);
    std::vector<int> Pn = SLE_FB<std::vector<int>>(n,rho);
    return 0;
}

int main(int argc, const char * argv[]) {
    int choice=1;
    if (choice==1)
        vector_Dn =initialize_hadamard(20);
    if (choice==0)
        vector_Dn =initialize_tchebechev(20);
    
    // vector_Dn is a global variable, see chebyshev.h
    // use n=1 if you want to enumerate Hadamard lattice
    // use n=0 if you want to enumerate Chebechev lattice
    
    size_t n=4;// we enumerate in dimension d=2^n=16
    std::vector<size_t> vec_m(26);
    std::iota(vec_m.begin(), vec_m.end(), 1);
    
    std::vector<float> runtimes;
    std::vector<size_t> num_enumerated_Points;
    
    for (size_t m:vec_m) {
        double N = pow(2,m);
        std::cout <<"   -----" <<" m=" << m <<", N=" << N <<"   -----" << std::endl;
        timer t ("for current enumeration");
        
        // reset counts,these are static signed integer, see in SLE_FB
        count_enumerated= 0;
        count_dead_ends = 0;
        
        // main enumeration
        if (choice==1)
            enumerate_hadamard(n,N);
        if (choice==0)
            enumerate_tchebechev(n, N);
        
        runtimes.push_back(t.timer_value());
        num_enumerated_Points.push_back(count_enumerated);
    
        std::cout << "Number of enumerated lattice points: " << count_enumerated << std::endl;
    }
    return 0;
}
