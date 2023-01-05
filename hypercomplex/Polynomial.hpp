// Copyright 2022 <Maciej Bak>
/*! \file */
/*
###############################################################################
#
#   Polynomial helper class and functions for the cryptosystems.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 02-12-2022
#   LICENSE: Apache 2.0
#
###############################################################################
*/

// mark that we included this library
#ifndef HYPERCOMPLEX_POLYNOMIAL_HPP_
#define HYPERCOMPLEX_POLYNOMIAL_HPP_

#include <cassert>
#include <cstdlib>
#include <iostream>

template <const unsigned int MaxDeg>
class Polynomial {
 private:
    int *coefficients;

 public:
    // Polynomial main constructor
    explicit Polynomial(const int* arr) {
        coefficients = new int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = arr[i];
    }

    // Polynomial copy constructor
    Polynomial(const Polynomial &P) {
        coefficients = new int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = P[i];
    }

    // Polynomial default constructor
    Polynomial() {
        coefficients = new int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = 0;
    }

    // Polynomial destructor
    ~Polynomial() {
        delete[] coefficients;
    }
};

#endif  // HYPERCOMPLEX_POLYNOMIAL_HPP_
