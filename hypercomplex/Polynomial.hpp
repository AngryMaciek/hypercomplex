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

long int RingInverse(const long int &x, const long int &mod) {
    long int y = x % mod;
    if (y < 0) y += mod;
    for (unsigned int i=1; i < mod; i++) {
        if ((y*i) % mod == 1) return i;
    }
    throw std::invalid_argument("non-invertible element");
}

template <const unsigned int MaxDeg>
class Polynomial {
 private:
    long int *coefficients;

 public:
    // Polynomial main constructor
    explicit Polynomial(const long int* arr) {
        coefficients = new long int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = arr[i];
    }

    // Polynomial copy constructor
    Polynomial(const Polynomial &P) {
        coefficients = new long int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = P[i];
    }

    // Polynomial default constructor
    Polynomial() {
        coefficients = new long int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = 0;
    }

    // Polynomial destructor
    ~Polynomial() {
        delete[] coefficients;
    }

    // overloaded = operator
    Polynomial& operator= (const Polynomial &P) {
        // self-assignment guard
        if (this == &P) return *this;
        // reassign
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = P[i];
        // return the existing object so we can chain this operator
        return *this;
    }

    // overloaded - unary operator
    Polynomial operator-() const {
        long int* temparr = new long int[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = -coefficients[i];
        Polynomial<MaxDeg> P(temparr);
        delete[] temparr;
        return P;
    }

    // overloaded [] operator
    long int& operator[](const unsigned int i) const {
        assert(0 <= i && i <= MaxDeg);
        return coefficients[i];
    }
};

// overloaded == operator
template <const unsigned int MaxDeg>
bool operator==(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    for (unsigned int i=0; i <= MaxDeg; i++) {
        if (P1[i] != P2[i]) return false;
    }
    return true;
}

// overloaded != operator
template <const unsigned int MaxDeg>
bool operator!=(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    return !(P1 == P2);
}

// overloaded << operator
template <const unsigned int MaxDeg>
std::ostream& operator<< (std::ostream &os, const Polynomial<MaxDeg> &P) {
    for (unsigned int i=0; i < MaxDeg; i++) os << P[i] << ",";
    os << P[MaxDeg];
    return os;
}

// overloaded * operator: multiplication by a scalar
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator*(const long int &x, const Polynomial<MaxDeg> &P) {
    long int *temparr = new long int[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P[i] * x;
    Polynomial<MaxDeg> p(temparr);
    delete[] temparr;
    return p;
}

// overloaded + binary operator
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator+(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    long int *temparr = new long int[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P1[i] + P2[i];
    Polynomial<MaxDeg> p(temparr);
    delete[] temparr;
    return p;
}

// overloaded - binary operator
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator-(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    long int *temparr = new long int[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P1[i] - P2[i];
    Polynomial<MaxDeg> p(temparr);
    delete[] temparr;
    return p;
}

// overloaded * binary operator:
// convolution multiplication in a polynomial quotient ring Z/(x^N-1)
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator*(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    long int *prod = new long int[2*MaxDeg+1]();
    long int *conv = new long int[MaxDeg+1]();
    for (unsigned int i=0; i <= MaxDeg; i++) {
        for (unsigned int j=0; j <= MaxDeg; j++)
            prod[i+j] += P1[i]*P2[j];
    }
    for (unsigned int i=0; i < 2*MaxDeg+1; i++) conv[i%(MaxDeg+1)] += prod[i];
    Polynomial<MaxDeg> p(conv);
    delete[] prod;
    delete[] conv;
    return p;
}

// overloaded % operator: reduce coefficients modulo a scalar
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator%(const Polynomial<MaxDeg> &P, const long int &x) {
    long int *temparr = new long int[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) {
        temparr[i] = P[i] % x;
        if (temparr[i] < 0) temparr[i] += x;
    }
    Polynomial<MaxDeg> p(temparr);
    delete[] temparr;
    return p;
}

// centered lift of a polynomial in the quotient ring Z/nZ/ / (x^N-1)
template <const unsigned int MaxDeg>
void CenteredLift(const Polynomial<MaxDeg> &P, const long int &mod) {
    long int lower = -mod/2;
    long int upper = mod/2;
    for (unsigned int i = 0; i <= MaxDeg; i++) {
        if (mod % 2) {  // odd: <lower, upper>
            if (P[i] < lower) P[i] = P[i] + mod;
            if (P[i] > upper) P[i] = P[i] - mod;
        } else {  // even: (lower, upper>
            if (P[i] <= lower) P[i] = P[i] + mod;
            if (P[i] > upper) P[i] = P[i] - mod;
        }
    }
}

#endif  // HYPERCOMPLEX_POLYNOMIAL_HPP_
