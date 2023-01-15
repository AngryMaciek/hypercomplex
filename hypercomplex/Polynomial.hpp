// Copyright 2022 <Maciej Bak>
/*! \file */
/*
###############################################################################
#
#   Polynomial helper class and functions for the cryptosystems.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Department_of_Mathematics_City_University_of_London
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
#include <iostream>

int64_t RingInverse(const int64_t x, const int64_t mod) {
    int64_t y = x % mod;
    if (y < 0) y += mod;
    for (unsigned int i=1; i < mod; i++) {
        if ((y*i) % mod == 1) return i;
    }
    throw std::invalid_argument("non-invertible element");
}

template <const unsigned int MaxDeg>
class Polynomial {
 private:
    int64_t coefficients[MaxDeg+1];

 public:
    // Polynomial main constructor
    explicit Polynomial(const int64_t* arr) {
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = arr[i];
    }

    // Polynomial copy constructor
    Polynomial(const Polynomial &P) {
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = P[i];
    }

    // Polynomial default constructor
    Polynomial() {
        for (unsigned int i=0; i <= MaxDeg; i++) coefficients[i] = 0;
    }

    // Polynomial destructor
    ~Polynomial() {
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
        int64_t temparr[MaxDeg+1];
        for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = -coefficients[i];
        Polynomial<MaxDeg> P(temparr);
        return P;
    }

    // overloaded [] operator (const)
    int64_t const & operator[](const unsigned int i) const {
        assert(0 <= i && i <= MaxDeg);
        return coefficients[i];
    }

    // overloaded [] operator (non-const)
    int64_t & operator[](const unsigned int i) {
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
Polynomial<MaxDeg> operator*(const int64_t x, const Polynomial<MaxDeg> &P) {
    int64_t temparr[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P[i] * x;
    Polynomial<MaxDeg> p(temparr);
    return p;
}

// overloaded + binary operator
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator+(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    int64_t temparr[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P1[i] + P2[i];
    Polynomial<MaxDeg> p(temparr);
    return p;
}

// overloaded - binary operator
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator-(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    int64_t temparr[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) temparr[i] = P1[i] - P2[i];
    Polynomial<MaxDeg> p(temparr);
    return p;
}

// overloaded * binary operator:
// convolution multiplication in a polynomial quotient ring Z/(x^N-1)
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator*(
    const Polynomial<MaxDeg> &P1,
    const Polynomial<MaxDeg> &P2
) {
    int64_t prod[2*MaxDeg+1];
    int64_t conv[MaxDeg+1];
    for (unsigned int i=0; i < 2*MaxDeg+1; i++) prod[i] = 0;
    for (unsigned int i=0; i <= MaxDeg; i++) conv[i] = 0;
    for (unsigned int i=0; i <= MaxDeg; i++) {
        for (unsigned int j=0; j <= MaxDeg; j++)
            prod[i+j] += P1[i]*P2[j];
    }
    for (unsigned int i=0; i < 2*MaxDeg+1; i++) conv[i%(MaxDeg+1)] += prod[i];
    Polynomial<MaxDeg> p(conv);
    return p;
}

// overloaded % operator: reduce coefficients modulo a scalar
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> operator%(const Polynomial<MaxDeg> &P, const int64_t x) {
    int64_t temparr[MaxDeg+1];
    for (unsigned int i=0; i <= MaxDeg; i++) {
        temparr[i] = P[i] % x;
        if (temparr[i] < 0) temparr[i] += x;
    }
    Polynomial<MaxDeg> p(temparr);
    return p;
}

// centered lift of a polynomial in a modular quotient ring Z/nZ/ / (x^N-1)
template <const unsigned int MaxDeg>
void CenteredLift(Polynomial<MaxDeg> *P, const int64_t mod) {
    int64_t lower = -mod/2;
    int64_t upper = mod/2;
    for (unsigned int i = 0; i <= MaxDeg; i++) {
        if (mod % 2) {  // odd: <lower, upper>
            if ((*P)[i] < lower) (*P)[i] = (*P)[i] + mod;
            if ((*P)[i] > upper) (*P)[i] = (*P)[i] - mod;
        } else {  // even: (lower, upper>
            if ((*P)[i] <= lower) (*P)[i] = (*P)[i] + mod;
            if ((*P)[i] > upper) (*P)[i] = (*P)[i] - mod;
        }
    }
}

// polynomial inverse in a modular quotient ring Z/nZ/ / (x^N-1)
// adapted from: https://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
template <const unsigned int MaxDeg>
Polynomial<MaxDeg> RingInverse(
    const Polynomial<MaxDeg> &P,
    const int64_t &mod
) {
    // make sure to shift coefficients to: [0,mod-1]
    Polynomial<MaxDeg> P_ = P % mod;
    //
    // define placeholders & init values
    //
    Polynomial<MaxDeg+1> q;
    Polynomial<MaxDeg+1> temp;
    //
    Polynomial<MaxDeg+1> t;
    Polynomial<MaxDeg+1> newt;
    newt[0] = 1;
    //
    Polynomial<MaxDeg+1> r;
    r[0] = mod-1;
    r[MaxDeg+1] = 1;
    unsigned int deg_r = MaxDeg+1;
    for (unsigned int i=0; i <= MaxDeg; i++) temp[i] = P_[i];
    Polynomial<MaxDeg+1> newr(temp);
    unsigned int deg_newr = 0;
    for (unsigned int i=0; i <= MaxDeg+1; i++) {
        if (newr[i]) deg_newr = i;
    }

    // algorithm loop
    while (deg_newr > 0) {
        // division loop
        while (deg_r >= deg_newr) {
            for (unsigned int i=0; i <= MaxDeg+1; i++) temp[i] = 0;
            for (unsigned int i=0; i <= deg_newr; i++)
                temp[i + deg_r - deg_newr] = newr[i];
            q[deg_r - deg_newr] =
                (r[deg_r] * RingInverse(temp[deg_r], mod)) % mod;
            for (unsigned int i=0; i <= deg_r; i++)
                temp[i] = (temp[i] * q[deg_r - deg_newr]) % mod;
            temp = temp % mod;
            for (unsigned int i=0; i <= deg_r; i++) r[i] = r[i] - temp[i];
            r = r % mod;
            deg_r -= 1;
        }

        // re-assign labels after one round of the algorithm:
        temp = newr;
        newr = r;
        r = temp;
        //
        temp = newt;
        newt = (t - (q * newt)) % mod;
        t = temp;

        // reset variables:
        for (unsigned int i=0; i <= MaxDeg+1; i++) q[i] = 0;
        for (unsigned int i=0; i <= MaxDeg+1; i++) {
            if (newr[i]) deg_newr = i;
        }
        for (unsigned int i=0; i <= MaxDeg+1; i++) {
            if (r[i]) deg_r = i;
        }
    }

    // after the algorithm: deg of the last remainder < deg(P)
    assert(newr[MaxDeg+1] == 0);

    // construct the inverse polynomial
    int64_t multiplier = RingInverse(newr[0], mod);
    Polynomial<MaxDeg> inverse;
    for (unsigned int i=0; i <= MaxDeg; i++) inverse[i] = newt[i];
    inverse = (multiplier * inverse) % mod;

    // test: P * 1/P == 1
    Polynomial<MaxDeg> unity;
    unity[0] = 1;
    assert((inverse * P) % mod == unity);

    return inverse;
}

#endif  // HYPERCOMPLEX_POLYNOMIAL_HPP_
