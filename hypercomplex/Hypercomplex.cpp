// Copyright 2020 <Maciej Bak>
/*
###############################################################################
#
#   Hypercomplex library: implementation file
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include "Hypercomplex.h" // NOLINT

// Hypercomplex main constructor
Hypercomplex::Hypercomplex(const unsigned int arg_d, const float* arg_arr) {
    if (arg_d == 0) throw std::invalid_argument("invalid dimension");
    if ((arg_d & (arg_d - 1)) != 0) {
        throw std::invalid_argument("invalid dimension");
    }
    d = arg_d;
    arr = new float[arg_d];
    for (unsigned int i=0; i < arg_d; i++) arr[i] = arg_arr[i];
}

// Hypercomplex copy constructor
Hypercomplex::Hypercomplex(const Hypercomplex& H) {
    d = H.d;
    arr = new float[H.d];
    for (unsigned int i=0; i < H.d; i++) arr[i] = H.arr[i];
}

// Hypercomplex destructor
Hypercomplex::~Hypercomplex() {
    delete[] arr;
}

// calculate norm of the number
float Hypercomplex::norm() const {
    float result = 0.0;
    for (unsigned int i=0; i < d; i++) result += arr[i] * arr[i];
    return sqrt(result);
}

// calculate inverse of the number
Hypercomplex Hypercomplex::inv() const {
    if ((*this).norm() == 0) {
        throw std::invalid_argument("division by zero");
    } else {
        float norm2 = pow((*this).norm(), 2);
        float *temparr = new float[d];
        temparr[0] = arr[0] / norm2;
        for (unsigned int i=1; i < d; i++) temparr[i] = -arr[i] / norm2;
        Hypercomplex H = Hypercomplex(d, temparr);
        delete[] temparr;
        return H;
    }
}

// expand object to a higher dimension
Hypercomplex Hypercomplex::expand(const unsigned int arg_d) const {
    if (arg_d <= d) throw std::invalid_argument("invalid dimension");
    float *temparr = new float[arg_d]();  // zero-init
    for (unsigned int i=0; i < d; i++) temparr[i] = arr[i];
    Hypercomplex H = Hypercomplex(arg_d, temparr);
    delete[] temparr;
    return H;
}

// overloaded ~ operator
Hypercomplex Hypercomplex::operator~() const {
    Hypercomplex H = -(*this);
    H[0] = arr[0];
    return H;
}

// overloaded - unary operator
Hypercomplex Hypercomplex::operator-() const {
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) temparr[i] = -arr[i];
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded = operator
Hypercomplex& Hypercomplex::operator=(const Hypercomplex &H) {
    // self-assignment guard
    if (this == &H) return *this;
    // reassign
    d = H.d;
    for (unsigned int i=0; i < d; i++) arr[i] = H.arr[i];
    // return the existing object so we can chain this operator
    return *this;
}

// overloaded [] operator
float& Hypercomplex::operator[](const unsigned int i) const {
    assert(0 <= i && i < d);
    return arr[i];
}

// overloaded == operator
bool operator==(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) return false;
    for (unsigned int i=0; i < H1._(); i++) {
        if (H1[i] != H2[i]) return false;
    }
    return true;
}

// overloaded != operator
bool operator!=(const Hypercomplex &H1, const Hypercomplex &H2) {
    return !(H1 == H2);
}

// overloaded + binary operator
Hypercomplex operator+(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) throw std::invalid_argument("operand mismatch");
    unsigned int d = H1._();
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) temparr[i] = H1[i] + H2[i];
    Hypercomplex H = Hypercomplex(H1._(), temparr);
    delete[] temparr;
    return H;
}

// overloaded - binary operator
Hypercomplex operator-(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) throw std::invalid_argument("operand mismatch");
    unsigned int d = H1._();
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) temparr[i] = H1[i] - H2[i];
    Hypercomplex H = Hypercomplex(H1._(), temparr);
    delete[] temparr;
    return H;
}

// overloaded * binary operator
Hypercomplex operator*(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) throw std::invalid_argument("operand mismatch");
    unsigned int d = H1._();
    // recursion base:
    if (d == 1) {
        float temparr[] = { H1[0] * H2[0] };
        return Hypercomplex(1, temparr);
    }
    // shared objects:
    unsigned int halfd = d / 2;
    float *temparr = new float[d];
    // construct helper objects:
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i];
    Hypercomplex H1a = Hypercomplex(halfd, temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i+halfd];
    Hypercomplex H1b = Hypercomplex(halfd, temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i];
    Hypercomplex H2a = Hypercomplex(halfd, temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i+halfd];
    Hypercomplex H2b = Hypercomplex(halfd, temparr);
    // multiply recursively:
    Hypercomplex H1a2a = H1a * H2a;
    Hypercomplex H2b_1b = ~H2b * H1b;
    Hypercomplex H2b1a = H2b * H1a;
    Hypercomplex H1b2a_ = H1b * ~H2a;
    // construct the final object
    Hypercomplex Ha = H1a2a - H2b_1b;
    Hypercomplex Hb = H2b1a + H1b2a_;
    for (unsigned int i=0; i < halfd; i++) temparr[i] = Ha[i];
    for (unsigned int i=0; i < halfd; i++) temparr[i+halfd] = Hb[i];
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded ^ binary operator
Hypercomplex operator^(const Hypercomplex &H, const unsigned int x) {
    if (!(x)) {
        throw std::invalid_argument("zero is not a valid argument");
    } else {
        Hypercomplex Hx = Hypercomplex(H);
        for (unsigned int i=0; i < x-1; i++) Hx = Hx * H;
        return Hx;
    }
}

// overloaded / binary operator
Hypercomplex operator/(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) throw std::invalid_argument("operand mismatch");
    // division H1 / H2 is implemented as H1 * 1/H2
    Hypercomplex H = H1 * H2.inv();
    return(H);
}

// overloaded += operator
Hypercomplex& Hypercomplex::operator+=(const Hypercomplex &H) {
    Hypercomplex result = (*this) + H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded -= operator
Hypercomplex& Hypercomplex::operator-=(const Hypercomplex &H) {
    Hypercomplex result = (*this) - H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded *= operator
Hypercomplex& Hypercomplex::operator*=(const Hypercomplex &H) {
    Hypercomplex result = (*this) * H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded ^= operator
Hypercomplex& Hypercomplex::operator^=(const unsigned int x) {
    Hypercomplex result = (*this) ^ x;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded /= operator
Hypercomplex& Hypercomplex::operator/=(const Hypercomplex &H) {
    Hypercomplex result = (*this) / H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overload << operator
std::ostream& operator<< (std::ostream &os, const Hypercomplex &H) {
    for (unsigned int i=0; i < H._() - 1; i++) os << H[i] << " ";
    os << H[H._() - 1];
    return os;
}

// return the real part of the number
Hypercomplex Re(const Hypercomplex &H) {
    Hypercomplex result = H;
    for (unsigned int i=1; i < H._(); i++) result[i] = 0;
    return result;
}

// return the imaginary part of the number
Hypercomplex Im(const Hypercomplex &H) {
    Hypercomplex result = H;
    result[0] = 0;
    return result;
}

// calculate e^H
Hypercomplex exp(const Hypercomplex &H) {
    unsigned int dim = H._();
    float *temparr = new float[1];
    temparr[0] = exp(H[0]);
    Hypercomplex term1 = Hypercomplex(1, temparr).expand(dim);
    Hypercomplex ImH = Im(H);
    temparr[0] = cos(ImH.norm());
    Hypercomplex term2 = Hypercomplex(1, temparr).expand(dim);
    temparr[0] = sin(ImH.norm()) / ImH.norm();
    Hypercomplex term3 = ImH * Hypercomplex(1, temparr).expand(dim);
    Hypercomplex result = term1 * (term2 + term3);
    delete[] temparr;
    return result;
}

// calculate ln(H)
Hypercomplex log(const Hypercomplex &H) {
    unsigned int dim = H._();
    float *temparr = new float[1];
    Hypercomplex result = H;
    result[0] = log(H.norm());
    Hypercomplex ImH = Im(H);
    temparr[0] = acos(H[0] / H.norm()) / ImH.norm();
    Hypercomplex product = Hypercomplex(1, temparr).expand(dim) * Im(H);
    for (unsigned int i=1; i < dim; i++) result[i] = product[i];
    delete[] temparr;
    return result;
}
