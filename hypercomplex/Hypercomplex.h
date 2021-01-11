// Copyright 2020 <Maciej Bak>
/*
###############################################################################
#
#   Hypercomplex header-only library.
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

// mark that we included this library
#ifndef HYPERCOMPLEX_HYPERCOMPLEX_H_
#define HYPERCOMPLEX_HYPERCOMPLEX_H_

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

/*
###############################################################################
#
#   Header section
#
###############################################################################
*/

// Main class of the library
template <typename T, const unsigned int dim>
class Hypercomplex {
 private:
    unsigned int d;
    T* arr;
 public:
    explicit Hypercomplex(const T* ARR);
    Hypercomplex(const Hypercomplex& H);
    Hypercomplex() = delete;  // forbid default constructor | c++11
    ~Hypercomplex();
    T _() const { return d; }
    T norm() const;
    Hypercomplex inv() const;
    Hypercomplex expand(const unsigned int D) const;
    Hypercomplex operator~ () const;
    Hypercomplex operator- () const;
    Hypercomplex& operator= (const Hypercomplex &H);
    T& operator[] (const unsigned int i) const;
    Hypercomplex& operator+= (const Hypercomplex &H);
    Hypercomplex& operator-= (const Hypercomplex &H);
    Hypercomplex& operator*= (const Hypercomplex &H);
    Hypercomplex& operator^= (const unsigned int x);
    Hypercomplex& operator/= (const Hypercomplex &H);
};

// Operators
bool operator== (const Hypercomplex &H1, const Hypercomplex &H2);
bool operator!= (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator+ (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator- (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator* (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator^ (const Hypercomplex &H, const unsigned int x);
Hypercomplex operator/ (const Hypercomplex &H1, const Hypercomplex &H2);
std::ostream& operator<< (std::ostream &os, const Hypercomplex &H);


// Global functions
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Re(const Hypercomplex<T, dim> &H);
Hypercomplex Im(const Hypercomplex &H);
Hypercomplex exp(const Hypercomplex &H);

/*
###############################################################################
#
#   Implementation section
#
###############################################################################
*/

// Hypercomplex main constructor
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>::Hypercomplex(const T* ARR) {
    if (dim == 0) throw std::invalid_argument("invalid dimension");
    if ((dim & (dim - 1)) != 0) {
        throw std::invalid_argument("invalid dimension");
    }
    d = dim;
    arr = new T[dim];
    for (unsigned int i=0; i < dim; i++) arr[i] = ARR[i];
}

// Hypercomplex copy constructor
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>::Hypercomplex(const Hypercomplex& H) {
    d = H.d;
    arr = new T[H.d];
    for (unsigned int i=0; i < H.d; i++) arr[i] = H.arr[i];
}

// Hypercomplex destructor
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>::~Hypercomplex() {
    delete[] arr;
}

// calculate norm of the number
template <typename T, const unsigned int dim>
inline T Hypercomplex<T, dim>::norm() const {
    T result = 0.0;
    for (unsigned int i=0; i < d; i++) result += arr[i] * arr[i];
    return sqrt(result);
}

// calculate inverse of the number
template <typename T, const unsigned int dim>
Hypercomplex Hypercomplex<T, dim>::inv() const {
    if ((*this).norm() == 0) {
        throw std::invalid_argument("division by zero");
    } else {
        T norm2 = pow((*this).norm(), 2);
        T* temparr = new T[d];
        temparr[0] = arr[0] / norm2;
        for (unsigned int i=1; i < d; i++) temparr[i] = -arr[i] / norm2;
        Hypercomplex H = Hypercomplex(d, temparr);
        delete[] temparr;
        return H;
    }
}

// expand object to a higher dimension
template <typename T, const unsigned int dim>
Hypercomplex Hypercomplex<T, dim>::expand(const unsigned int D) const {
    if (D <= d) throw std::invalid_argument("invalid dimension");
    T* temparr = new T[D]();  // zero-init
    for (unsigned int i=0; i < d; i++) temparr[i] = arr[i];
    Hypercomplex H = Hypercomplex(D, temparr);
    delete[] temparr;
    return H;
}

// overloaded ~ operator
template <typename T, const unsigned int dim>
inline Hypercomplex Hypercomplex<T, dim>::operator~() const {
    Hypercomplex H = -(*this);
    H[0] = arr[0];
    return H;
}

// overloaded - unary operator
template <typename T, const unsigned int dim>
Hypercomplex Hypercomplex<T, dim>::operator-() const {
    T* temparr = new T[d];
    for (unsigned int i=0; i < d; i++) temparr[i] = -arr[i];
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded = operator
template <typename T, const unsigned int dim>
inline Hypercomplex& Hypercomplex<T, dim>::operator=(const Hypercomplex &H) {
    // self-assignment guard
    if (this == &H) return *this;
    // reassign
    d = H.d;
    for (unsigned int i=0; i < d; i++) arr[i] = H.arr[i];
    // return the existing object so we can chain this operator
    return *this;
}

// overloaded [] operator
template <typename T, const unsigned int dim>
inline T& Hypercomplex<T, dim>::operator[](const unsigned int i) const {
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
    T *temparr = new T[d];
    for (unsigned int i=0; i < d; i++) temparr[i] = H1[i] + H2[i];
    Hypercomplex H = Hypercomplex(H1._(), temparr);
    delete[] temparr;
    return H;
}

// overloaded - binary operator
Hypercomplex operator-(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) throw std::invalid_argument("operand mismatch");
    unsigned int d = H1._();
    T* temparr = new T[d];
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
        T temparr[] = { H1[0] * H2[0] };
        return Hypercomplex(1, temparr);
    }
    // shared objects:
    unsigned int halfd = d / 2;
    T* temparr = new T[d];
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
template <typename T, const unsigned int dim>
Hypercomplex& Hypercomplex<T, dim>::operator+=(const Hypercomplex &H) {
    Hypercomplex result = (*this) + H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded -= operator
template <typename T, const unsigned int dim>
Hypercomplex& Hypercomplex<T, dim>::operator-=(const Hypercomplex &H) {
    Hypercomplex result = (*this) - H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded *= operator
template <typename T, const unsigned int dim>
Hypercomplex& Hypercomplex<T, dim>::operator*=(const Hypercomplex &H) {
    Hypercomplex result = (*this) * H;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded ^= operator
template <typename T, const unsigned int dim>
Hypercomplex& Hypercomplex<T, dim>::operator^=(const unsigned int x) {
    Hypercomplex result = (*this) ^ x;
    for (unsigned int i=0; i < d; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded /= operator
template <typename T, const unsigned int dim>
Hypercomplex& Hypercomplex<T, dim>::operator/=(const Hypercomplex &H) {
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
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Re(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = H;
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
    Hypercomplex result = Im(H);
    T norm = result.norm();
    if (norm == 0.0) {
        result[0] = exp(H[0]);
        for (unsigned int i=1; i < dim; i++) result[i] = 0;
    } else {
        T sinv_v = sin(norm) / norm;
        for (unsigned int i=0; i < dim; i++) result[i] *= sinv_v;
        result[0] += cos(norm);
        for (unsigned int i=0; i < dim; i++) result[i] *= exp(H[0]);
    }
    return result;
}

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_H_
