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
#ifndef HYPERCOMPLEX_HYPERCOMPLEX_HPP_
#define HYPERCOMPLEX_HYPERCOMPLEX_HPP_

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
    T* arr;
 public:
    explicit Hypercomplex(const T* ARR);
    Hypercomplex(const Hypercomplex &H);
    Hypercomplex() = delete;  // forbid default constructor | c++11
    ~Hypercomplex();
    T _() const { return dim; }
    T norm() const;
    Hypercomplex inv() const;
    template <const unsigned int newdim>
    Hypercomplex<T, newdim> expand() const;
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
template <typename T, const unsigned int dim>
bool operator== (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
bool operator!= (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator+ (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator- (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator* (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator^ (
    const Hypercomplex<T, dim> &H,
    const unsigned int x
);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator/ (
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
);
template <typename T, const unsigned int dim>
std::ostream& operator<< (std::ostream &os, const Hypercomplex<T, dim> &H);


// Global functions
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Re(const Hypercomplex<T, dim> &H);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Im(const Hypercomplex<T, dim> &H);
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> exp(const Hypercomplex<T, dim> &H);

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
    arr = new T[dim];
    for (unsigned int i=0; i < dim; i++) arr[i] = ARR[i];
}

// Hypercomplex copy constructor
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>::Hypercomplex(const Hypercomplex<T, dim> &H) {
    arr = new T[dim];
    for (unsigned int i=0; i < dim; i++) arr[i] = H[i];
}

// Hypercomplex destructor
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>::~Hypercomplex() {
    delete[] arr;
}

// calculate norm of the number
template <typename T, const unsigned int dim>
inline T Hypercomplex<T, dim>::norm() const {
    T result = T();
    for (unsigned int i=0; i < dim; i++) result = result + arr[i] * arr[i];
    return sqrt(result);
}

// calculate inverse of the number
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Hypercomplex<T, dim>::inv() const {
    T zero = T();
    T norm = (*this).norm();
    if (norm == zero) {
        throw std::invalid_argument("division by zero");
    } else {
        T* temparr = new T[dim];
        temparr[0] = arr[0] / (norm * norm);
        for (unsigned int i=1; i < dim; i++)
            temparr[i] = -arr[i] / (norm * norm);
        Hypercomplex<T, dim> H(temparr);
        delete[] temparr;
        return H;
    }
}

// cast object to a higher dimension
template <typename T, const unsigned int dim>
template <const unsigned int newdim>
Hypercomplex<T, newdim> Hypercomplex<T, dim>::expand() const {
    if (newdim <= dim) throw std::invalid_argument("invalid dimension");
    T* temparr = new T[newdim]();  // zero-init
    for (unsigned int i=0; i < dim; i++) temparr[i] = arr[i];
    Hypercomplex<T, newdim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded ~ operator
template <typename T, const unsigned int dim>
inline Hypercomplex<T, dim> Hypercomplex<T, dim>::operator~() const {
    T* temparr = new T[dim];
    temparr[0] = arr[0];
    for (unsigned int i=1; i < dim; i++) temparr[i] = -arr[i];
    Hypercomplex<T, dim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded - unary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Hypercomplex<T, dim>::operator-() const {
    T* temparr = new T[dim];
    for (unsigned int i=0; i < dim; i++) temparr[i] = -arr[i];
    Hypercomplex<T, dim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded = operator
template <typename T, const unsigned int dim>
inline Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator=(
    const Hypercomplex &H
) {
    // self-assignment guard
    if (this == &H) return *this;
    // reassign
    for (unsigned int i=0; i < dim; i++) arr[i] = H[i];
    // return the existing object so we can chain this operator
    return *this;
}

// overloaded [] operator
template <typename T, const unsigned int dim>
inline T& Hypercomplex<T, dim>::operator[](const unsigned int i) const {
    assert(0 <= i && i < dim);
    return arr[i];
}

// overloaded == operator
template <typename T, const unsigned int dim>
bool operator==(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    for (unsigned int i=0; i < dim; i++) {
        if (H1[i] != H2[i]) return false;
    }
    return true;
}

// overloaded != operator
template <typename T, const unsigned int dim>
bool operator!=(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    return !(H1 == H2);
}

// overloaded + binary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator+(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    T *temparr = new T[dim];
    for (unsigned int i=0; i < dim; i++) temparr[i] = H1[i] + H2[i];
    Hypercomplex<T, dim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded - binary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator-(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    T* temparr = new T[dim];
    for (unsigned int i=0; i < dim; i++) temparr[i] = H1[i] - H2[i];
    Hypercomplex<T, dim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded * binary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator*(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    // recursion base:
    if (dim == 1) {
        T temparr[] = { H1[0] * H2[0] };
        return Hypercomplex<T, 1>(temparr);
    }
    // shared objects:
    const unsigned int halfd = dim / 2;
    T* temparr = new T[dim];
    // construct helper objects:
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i];
    Hypercomplex<T, halfd> H1a(temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i+halfd];
    Hypercomplex<T, halfd> H1b(temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i];
    Hypercomplex<T, halfd> H2a(temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i+halfd];
    Hypercomplex<T, halfd> H2b(temparr);
    // multiply recursively:
    Hypercomplex<T, halfd> H1a2a = H1a * H2a;
    Hypercomplex<T, halfd> H2b_1b = ~H2b * H1b;
    Hypercomplex<T, halfd> H2b1a = H2b * H1a;
    Hypercomplex<T, halfd> H1b2a_ = H1b * ~H2a;
    // construct the final object
    Hypercomplex<T, halfd> Ha = H1a2a - H2b_1b;
    Hypercomplex<T, halfd> Hb = H2b1a + H1b2a_;
    for (unsigned int i=0; i < halfd; i++) temparr[i] = Ha[i];
    for (unsigned int i=0; i < halfd; i++) temparr[i+halfd] = Hb[i];
    Hypercomplex<T, dim> H(temparr);
    delete[] temparr;
    return H;
}

// overloaded ^ binary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator^(
    const Hypercomplex<T, dim> &H,
    const unsigned int x
) {
    if (!(x)) {
        throw std::invalid_argument("zero is not a valid argument");
    } else {
        Hypercomplex<T, dim> Hx(H);
        for (unsigned int i=0; i < x-1; i++) Hx = Hx * H;
        return Hx;
    }
}

// overloaded / binary operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> operator/(
    const Hypercomplex<T, dim> &H1,
    const Hypercomplex<T, dim> &H2
) {
    // division H1 / H2 is implemented as H1 * 1/H2
    Hypercomplex<T, dim> H = H1 * H2.inv();
    return(H);
}

// overloaded += operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator+=(
    const Hypercomplex<T, dim> &H
) {
    Hypercomplex<T, dim> result = (*this) + H;
    for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded -= operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator-=(
    const Hypercomplex<T, dim> &H
) {
    Hypercomplex<T, dim> result = (*this) - H;
    for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded *= operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator*=(
    const Hypercomplex<T, dim> &H
) {
    Hypercomplex<T, dim> result = (*this) * H;
    for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded ^= operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator^=(
    const unsigned int x
) {
    Hypercomplex<T, dim> result = (*this) ^ x;
    for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
    return *this;
}

// overloaded /= operator
template <typename T, const unsigned int dim>
Hypercomplex<T, dim>& Hypercomplex<T, dim>::operator/=(
    const Hypercomplex<T, dim> &H
) {
    Hypercomplex<T, dim> result = (*this) / H;
    for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
    return *this;
}

// overload << operator
template <typename T, const unsigned int dim>
std::ostream& operator<< (std::ostream &os, const Hypercomplex<T, dim> &H) {
    for (unsigned int i=0; i < dim - 1; i++) os << H[i] << " ";
    os << H[dim - 1];
    return os;
}

// return the real part of the number
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Re(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = H;
    for (unsigned int i=1; i < dim; i++) result[i] = 0;
    return result;
}

// return the imaginary part of the number
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Im(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = H;
    result[0] = 0;
    return result;
}

// calculate e^H
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> exp(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = Im(H);
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

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_HPP_
