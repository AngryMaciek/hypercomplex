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

#include <mpfr.h>
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
    unsigned int _() const { return dim; }
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
    if constexpr (dim == 1) {
        T temparr[] = { H1[0] * H2[0] };
        Hypercomplex<T, 1> H_(temparr);
        return H_;
    // recursion step:
    } else {
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
    for (unsigned int i=1; i < dim; i++) result[i] = T();
    return result;
}

// return the imaginary part of the number
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> Im(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = H;
    result[0] = T();
    return result;
}

// calculate e^H
template <typename T, const unsigned int dim>
Hypercomplex<T, dim> exp(const Hypercomplex<T, dim> &H) {
    Hypercomplex<T, dim> result = Im(H);
    T zero = T();
    T norm = result.norm();
    if (norm == zero) {
        result[0] = exp(H[0]);
        for (unsigned int i=1; i < dim; i++) result[i] = zero;
    } else {
        T sinv_v = sin(norm) / norm;
        for (unsigned int i=0; i < dim; i++) result[i] = result[i] * sinv_v;
        result[0] = result[0] + cos(norm);
        for (unsigned int i=0; i < dim; i++) result[i] = result[i] * exp(H[0]);
    }
    return result;
}

/*
###############################################################################
#
#   Explicit template specialisation & function overloading for mpfr_t type
#
###############################################################################
*/

static unsigned int MPFR_global_precision;

unsigned int get_mpfr_precision() {
    return MPFR_global_precision;
}

void set_mpfr_precision(unsigned int n) {
    if (n>0) {
        MPFR_global_precision = n;
    } else {
        throw std::invalid_argument("invalid argument");
    }
}

void clear_mpfr_memory() {
    mpfr_free_cache();
    assert(!mpfr_mp_memory_cleanup());
}

template <const unsigned int dim>
class Hypercomplex<mpfr_t, dim> {

 private:
    mpfr_ptr arr;

 public:

    explicit Hypercomplex(const mpfr_ptr ARR) {
        if (dim == 0) throw std::invalid_argument("invalid dimension");
        if ((dim & (dim - 1)) != 0) {
            throw std::invalid_argument("invalid dimension");
        }
        arr = new mpfr_t[dim];
        for (unsigned int i=0; i < dim; i++) arr[i] = ARR[i];
    }

    Hypercomplex(const Hypercomplex &H) {
        arr = new mpfr_t[dim];
        for (unsigned int i=0; i < dim; i++) arr[i] = H[i];
    }

    Hypercomplex() = delete;

    ~Hypercomplex() {
        for (unsigned int i=0; i < dim; i++) mpfr_clear(arr[i]);
        delete[] arr;
    }

    unsigned int _() const { return dim; }

    mpfr_t norm() const {
        mpfr_t result, temp;
        mpfr_init2(result, MPFR_global_precision);
        mpfr_init2(temp, MPFR_global_precision);
        mpfr_set_zero(result, 0);
        for (unsigned int i=0; i < dim; i++) {
            mpfr_mul(temp, arr[i], arr[i], MPFR_RNDN);
            mpfr_add(result, result, temp, MPFR_RNDN);
        }
        mpfr_sqrt(result, result, MPFR_RNDN);
        mpfr_clear(temp);
        return result;
    }

    Hypercomplex inv() const {
        mpfr_t zero, norm;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_init2(norm, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        norm = (*this).norm();
        if (mpfr_equal_p(norm, zero)) {
            mpfr_clear(zero);
            mpfr_clear(norm);
            throw std::invalid_argument("division by zero");
        } else {
            mpfr_ptr temparr = new mpfr_t[dim];
            mpfr_mul(norm, norm, norm, MPFR_RNDN);
            mpfr_div(temparr[0], arr[0], norm, MPFR_RNDN);
            for (unsigned int i=1; i < dim; i++) {
                mpfr_div(temparr[i], arr[i], norm, MPFR_RNDN);
                mpfr_sub(temparr[i], zero, temparr[i], MPFR_RNDN);
            }
            Hypercomplex<mpfr_t, dim> H(temparr);
            mpfr_clear(zero);
            mpfr_clear(norm);
            for (unsigned int i=0; i < dim; i++) mpfr_clear(temparr[i]);
            delete[] temparr;
            return H;
        }
    }

    template <const unsigned int newdim>
    Hypercomplex<mpfr_t, newdim> expand() const {
        if (newdim <= dim) throw std::invalid_argument("invalid dimension");
        mpfr_ptr temparr = new mpfr_t[newdim];
        for (unsigned int i=0; i < dim; i++) temparr[i] = arr[i];
        for (unsigned int i=dim; i < newdim; i++) mpfr_set_zero(temparr[i], 0);
        Hypercomplex<mpfr_t, newdim> H(temparr);
        for (unsigned int i=0; i < newdim; i++) mpfr_clear(temparr[i]);
        delete[] temparr;
        return H;
    }

    Hypercomplex operator~ () const {
        mpfr_t zero;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        mpfr_ptr temparr = new mpfr_t[dim];
        temparr[0] = arr[0];
        for (unsigned int i=1; i < dim; i++)
            mpfr_sub(temparr[i], zero, arr[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (unsigned int i=0; i < newdim; i++) mpfr_clear(temparr[i]);
        delete[] temparr;
        mpfr_clear(zero);
        return H;
    }

    Hypercomplex operator- () const {
        mpfr_t zero;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        mpfr_ptr temparr = new mpfr_t[dim];
        for (unsigned int i=0; i < dim; i++)
            mpfr_sub(temparr[i], zero, arr[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (unsigned int i=0; i < newdim; i++) mpfr_clear(temparr[i]);
        delete[] temparr;
        mpfr_clear(zero);
        return H;
    }

    Hypercomplex& operator= (const Hypercomplex &H) {
        if (this == &H) return *this;
        for (unsigned int i=0; i < dim; i++) arr[i] = H[i];
        return *this;
    }

    mpfr_t& operator[] (const unsigned int i) const {
        assert(0 <= i && i < dim);
        return arr[i];
    }

    Hypercomplex& operator+= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) + H;
        for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
        return *this;
    }

    Hypercomplex& operator-= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) - H;
        for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
        return *this;
    }

    Hypercomplex& operator*= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) * H;
        for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
        return *this;
    }

    Hypercomplex& operator^= (const unsigned int x) {
        Hypercomplex<mpfr_t, dim> result = (*this) ^ x;
        for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
        return *this;
    }

    Hypercomplex& operator/= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) / H;
        for (unsigned int i=0; i < dim; i++) (*this)[i] = result[i];
        return *this;
    }

};

bool operator==(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    for (unsigned int i=0; i < dim; i++) {
        if (!mpfr_equal_p(H1[i], H2[i])) return false;
    }
    return true;
}

bool operator!=(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    return !(H1 == H2);
}

Hypercomplex<mpfr_t, dim> operator+(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    mpfr_ptr temparr = new mpfr_t[dim];
    for (unsigned int i=0; i < dim; i++)
        mpfr_add(temparr[i], H1[i], H2[i], MPFR_RNDN);
    Hypercomplex<mpfr_t, dim> H(temparr);
    for (unsigned int i=0; i < dim; i++) mpfr_clear(temparr[i]);
    delete[] temparr;
    return H;
}

Hypercomplex<mpfr_t, dim> operator-(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    mpfr_ptr temparr = new mpfr_t[dim];
    for (unsigned int i=0; i < dim; i++)
        mpfr_sub(temparr[i], H1[i], H2[i], MPFR_RNDN);
    Hypercomplex<mpfr_t, dim> H(temparr);
    for (unsigned int i=0; i < dim; i++) mpfr_clear(temparr[i]);
    delete[] temparr;
    return H;
}

Hypercomplex<mpfr_t, dim> operator*(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    // recursion base:
    if constexpr (dim == 1) {
        mpfr_t result;
        mpfr_init2(result, MPFR_global_precision);
        mpfr_mul(result, H1[0], H2[0], MPFR_RNDN);
        mpfr_t temparr[] = { result };
        Hypercomplex<mpfr_t, 1> H_(temparr);
        mpfr_clear(result);
        return H_;
    // recursion step:
    } else {
        // shared objects:
        const unsigned int halfd = dim / 2;
        mpfr_ptr temparr = new mpfr_t[dim];
        // construct helper objects:
        for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i];
        Hypercomplex<mpfr_t, halfd> H1a(temparr);
        for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i+halfd];
        Hypercomplex<mpfr_t, halfd> H1b(temparr);
        for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i];
        Hypercomplex<mpfr_t, halfd> H2a(temparr);
        for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i+halfd];
        Hypercomplex<mpfr_t, halfd> H2b(temparr);
        // multiply recursively:
        Hypercomplex<mpfr_t, halfd> H1a2a = H1a * H2a;
        Hypercomplex<mpfr_t, halfd> H2b_1b = ~H2b * H1b;
        Hypercomplex<mpfr_t, halfd> H2b1a = H2b * H1a;
        Hypercomplex<mpfr_t, halfd> H1b2a_ = H1b * ~H2a;
        // construct the final object
        Hypercomplex<mpfr_t, halfd> Ha = H1a2a - H2b_1b;
        Hypercomplex<mpfr_t, halfd> Hb = H2b1a + H1b2a_;
        for (unsigned int i=0; i < halfd; i++) temparr[i] = Ha[i];
        for (unsigned int i=0; i < halfd; i++) temparr[i+halfd] = Hb[i];
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (unsigned int i=0; i < dim; i++) mpfr_clear(temparr[i]);
        delete[] temparr;
        return H;
    }
}

Hypercomplex<mpfr_t, dim> operator^(
    const Hypercomplex<mpfr_t, dim> &H,
    const unsigned int x
) {
    if (!(x)) {
        throw std::invalid_argument("zero is not a valid argument");
    } else {
        Hypercomplex<mpfr_t, dim> Hx(H);
        for (unsigned int i=0; i < x-1; i++) Hx = Hx * H;
        return Hx;
    }
}

Hypercomplex<mpfr_t, dim> operator/(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    Hypercomplex<mpfr_t, dim> H = H1 * H2.inv();
    return(H);
}

std::ostream& operator<<(
    std::ostream &os,
    const Hypercomplex<mpfr_t, dim> &H
) {
    for (unsigned int i=0; i < dim - 1; i++) {
        mpfr_out_str(os, 10, 0, H[i], MPFR_RNDN);
        os << " ";
    }
    mpfr_out_str(os, 10, 0, H[dim - 1]], MPFR_RNDN);
    return os;
}

Hypercomplex<mpfr_t, dim> Re(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = H;
    for (unsigned int i=1; i < dim; i++) mpfr_set_zero(result[i], 0);
    return result;
}

Hypercomplex<mpfr_t, dim> Im(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = H;
    mpfr_set_zero(result[0], 0);
    return result;
}

Hypercomplex<mpfr_t, dim> exp(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = Im(H);
    mpfr_t zero, norm, expreal;
    mpfr_init2(zero, MPFR_global_precision);
    mpfr_init2(norm, MPFR_global_precision);
    mpfr_init2(expreal, MPFR_global_precision);
    mpfr_set_zero(zero, 0);
    norm = result.norm();
    mpfr_exp(expreal, H[0], MPFR_RNDN);

    if (mpfr_equal_p(norm, zero)) {
        result[0] = expreal;
        for (unsigned int i=1; i < dim; i++) mpfr_set_zero(result[i], 0);
    } else {
        mpfr_t sinv_v;
        mpfr_init2(sinv_v, MPFR_global_precision);
        mpfr_sin(sinv_v, norm, MPFR_RNDN);
        mpfr_div(sinv_v, sinv_v, norm, MPFR_RNDN);
        for (unsigned int i=0; i < dim; i++) {
            mpfr_mul(result[i], result[i], sinv_v, MPFR_RNDN)
        }
        mpfr_cos(norm, norm, MPFR_RNDN);
        mpfr_add(result[0], result[0], norm, MPFR_RNDN);
        for (unsigned int i=0; i < dim; i++) {
            mpfr_mul(result[i], result[i], expreal, MPFR_RNDN);
        }
    }
    mpfr_clear(zero);
    mpfr_clear(norm);
    mpfr_clear(expreal);
    return result;
}


// function partial specialisation
// https://stackoverflow.com/questions/3768862/c-single-template-specialisation-with-multiple-template-parameters/3769025

// delete arr from constrctor?
// expand: () after array
// init2 for all temparr[i]
// assignment operator between mpfr_t?

// include at the end, check docs why
// last paragraph 4.7
// https://www.mpfr.org/mpfr-current/mpfr.pdf
// 1 function to set global precision at the start
// 1 function to clear all cache annd mem at the end

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_HPP_
