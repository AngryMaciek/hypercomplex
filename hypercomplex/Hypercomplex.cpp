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

#include <cstdlib>
#include <cassert>
#include <stdexcept>
#include "Hypercomplex.h" // NOLINT

// Hypercomplex main constructor
Hypercomplex::Hypercomplex(unsigned int arg_d, float* arg_arr) {
    if (arg_d == 0) {
        throw std::invalid_argument("zero is not a valid argument");
    }
    d = arg_d;
    arr = new float[arg_d];
    for (unsigned int i=0; i < arg_d; i++) {
        arr[i] = arg_arr[i];
    }
}

// Hypercomplex copy constructor
Hypercomplex::Hypercomplex(const Hypercomplex& H) {
    d = H.d;
    arr = new float[H.d];
    for (unsigned int i=0; i < H.d; i++) {
        arr[i] = H.arr[i];
    }
}

// Hypercomplex destructor
Hypercomplex::~Hypercomplex() {
    delete[] arr;
}

// overloaded ~ operator
Hypercomplex Hypercomplex::operator~() const {
    float *temparr = new float[d];
    temparr[0] = arr[0];
    for (unsigned int i=1; i < d; i++) {
        temparr[i] = -arr[i];
    }
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded - unary operator
Hypercomplex Hypercomplex::operator-() const {
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) {
        temparr[i] = -arr[i];
    }
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
    for (unsigned int i=0; i < d; i++) {
        arr[i] = H.arr[i];
    }
    // return the existing object so we can chain this operator
    return *this;
}

// overloaded [] operator
float& Hypercomplex::operator[](unsigned int i) {
    assert(0 <= i && i < d);
    return arr[i];
}

// overloaded [] operator
const float& Hypercomplex::operator[](unsigned int i) const {
    assert(0 <= i && i < d);
    return arr[i];
}

// overloaded == operator
bool operator==(const Hypercomplex &H1, const Hypercomplex &H2) {
    if (H1._() != H2._()) {
        return false;
    }
    for (unsigned int i=0; i < H1._(); i++) {
        if (H1[i] != H2[i]) {
            return false;
        }
    }
    return true;
}

// overloaded != operator
bool operator!=(const Hypercomplex &H1, const Hypercomplex &H2) {
    return !(H1 == H2);
}

// overloaded + binary operator
Hypercomplex operator+(const Hypercomplex &H1, const Hypercomplex &H2) {
    unsigned int d = H1._();
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) {
        temparr[i] = H1[i] + H2[i];
    }
    Hypercomplex H = Hypercomplex(H1._(), temparr);
    delete[] temparr;
    return H;
}

// overloaded - binary operator
Hypercomplex operator-(const Hypercomplex &H1, const Hypercomplex &H2) {
    unsigned int d = H1._();
    float *temparr = new float[d];
    for (unsigned int i=0; i < d; i++) {
        temparr[i] = H1[i] - H2[i];
    }
    Hypercomplex H = Hypercomplex(H1._(), temparr);
    delete[] temparr;
    return H;
}

// overloaded * binary operator
Hypercomplex operator*(const Hypercomplex &H1, const Hypercomplex &H2) {
    unsigned int d = H1._();
    // recursion base:
    if (d==1) {
        float temparr[] = { H1[0] * H2[0] };
        return Hypercomplex(1, temparr);
    }
    // shared objects:
    unsigned int halfd = d / 2;
    float *temparr;
    temparr = new float[halfd];
    // construct helper objects:
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H1[i];
    Hypercomplex H1a = Hypercomplex(halfd, temparr);
    for (unsigned int i=halfd; i < d; i++) temparr[i] = H1[i];
    Hypercomplex H1b = Hypercomplex(halfd, temparr);
    for (unsigned int i=0; i < halfd; i++) temparr[i] = H2[i];
    Hypercomplex H2a = Hypercomplex(halfd, temparr);
    for (unsigned int i=halfd; i < d; i++) temparr[i] = H2[i];
    Hypercomplex H2b = Hypercomplex(halfd, temparr);
    delete[] temparr;
    // multiply recursively:
    Hypercomplex H1a1b = H1a * H1b;
    Hypercomplex H1a2b = H1a * H2b;
    Hypercomplex H2a1b = H2a * H1b;
    Hypercomplex H2a2b = H2a * H2b;
    // helper addition/subtraction
    Hypercomplex Ha = H1a2b - H2a1b;
    Hypercomplex Hb = H1a1b + H2a2b;
    // construct the final object
    temparr = new float[d];
    for (unsigned int i=0; i < halfd; i++) temparr[i] = Ha[i];
    for (unsigned int i=0; i < halfd; i++) temparr[i+halfd] = Hb[i];
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded += operator
Hypercomplex& Hypercomplex::operator+=(const Hypercomplex &H) {
    for (unsigned int i=0; i < d; i++) {
        (*this)[i] = (*this)[i] + H[i];
    }
    return *this;
}

// overloaded -= operator
Hypercomplex& Hypercomplex::operator-=(const Hypercomplex &H) {
    for (unsigned int i=0; i < d; i++) {
        (*this)[i] = (*this)[i] - H[i];
    }
    return *this;
}
