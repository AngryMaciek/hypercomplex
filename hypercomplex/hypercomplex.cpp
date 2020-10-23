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
#include "hypercomplex.h" // NOLINT

// Hypercomplex constructor
Hypercomplex::Hypercomplex(unsigned int d, float* arr) {
    float * temparr = new float[d];
    for (unsigned int i=0; i < d; i++) {
        temparr[i] = arr[i];
    }
    this->d = d;
    this->arr = temparr;
}

// Hypercomplex destructor
Hypercomplex::~Hypercomplex() {
    delete[] arr;
}

// ~ operator
Hypercomplex Hypercomplex::operator~ () {
    float * temparr = new float[d];
    temparr[0] = arr[0];
    for (unsigned int i=1; i < d; i++) {
        temparr[i] = -arr[i];
    }
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// - unary operator
Hypercomplex Hypercomplex::operator- () {
    float * temparr = new float[d];
    for (unsigned int i=0; i < d; i++) {
        temparr[i] = -arr[i];
    }
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}

// overloaded == operator
bool Hypercomplex::operator ==(const Hypercomplex& H) {
    if (d != H.d) {
        return false;
    }
    for (unsigned int i=0; i < d; i++) {
        if (arr[i] != H.arr[i]) {
            return false;
        }
    }
    return true;
}

// overloaded != operator
bool Hypercomplex::operator !=(const Hypercomplex& H) {
    if (d != H.d) {
        return true;
    }
    for (unsigned int i=0; i < d; i++) {
        if (arr[i] != H.arr[i]) {
            return true;
        }
    }
    return false;
}

/*
Operators:
+ - += -=
* / *= /=
^ ^=
= (h1=h2 memory leak!)
*/
