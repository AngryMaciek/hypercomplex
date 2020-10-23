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

// overloaded operator
Hypercomplex Hypercomplex::operator! () {
    float * temparr = new float[d];
    for (unsigned int i=0; i < d/2; i++) {
        temparr[i] = arr[i];
    }
    for (unsigned int i=d/2; i < d; i++) {
        temparr[i] = -arr[i];
    }
    Hypercomplex H = Hypercomplex(d, temparr);
    delete[] temparr;
    return H;
}
