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
    this->d = d;
    this->arr = arr;
}

// overloaded operator
Hypercomplex Hypercomplex::operator! () {
    float temparr[d];
    for (unsigned int i=0; i<d; i++) {
        temparr[i] = -arr[i];
    }
    Hypercomplex H = Hypercomplex(d, temparr);
    return H;
}
