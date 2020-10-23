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
Hypercomplex::Hypercomplex(float v) {
    this->v=v;
}

// overloaded operator
Hypercomplex Hypercomplex::operator! () {
    return Hypercomplex(-v);
}
