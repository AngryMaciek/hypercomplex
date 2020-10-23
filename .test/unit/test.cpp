/*
###############################################################################
#
#   Test: unit tests
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 23-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include <iostream>
#include "hypercomplex.h"

int main(void){

    if (!mirror(1)){
        std::abort();
    }
    if (mirror(0)){
        std::abort();
    }

    float x = 1.2;
    Hypercomplex h = Hypercomplex(x);
    if (h._() != x){ std::abort(); }
    if ( (!h)._() != -x){ std::abort(); }

    std::cout << (!h)._() << std::endl;

    return 0;
}
