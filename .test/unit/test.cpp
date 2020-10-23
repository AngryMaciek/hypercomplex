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
#include "hypercomplex/hypercomplex.h"

int main(void){
    unsigned int x = 4;
    float A[] = {1.0, 2.0, 0.0, -1.0};

    Hypercomplex h = Hypercomplex(x, A);

    for (unsigned int i=0; i<x; i++) {
        std::cout << h.arr[i] << " ";
    } std::cout << std::endl;

    if (h._() != x){ std::abort(); }

    Hypercomplex h_ = !h;

    for (unsigned int i=0; i<x; i++) {
        std::cout << h_.arr[i] << " ";
    } std::cout << std::endl;

    if (h_.arr[0] != A[0]) {std::abort(); }
    if (h_.arr[1] != A[1]) {std::abort(); }
    if (h_.arr[2] != -A[2]) {std::abort(); }
    if (h_.arr[3] != -A[3]) {std::abort(); }

    return 0;
}
