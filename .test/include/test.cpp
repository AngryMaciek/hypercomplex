/*
###############################################################################
#
#   Test: raw library include
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include "hypercomplex/Hypercomplex.h"

int main(void){
    double A[] = {1.0, 2.0, 0.0, -1.0};
    Hypercomplex<double, 4> h(A);
    return 0;
}
