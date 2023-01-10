/*
###############################################################################
#
#   Test: raw library include
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Department_of_Mathematics_City_University_of_London
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include "hypercomplex/Hypercomplex.hpp"

int main(void){
    double A[] = {1.0, 2.0, 0.0, -1.0};
    Hypercomplex<double, 4> h(A);
    return 0;
}
