/*
###############################################################################
#
#   Test: include installed library
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include <Hypercomplex/Hypercomplex.h>

int main(void){
    unsigned int x = 4;
    float A[] = {1.0, 2.0, 0.0, -1.0};
    Hypercomplex h = Hypercomplex(x, A);
    return 0;
}
