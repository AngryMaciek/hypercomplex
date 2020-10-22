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

#include <hypercomplex.h>

int main(void){
    if (!mirror(1)){
        return -1;
    }
    if (mirror(0)){
        return -1;
    }
    return 0;
}
