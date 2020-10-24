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

#include "test.h"
#include <iostream>
#include "hypercomplex/hypercomplex.h"


TEST_CASE( "DemoTest", "Demo" ) {
REQUIRE( 2 == 1 );
}

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}



/*int main(void){
    unsigned int dim = 4;
    float A[] = {1.0, 2.0, 0.0, -1.0};
    float B[] = {-0.5, 1.0, 0.0, 6.0};

    Hypercomplex h1 = Hypercomplex(dim, A);
    Hypercomplex h2 = Hypercomplex(dim, B);

    for (unsigned int i=0; i<dim; i++) {
        std::cout << h1.arr[i] << " ";
    } std::cout << std::endl;

    // test attribute getter
    if (h1._() != dim){ std::abort(); }

    // test ~ operator
    Hypercomplex h1_ = ~h1;

    for (unsigned int i=0; i<dim; i++) {
        std::cout << h1_.arr[i] << " ";
    } std::cout << std::endl;

    if (h1_.arr[0] != A[0]) {std::abort(); }
    if (h1_.arr[1] != -A[1]) {std::abort(); }
    if (h1_.arr[2] != -A[2]) {std::abort(); }
    if (h1_.arr[3] != -A[3]) {std::abort(); }

    // test == and != operators
    if (h1 == h2){ std::abort(); }
    if (!(h1 != h2)){ std::abort(); }

    // test unary - operator
    // ?

    // test [] operator
    if (h1[0]!=1.0) {std::abort(); }
    if (h1[1]!=2.0) {std::abort(); }
    if (h1[2]!=0.0) {std::abort(); }
    if (h1[3]!=-1.0) {std::abort(); }

    return 0;
}
*/