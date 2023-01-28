/*
###############################################################################
#
#   Test: code presented in the documentation
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Department_of_Mathematics_City_University_of_London
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

#include <iostream>
#include "Hypercomplex.hpp"

int main(void){
    //
    double arr1[4] = {1.0,0.0,-0.5,5.0};
    Hypercomplex<double, 4> H1(arr1);
    std::cout << "H1 = " << H1 << std::endl;
    
    double arr2[4] = {-2.0,-4.0,-6.0,0.0};
    Hypercomplex<double, 4> H2(arr2);
    std::cout << "H2 = " << H2 << std::endl;

    std::cout << "Re(H1) = " << Re(H1) << std::endl;
    std::cout << "Im(H1) = " << Im(H1) << std::endl;

    std::cout << "dim(H1) = " << H1._() << std::endl;

    std::cout << "H2(2) = " << H2[2] << std::endl;

    std::cout << "Oct(H1) = " << H1.expand<8>() << std::endl;

    std::cout << "||H2|| = " << H2.norm() << std::endl;

    std::cout << "H2^-1 = " << H2.inv() << std::endl;

    std::cout << "H1 + H2 = " << H1 + H2 << std::endl;
    std::cout << "H1 - H2 = " << H1 - H2 << std::endl;
    std::cout << "H1 * H2 = " << H1 * H2 << std::endl;
    std::cout << "H1 / H2 = " << H1 / H2 << std::endl;

    std::cout << "H2^4 = " << (H2^4) << std::endl;

    std::cout << "e^H1 = " << exp(H1) << std::endl;

    set_mpfr_precision(200);

    mpfr_t A[8];
    mpfr_init2(A[0], get_mpfr_precision());
    mpfr_init2(A[1], get_mpfr_precision());
    mpfr_init2(A[2], get_mpfr_precision());
    mpfr_init2(A[3], get_mpfr_precision());
    mpfr_init2(A[4], get_mpfr_precision());
    mpfr_init2(A[5], get_mpfr_precision());
    mpfr_init2(A[6], get_mpfr_precision());
    mpfr_init2(A[7], get_mpfr_precision());
    mpfr_set_d(A[0], 1.5, MPFR_RNDN);
    mpfr_set_d(A[1], 2.5, MPFR_RNDN);
    mpfr_set_d(A[2], 0.0, MPFR_RNDN);
    mpfr_set_d(A[3], -1.5, MPFR_RNDN);
    mpfr_set_d(A[4], 0.5, MPFR_RNDN);
    mpfr_set_d(A[5], -0.5, MPFR_RNDN);
    mpfr_set_d(A[6], -0.5, MPFR_RNDN);
    mpfr_set_d(A[7], -1.5, MPFR_RNDN);

    Hypercomplex<mpfr_t, 8> Hx(A);

    std::cout << "Hx^30 = ";
    mpfr_out_str(stdout, 10, 0, (Hx^30)[0], MPFR_RNDN);
    std::cout << std::endl;

    mpfr_clear(A[0]);
    mpfr_clear(A[1]);
    mpfr_clear(A[2]);
    mpfr_clear(A[3]);
    mpfr_clear(A[4]);
    mpfr_clear(A[5]);
    mpfr_clear(A[6]);
    mpfr_clear(A[7]);

    clear_mpfr_memory();

    return 0;
}
