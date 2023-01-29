/*
###############################################################################
#
#   Test: code presented in the documentation & the paper
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
#include <cassert>
#include "Hypercomplex.hpp"

int main(void){
    //
    // DOCUMENTATION: main template code
    // 
    double arr1[4] = {1.0,0.0,-0.5,5.0};
    Hypercomplex<double, 4> H1(arr1);
    std::cout << "H1 = " << H1 << std::endl;
    double arr2[4] = {-2.0,-4.0,-6.0,0.0};
    Hypercomplex<double, 4> H2(arr2);
    //
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
    //
    // DOCUMENTATION: mpfr specialisation
    // 
    set_mpfr_precision(200);
    //
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
    //
    Hypercomplex<mpfr_t, 8> Hx(A);
    std::cout << "Hx^30 = ";
    mpfr_out_str(stdout, 10, 0, (Hx^30)[0], MPFR_RNDN);
    std::cout << std::endl;
    //
    mpfr_clear(A[0]);
    mpfr_clear(A[1]);
    mpfr_clear(A[2]);
    mpfr_clear(A[3]);
    mpfr_clear(A[4]);
    mpfr_clear(A[5]);
    mpfr_clear(A[6]);
    mpfr_clear(A[7]);
    clear_mpfr_memory();
    //
    // PUBLICATION: 7 magic numbers
    //
    unsigned int seedzero = 0;
    const unsigned int fig1a_dim = 64;
    const unsigned int fig1a_MaxDeg = 6; // N = 7
    const int64_t fig1a_p = 2;
    const int64_t fig1a_q = 1151;
    // Public Key
    Polynomial<fig1a_MaxDeg> F_coefficients[fig1a_dim];
    F_coefficients[1][0] = 0;
    F_coefficients[1][1] = 1;
    F_coefficients[1][2] = 0;
    F_coefficients[1][3] = 1;
    F_coefficients[1][4] = 1;
    F_coefficients[1][5] = 1;
    F_coefficients[1][6] = 1;
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> F(F_coefficients);
    CenteredLift(&F, fig1a_p);
    Polynomial<fig1a_MaxDeg> G_coefficients[fig1a_dim];
    for (unsigned int i=0; i < fig1a_dim; i++) {
        for (unsigned int j=0; j <= fig1a_MaxDeg; j++) {
            G_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> G(G_coefficients);
    CenteredLift(&G, fig1a_p);
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> H = PUBLICKEY(
        F, G, fig1a_q
    );
    // Encryption
    //
    // A: 72057611217901571
    // B: 1572864
    // C: 2308130010587988488
    // D: 16140901064495858689
    // E: 5
    // F: 9134425493490980523
    // G: 9007199254740992
    //
    // A B C D E F G
    Polynomial<fig1a_MaxDeg> M_coefficients[fig1a_dim];
    //
    // A
    M_coefficients[0][0] = 1;
    M_coefficients[1][0] = 1;
    M_coefficients[11][0] = 1;
    M_coefficients[12][0] = 1;
    M_coefficients[15][0] = 1;
    M_coefficients[16][0] = 1;
    M_coefficients[34][0] = 1;
    M_coefficients[56][0] = 1;
    // B
    M_coefficients[19][1] = 1;
    M_coefficients[20][1] = 1;
    // C
    M_coefficients[3][2] = 1;
    M_coefficients[9][2] = 1;
    M_coefficients[18][2] = 1;
    M_coefficients[23][2] = 1;
    M_coefficients[34][2] = 1;
    M_coefficients[45][2] = 1;
    M_coefficients[51][2] = 1;
    M_coefficients[61][2] = 1;
    // D
    M_coefficients[0][3] = 1;
    M_coefficients[10][3] = 1;
    M_coefficients[61][3] = 1;
    M_coefficients[62][3] = 1;
    M_coefficients[63][3] = 1;
    // E
    M_coefficients[0][4] = 1;
    M_coefficients[2][4] = 1;
    // F
    M_coefficients[0][5] = 1;
    M_coefficients[1][5] = 1;
    M_coefficients[3][5] = 1;
    M_coefficients[5][5] = 1;
    M_coefficients[7][5] = 1;
    M_coefficients[9][5] = 1;
    M_coefficients[10][5] = 1;
    M_coefficients[12][5] = 1;
    M_coefficients[14][5] = 1;
    M_coefficients[16][5] = 1;
    M_coefficients[18][5] = 1;
    M_coefficients[19][5] = 1;
    M_coefficients[22][5] = 1;
    M_coefficients[23][5] = 1;
    M_coefficients[25][5] = 1;
    M_coefficients[26][5] = 1;
    M_coefficients[27][5] = 1;
    M_coefficients[32][5] = 1;
    M_coefficients[33][5] = 1;
    M_coefficients[34][5] = 1;
    M_coefficients[36][5] = 1;
    M_coefficients[39][5] = 1;
    M_coefficients[40][5] = 1;
    M_coefficients[41][5] = 1;
    M_coefficients[42][5] = 1;
    M_coefficients[43][5] = 1;
    M_coefficients[44][5] = 1;
    M_coefficients[45][5] = 1;
    M_coefficients[46][5] = 1;
    M_coefficients[47][5] = 1;
    M_coefficients[48][5] = 1;
    M_coefficients[49][5] = 1;
    M_coefficients[54][5] = 1;
    M_coefficients[55][5] = 1;
    M_coefficients[57][5] = 1;
    M_coefficients[58][5] = 1;
    M_coefficients[59][5] = 1;
    M_coefficients[60][5] = 1;
    M_coefficients[61][5] = 1;
    M_coefficients[62][5] = 1;
    // G
    M_coefficients[53][6] = 1;
    //
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> M(M_coefficients);
    //
    Polynomial<fig1a_MaxDeg> PHI_coefficients[fig1a_dim];
    for (unsigned int i=0; i < fig1a_dim; i++) {
        for (unsigned int j=0; j <= fig1a_MaxDeg; j++) {
            PHI_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> PHI(PHI_coefficients);
    CenteredLift(&PHI, fig1a_p);
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> E = ENCRYPT(H, M, PHI, fig1a_p, fig1a_q);
    // Decryption
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> D = DECRYPT(F, E, fig1a_p, fig1a_q);
    CenteredLift(&M, fig1a_p);
    assert( D == M );
    //
    return 0;
}
