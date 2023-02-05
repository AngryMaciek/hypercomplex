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
#include<fstream>
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
    Polynomial<fig1a_MaxDeg> F1a_coefficients[fig1a_dim];
    F1a_coefficients[1][0] = 0;
    F1a_coefficients[1][1] = 1;
    F1a_coefficients[1][2] = 0;
    F1a_coefficients[1][3] = 1;
    F1a_coefficients[1][4] = 1;
    F1a_coefficients[1][5] = 1;
    F1a_coefficients[1][6] = 1;
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> F1a(F1a_coefficients);
    CenteredLift(&F1a, fig1a_p);
    Polynomial<fig1a_MaxDeg> G1a_coefficients[fig1a_dim];
    for (unsigned int i=0; i < fig1a_dim; i++) {
        for (unsigned int j=0; j <= fig1a_MaxDeg; j++) {
            G1a_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> G1a(G1a_coefficients);
    CenteredLift(&G1a, fig1a_p);
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> H1a = PUBLICKEY(
        F1a, G1a, fig1a_q
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
    Polynomial<fig1a_MaxDeg> M1a_coefficients[fig1a_dim];
    //
    // A
    M1a_coefficients[0][0] = 1;
    M1a_coefficients[1][0] = 1;
    M1a_coefficients[11][0] = 1;
    M1a_coefficients[12][0] = 1;
    M1a_coefficients[15][0] = 1;
    M1a_coefficients[16][0] = 1;
    M1a_coefficients[34][0] = 1;
    M1a_coefficients[56][0] = 1;
    // B
    M1a_coefficients[19][1] = 1;
    M1a_coefficients[20][1] = 1;
    // C
    M1a_coefficients[3][2] = 1;
    M1a_coefficients[9][2] = 1;
    M1a_coefficients[18][2] = 1;
    M1a_coefficients[23][2] = 1;
    M1a_coefficients[34][2] = 1;
    M1a_coefficients[45][2] = 1;
    M1a_coefficients[51][2] = 1;
    M1a_coefficients[61][2] = 1;
    // D
    M1a_coefficients[0][3] = 1;
    M1a_coefficients[10][3] = 1;
    M1a_coefficients[61][3] = 1;
    M1a_coefficients[62][3] = 1;
    M1a_coefficients[63][3] = 1;
    // E
    M1a_coefficients[0][4] = 1;
    M1a_coefficients[2][4] = 1;
    // F
    M1a_coefficients[0][5] = 1;
    M1a_coefficients[1][5] = 1;
    M1a_coefficients[3][5] = 1;
    M1a_coefficients[5][5] = 1;
    M1a_coefficients[7][5] = 1;
    M1a_coefficients[9][5] = 1;
    M1a_coefficients[10][5] = 1;
    M1a_coefficients[12][5] = 1;
    M1a_coefficients[14][5] = 1;
    M1a_coefficients[16][5] = 1;
    M1a_coefficients[18][5] = 1;
    M1a_coefficients[19][5] = 1;
    M1a_coefficients[22][5] = 1;
    M1a_coefficients[23][5] = 1;
    M1a_coefficients[25][5] = 1;
    M1a_coefficients[26][5] = 1;
    M1a_coefficients[27][5] = 1;
    M1a_coefficients[32][5] = 1;
    M1a_coefficients[33][5] = 1;
    M1a_coefficients[34][5] = 1;
    M1a_coefficients[36][5] = 1;
    M1a_coefficients[39][5] = 1;
    M1a_coefficients[40][5] = 1;
    M1a_coefficients[41][5] = 1;
    M1a_coefficients[42][5] = 1;
    M1a_coefficients[43][5] = 1;
    M1a_coefficients[44][5] = 1;
    M1a_coefficients[45][5] = 1;
    M1a_coefficients[46][5] = 1;
    M1a_coefficients[47][5] = 1;
    M1a_coefficients[48][5] = 1;
    M1a_coefficients[49][5] = 1;
    M1a_coefficients[54][5] = 1;
    M1a_coefficients[55][5] = 1;
    M1a_coefficients[57][5] = 1;
    M1a_coefficients[58][5] = 1;
    M1a_coefficients[59][5] = 1;
    M1a_coefficients[60][5] = 1;
    M1a_coefficients[61][5] = 1;
    M1a_coefficients[62][5] = 1;
    // G
    M1a_coefficients[53][6] = 1;
    //
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> M1a(M1a_coefficients);
    //
    Polynomial<fig1a_MaxDeg> PHI1a_coefficients[fig1a_dim];
    for (unsigned int i=0; i < fig1a_dim; i++) {
        for (unsigned int j=0; j <= fig1a_MaxDeg; j++) {
            PHI1a_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> PHI1a(
        PHI1a_coefficients
    );
    CenteredLift(&PHI1a, fig1a_p);
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> E1a = ENCRYPT(
        H1a, M1a, PHI1a, fig1a_p, fig1a_q
    );
    // Decryption
    Hypercomplex<Polynomial<fig1a_MaxDeg>, fig1a_dim> D1a = DECRYPT(
        F1a, E1a, fig1a_p, fig1a_q
    );
    CenteredLift(&M1a, fig1a_p);
    assert( D1a == M1a );
    //
    // PUBLICATION: QR code
    //
    const unsigned int fig1b_dim = 32;
    const unsigned int fig1b_MaxDeg = 28; // N = 29
    const int64_t fig1b_p = 3;
    const int64_t fig1b_q = 1723;
    // Public Key
    Polynomial<fig1b_MaxDeg> F1b_coefficients[fig1b_dim];
    F1b_coefficients[1][0] = 1;
    F1b_coefficients[1][1] = 1;
    F1b_coefficients[1][2] = 1;
    F1b_coefficients[1][8] = 1;
    F1b_coefficients[1][13] = 1;
    F1b_coefficients[1][16] = 1;
    F1b_coefficients[1][22] = 1;
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> F1b(F1b_coefficients);
    CenteredLift(&F1b, fig1b_p);
    Polynomial<fig1b_MaxDeg> G1b_coefficients[fig1b_dim];
    for (unsigned int i=0; i < fig1b_dim; i++) {
        for (unsigned int j=0; j <= fig1b_MaxDeg; j++) {
            G1b_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> G1b(G1b_coefficients);
    CenteredLift(&G1b, fig1b_p);
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> H1b = PUBLICKEY(
        F1b, G1b, fig1b_q
    );
    // Encryption
    int temp_M_arr[32][29] = {
        {1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1},
        {1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1},
        {1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1},
        {1, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1},
        {1, 0, 1, 1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1},
        {1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1},
        {1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
        {0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1},
        {0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1},
        {1, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0},
        {0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 0},
        {1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1},
        {1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1},
        {0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0},
        {1, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1},
        {1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1},
        {1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 0, 1, 1, 1, 1, 0, 0, 0},
        {1, 0, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1},
        {1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0},
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1},
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 1, 1},
        {1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0},
        {1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 0, 0, 0, 0, 0},
        {1, 0, 1, 1, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
        {1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0},
        {1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0},
        {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
        {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
        {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2},
    };
    Polynomial<fig1b_MaxDeg> M1b_coefficients[fig1b_dim];
    for (unsigned int i=0; i < fig1b_dim; i++) {
        for (unsigned int j=0; j <= fig1b_MaxDeg; j++) {
            M1b_coefficients[i][j] = temp_M_arr[i][j];
        }
    }
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> M1b(M1b_coefficients);
    //
    Polynomial<fig1b_MaxDeg> PHI1b_coefficients[fig1b_dim];
    for (unsigned int i=0; i < fig1b_dim; i++) {
        for (unsigned int j=0; j <= fig1b_MaxDeg; j++) {
            PHI1b_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> PHI1b(
        PHI1b_coefficients
    );
    CenteredLift(&PHI1b, fig1b_p);
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> E1b = ENCRYPT(
        H1b, M1b, PHI1b, fig1b_p, fig1b_q
    );
    std::cout << "E[QR]:" << std::endl;
    std::cout << E1b << std::endl;
    // Decryption
    Hypercomplex<Polynomial<fig1b_MaxDeg>, fig1b_dim> D1b = DECRYPT(
        F1b, E1b, fig1b_p, fig1b_q
    );
    CenteredLift(&M1b, fig1b_p);
    assert( D1b == M1b );
    //
    // PUBLICATION: MEME
    //
    const unsigned int fig1c_dim = 128;
    const unsigned int fig1c_MaxDeg = 126; // N = 127
    const int64_t fig1c_p = 17;
    const int64_t fig1c_q = 16777213;
    // Public Key
    Polynomial<fig1c_MaxDeg> F1c_coefficients[fig1c_dim];
    F1c_coefficients[1][0] = 1;
    F1c_coefficients[1][1] = 1;
    F1c_coefficients[1][2] = 1;
    F1c_coefficients[1][8] = 1;
    F1c_coefficients[1][13] = 1;
    F1c_coefficients[1][16] = 1;
    F1c_coefficients[1][22] = 1;
    F1c_coefficients[1][24] = 1;
    F1c_coefficients[1][29] = 1;
    F1c_coefficients[1][31] = 1;
    F1c_coefficients[1][76] = 1;
    F1c_coefficients[1][87] = 1;
    F1c_coefficients[1][88] = 1;
    F1c_coefficients[1][94] = 1;
    F1c_coefficients[1][111] = 1;
    F1c_coefficients[1][122] = 1;
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> F1c(F1c_coefficients);
    CenteredLift(&F1c, fig1c_p);
    Polynomial<fig1c_MaxDeg> G1c_coefficients[fig1c_dim];
    for (unsigned int i=0; i < fig1c_dim; i++) {
        for (unsigned int j=0; j <= fig1c_MaxDeg; j++) {
            G1c_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> G1c(G1c_coefficients);
    CenteredLift(&G1c, fig1c_p);
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> H1c = PUBLICKEY(
        F1c, G1c, fig1c_q
    );
    // Encryption
    int64_t nyanarr[fig1c_dim][fig1c_MaxDeg+1];
    Polynomial<fig1c_MaxDeg> M1c_coefficients[fig1c_dim];
    std::ifstream inputfile("nyan.txt");    
    for (unsigned int i=0; i < fig1c_dim; i++) {
        for (unsigned int j=0; j <= fig1c_MaxDeg; j++) {
            inputfile >> M1c_coefficients[i][j];
        }
    }
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> M1c(M1c_coefficients);
    //
    Polynomial<fig1c_MaxDeg> PHI1c_coefficients[fig1c_dim];
    for (unsigned int i=0; i < fig1c_dim; i++) {
        for (unsigned int j=0; j <= fig1c_MaxDeg; j++) {
            PHI1c_coefficients[i][j] = rand_r(&seedzero) % 3;
        }
    }
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> PHI1c(
        PHI1c_coefficients
    );
    CenteredLift(&PHI1c, fig1c_p);
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> E1c = ENCRYPT(
        H1c, M1c, PHI1c, fig1c_p, fig1c_q
    );
    std::cout << "E[NYAN]:" << std::endl;
    std::cout << E1c << std::endl;
    // Decryption
    Hypercomplex<Polynomial<fig1c_MaxDeg>, fig1c_dim> D1c = DECRYPT(
        F1c, E1c, fig1c_p, fig1c_q
    );
    CenteredLift(&M1c, fig1c_p);
    assert( D1c == M1c );
    //
    return 0;
}
