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

#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "hypercomplex/Hypercomplex.hpp"
#include <tuple>
#include <stdexcept>
#include <iostream>

template<typename T>
using Hypercomplex0 = Hypercomplex<T, 0>;
template<typename T>
using Hypercomplex1 = Hypercomplex<T, 1>;
template<typename T>
using Hypercomplex2 = Hypercomplex<T, 2>;
template<typename T>
using Hypercomplex3 = Hypercomplex<T, 3>;

using MPFR_Hypercomplex0 = Hypercomplex<mpfr_t, 0>;
using MPFR_Hypercomplex1 = Hypercomplex<mpfr_t, 1>;
using MPFR_Hypercomplex2 = Hypercomplex<mpfr_t, 2>;
using MPFR_Hypercomplex3 = Hypercomplex<mpfr_t, 3>;

using TestTypes = std::tuple<float, double, long double>;

TEMPLATE_LIST_TEST_CASE( "Class Structure", "[unit]", TestTypes ) {
    //
    SECTION( "Main constructor & functions" ) {
        const unsigned int dim = 4;
        TestType A[] = {1.0, 2.0, 0.0, -1.0};
        TestType invalidA[] = {1.0, 2.0, 0.0};
        Hypercomplex<TestType, dim> h1(A);
        REQUIRE_THROWS_AS(
            Hypercomplex3<TestType>(invalidA),
            std::invalid_argument
        );

        SECTION( "Getters" ) {
            REQUIRE( h1._() == dim );
        }

        SECTION( "Norm" ) {
            Approx target = Approx(2.45).epsilon(0.01);
            REQUIRE( h1.norm() == target );
        }

        SECTION( "Inverse" ) {
            Approx target1 = Approx(0.166).epsilon(0.01);
            Approx target2 = Approx(-0.333).epsilon(0.01);
            TestType target3 = 0.0;
            Approx target4 = Approx(0.166).epsilon(0.01);
            Hypercomplex<TestType, dim> invh1 = h1.inv();
            REQUIRE( invh1[0] == target1 );
            REQUIRE( invh1[1] == target2 );
            REQUIRE( invh1[2] == target3 );
            REQUIRE( invh1[3] == target4 );
            TestType A0[] = {0.0,0.0};
            REQUIRE_THROWS_AS(
                Hypercomplex2<TestType>(A0).inv(),
                std::invalid_argument
            );
        }

        SECTION( "Real part" ) {
            Hypercomplex<TestType, dim> real_h1 = Re(h1);
            REQUIRE( real_h1[0] == h1[0] );
            REQUIRE( real_h1[1] == 0.0 );
            REQUIRE( real_h1[2] == 0.0 );
            REQUIRE( real_h1[3] == 0.0 );
        }

        SECTION( "Imaginary part" ) {
            Hypercomplex<TestType, dim> imaginary_h1 = Im(h1);
            REQUIRE( imaginary_h1[0] == 0.0 );
            REQUIRE( imaginary_h1[1] == h1[1] );
            REQUIRE( imaginary_h1[2] == h1[2] );
            REQUIRE( imaginary_h1[3] == h1[3] );
        }

        SECTION( "Hypercomplex exponentiation" ) {
            Approx target1 = Approx(-1.678).epsilon(0.01);
            Approx target2 = Approx(1.913).epsilon(0.01);
            TestType target3 = 0.0;
            Approx target4 = Approx(-0.956).epsilon(0.01);
            Hypercomplex<TestType, dim> exp_h1 = exp(h1);
            REQUIRE( exp_h1[0] == target1 );
            REQUIRE( exp_h1[1] == target2 );
            REQUIRE( exp_h1[2] == target3 );
            REQUIRE( exp_h1[3] == target4 );
            TestType B[] = {5.0, 0.0, 0.0, 0.0};
            Hypercomplex<TestType, dim> h2(B);
            Hypercomplex<TestType, dim> exp_h2 = exp(h2);
            Approx target5 = Approx(148.413).epsilon(0.01);
            REQUIRE( exp_h2[0] == target5 );
            REQUIRE( exp_h2[1] == 0.0 );
            REQUIRE( exp_h2[2] == 0.0 );
            REQUIRE( exp_h2[3] == 0.0 );
        }
    }

    SECTION( "Main constructor: exception" ) {
        TestType A1[] = {10.10};
        TestType A0[] = {};
        REQUIRE_NOTHROW(Hypercomplex1<TestType>(A1));
        REQUIRE_THROWS_AS(
            Hypercomplex0<TestType>(A0),
            std::invalid_argument
        );
    }

    SECTION( "Copy constructor" ) {
        const unsigned int dim = 4;
        TestType A[] = {1.0, 2.0, 0.0, -1.0};
        Hypercomplex<TestType, dim> h1(A);
        Hypercomplex<TestType, dim> h2(h1);
        Hypercomplex<TestType, dim> h3 = h2;
        REQUIRE( &h1 != &h2 );
        REQUIRE( &h2 != &h3 );
        REQUIRE( &h3 != &h1 );
        REQUIRE( h1._() == h2._() );
        REQUIRE( h2._() == h3._() );
        REQUIRE( h3._() == h1._() );
        REQUIRE( h1[0] == h2[0] );
        REQUIRE( h2[0] == h3[0] );
        REQUIRE( h3[0] == h1[0] );
    }

    SECTION( "Destructor" ) {
        const unsigned int dim = 4;
        TestType A[] = {1.0, 2.0, 0.0, -1.0};
        // dynamic memory allocation for memory leak test:
        Hypercomplex<TestType, dim>* h = new Hypercomplex<TestType, dim>(A);
        delete h;
        REQUIRE( true );
    }
}

TEMPLATE_LIST_TEST_CASE( "Overloading Operators", "[unit]", TestTypes ) {
    //
    const unsigned int dim2 = 2;
    const unsigned int dim4 = 4;
    TestType A[] = {1.0, 2.0, 0.0, -1.0};
    TestType B[] = {-0.5, 1.0, 0.0, 6.0};
    TestType C[] = {10.0, -10.0};

    Hypercomplex<TestType, dim4> h1(A);
    Hypercomplex<TestType, dim4> h2(B);
    Hypercomplex<TestType, dim2> h3(C);

    SECTION( "Conjugate operator" ) {
        Hypercomplex<TestType, dim4> h1_ = ~h1;
        REQUIRE( &h1 != &h1_ );
        REQUIRE( h1_[0] == A[0] );
        REQUIRE( h1_[1] == -A[1] );
        REQUIRE( h1_[2] == -A[2] );
        REQUIRE( h1_[3] == -A[3] );
        unsigned int dim = (~h1)._();
        REQUIRE( dim == dim4 );
    }

    SECTION( "Access operator" ) {
        REQUIRE( h1[0] == A[0] );
        REQUIRE( h1[1] == A[1] );
        REQUIRE( h1[2] == A[2] );
        REQUIRE( h1[3] == A[3] );
        h1[0] = 100;
        REQUIRE( h1[0] == 100 );
    }

    SECTION( "Equality operator" ) {
        bool result;
        result = h1 == h2;
        REQUIRE( result == false );
        result = h1 == h1;
        REQUIRE( result == true );
    }

    SECTION( "Inequality operator" ) {
        bool result;
        result = h1 != h2;
        REQUIRE( result == true );
        result = h1 != h1;
        REQUIRE( result == false );
    }

    SECTION( "Negation operator" ) {
        Hypercomplex<TestType, dim4> h1_ = -h1;
        REQUIRE( &h1 != &h1_ );
        REQUIRE( h1_[0] == -A[0] );
        REQUIRE( h1_[1] == -A[1] );
        REQUIRE( h1_[2] == -A[2] );
        REQUIRE( h1_[3] == -A[3] );
        unsigned int dim = (-h1)._();
        REQUIRE( dim == dim4 );
    }

    SECTION( "Assignment operator" ) {
        TestType a[] = {-3.0, 5.0, 2.0, 1.0};
        TestType b[] = {9.0, 0.0, -4.0, 1.0};
        TestType c[] = {5.0, 8.0, 0.0, -8.0};
        Hypercomplex<TestType, dim4> ha(a);
        Hypercomplex<TestType, dim4> hb(b);
        Hypercomplex<TestType, dim4> hc(c);
        REQUIRE( &h1 != &ha );
        REQUIRE( h1[0] != ha[0] );
        ha = h1;
        REQUIRE( &h1 != &ha );
        REQUIRE( h1[0] == ha[0] );
        // chain assignment:
        hc = hb = ha;
        REQUIRE( &ha != &hb );
        REQUIRE( &hb != &hc );
        REQUIRE( &hc != &ha );    
        REQUIRE( ha[0] == hb[0] );
        REQUIRE( hb[0] == hc[0] );
        REQUIRE( hc[0] == ha[0] );    
        // test self-assignment:
        h1 = h1;
    }

    SECTION( "Addition operator" ) {
        Hypercomplex<TestType, dim4> h = h1 + h2;
        REQUIRE( h[0] == 0.5 );
        REQUIRE( h[1] == 3.0 );
        REQUIRE( h[2] == 0.0 );
        REQUIRE( h[3] == 5.0 );
    }

    SECTION( "Subtraction operator" ) {
        Hypercomplex<TestType, dim4> h = h1 - h2;
        REQUIRE( h[0] == 1.5 );
        REQUIRE( h[1] == 1.0 );
        REQUIRE( h[2] == 0.0 );
        REQUIRE( h[3] == -7.0 );
    }

    SECTION( "Addition-Assignment operator" ) {
        h1 += h2;
        REQUIRE( h1[0] == 0.5 );
        REQUIRE( h1[1] == 3.0 );
        REQUIRE( h1[2] == 0.0 );
        REQUIRE( h1[3] == 5.0 );
    }

    SECTION( "Subtraction-Assignment operator" ) {
        h1 -= h2;
        REQUIRE( h1[0] == 1.5 );
        REQUIRE( h1[1] == 1.0 );
        REQUIRE( h1[2] == 0.0 );
        REQUIRE( h1[3] == -7.0 );
    }

    SECTION( "Multiplication operator" ) {
        Hypercomplex<TestType, dim4> h = h1 * h2;
        REQUIRE( h[0] == 3.5 );
        REQUIRE( h[1] == 0.0 );
        REQUIRE( h[2] == -13.0 );
        REQUIRE( h[3] == 6.5 );
    }

    SECTION( "Multiplication-Assignment operator" ) {
        h1 *= h2;
        REQUIRE( h1[0] == 3.5 );
        REQUIRE( h1[1] == 0.0 );
        REQUIRE( h1[2] == -13.0 );
        REQUIRE( h1[3] == 6.5 );
    }

    SECTION( "Power operator" ) {
        REQUIRE_THROWS_AS(h1 ^ 0, std::invalid_argument);
        REQUIRE_NOTHROW(h1 ^ 1);
        Hypercomplex<TestType, dim4> h = h1 ^ 2;
        REQUIRE( h[0] == -4.0 );
        REQUIRE( h[1] == 4.0 );
        REQUIRE( h[2] == 0.0 );
        REQUIRE( h[3] == -2.0 );
        // test implicit type conversion
        short int si = 2;
        unsigned short int usi = 2;
        int i = 2;
        unsigned int ui = 2;
        long int li = 2;
        unsigned long int uli = 2;
        long long int lli = 2;
        unsigned long long int ulli = 2;
        REQUIRE_NOTHROW(h1 ^ si);
        REQUIRE_NOTHROW(h1 ^ usi);
        REQUIRE_NOTHROW(h1 ^ i);
        REQUIRE_NOTHROW(h1 ^ ui);
        REQUIRE_NOTHROW(h1 ^ li);
        REQUIRE_NOTHROW(h1 ^ uli);
        REQUIRE_NOTHROW(h1 ^ lli);
        REQUIRE_NOTHROW(h1 ^ ulli);
    }

    SECTION( "Power-Assignment operator" ) {
        REQUIRE_THROWS_AS(h1 ^= 0, std::invalid_argument);
        REQUIRE_NOTHROW(h1 ^= 1);
        h1 ^= 2;
        REQUIRE( h1[0] == -4.0 );
        REQUIRE( h1[1] == 4.0 );
        REQUIRE( h1[2] == 0.0 );
        REQUIRE( h1[3] == -2.0 );
        // test implicit type conversion
        short int si = 2;
        unsigned short int usi = 2;
        int i = 2;
        unsigned int ui = 2;
        long int li = 2;
        unsigned long int uli = 2;
        long long int lli = 2;
        unsigned long long int ulli = 2;
        REQUIRE_NOTHROW(h1 ^= si);
        REQUIRE_NOTHROW(h1 ^= usi);
        REQUIRE_NOTHROW(h1 ^= i);
        REQUIRE_NOTHROW(h1 ^= ui);
        REQUIRE_NOTHROW(h1 ^= li);
        REQUIRE_NOTHROW(h1 ^= uli);
        REQUIRE_NOTHROW(h1 ^= lli);
        REQUIRE_NOTHROW(h1 ^= ulli);
    }

    SECTION( "Division operator" ) {
        Approx target1 = Approx(-0.121).epsilon(0.01);
        Approx target2 = Approx(-0.054).epsilon(0.01);
        Approx target3 = Approx(0.350).epsilon(0.01);
        Approx target4 = Approx(-0.148).epsilon(0.01);
        Hypercomplex<TestType, dim4> h = h1 / h2;
        REQUIRE( h[0] == target1 );
        REQUIRE( h[1] == target2 );
        REQUIRE( h[2] == target3 );
        REQUIRE( h[3] == target4 );
        TestType D[] = {0.0, 0.0, 0.0, 0.0};
        Hypercomplex<TestType, dim4> h4(D);
        REQUIRE_THROWS_AS(h1 / h4, std::invalid_argument);
    }

    SECTION( "Division-Assignment operator" ) {
        Approx target1 = Approx(-0.121).epsilon(0.01);
        Approx target2 = Approx(-0.054).epsilon(0.01);
        Approx target3 = Approx(0.350).epsilon(0.01);
        Approx target4 = Approx(-0.148).epsilon(0.01);
        h1 /= h2;
        REQUIRE( h1[0] == target1 );
        REQUIRE( h1[1] == target2 );
        REQUIRE( h1[2] == target3 );
        REQUIRE( h1[3] == target4 );
        TestType D[] = {0.0, 0.0, 0.0, 0.0};
        Hypercomplex<TestType, dim4> h4(D);
        REQUIRE_THROWS_AS(h1 /= h4, std::invalid_argument);
    }

    SECTION( "Output stream operator" ) {
        REQUIRE_NOTHROW(std::cout << h1 << std::endl);
    }
}

TEMPLATE_LIST_TEST_CASE( "Special", "[usecase]", TestTypes ) {
    //
    SECTION( "Multiplication optimization" ) {
        TestType A[] = {1.51, -1.13, 2.28, -10.77, -2.63, -9.11, 0.01, 4.02};
        TestType B[] = {-7.32, -0.70, 0.91, 99.32, 8.09, -9.33, 0.84, -5.32};
        TestType C[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Hypercomplex<TestType, 8> h1(A);
        Hypercomplex<TestType, 8> h2(B);
        Hypercomplex<TestType, 8> result(C);
        for (unsigned int i=0; i < 10000; i++) result = h1 * h2;
        REQUIRE( true == true );
    }

    SECTION( "Const objects" ) {
        const unsigned int dim = 4;
        const unsigned int cui = 2;
        const TestType A[] = {1.0, 2.0, 0.0, -1.0};
        const TestType B[] = {-0.5, 1.0, 0.0, 6.0};
        const Hypercomplex<TestType, dim> const_h1(A);
        const Hypercomplex<TestType, dim> const_h2(B);
        REQUIRE_NOTHROW(const_h1._());
        REQUIRE_NOTHROW(const_h1.norm());
        REQUIRE_NOTHROW(const_h1.inv());
        REQUIRE_NOTHROW(~const_h1);
        REQUIRE_NOTHROW(-const_h1);
        REQUIRE_NOTHROW(const_h1[0]);
        REQUIRE_NOTHROW(const_h1 == const_h2);
        REQUIRE_NOTHROW(const_h1 != const_h2);
        REQUIRE_NOTHROW(const_h1 + const_h2);
        REQUIRE_NOTHROW(const_h1 - const_h2);
        REQUIRE_NOTHROW(const_h1 * const_h2);
        REQUIRE_NOTHROW(const_h1 / const_h2);
        REQUIRE_NOTHROW(const_h1 ^ cui);
        REQUIRE_NOTHROW(std::cout << const_h1 << std::endl);
        REQUIRE_NOTHROW(Re(const_h1));
        REQUIRE_NOTHROW(Im(const_h1));
        REQUIRE_NOTHROW(exp(const_h1));
    }
}

TEST_CASE( "Expansion", "[unit]" ) {
    // expand method is a template member function of a template class
    // as such it cannot be tested within TEMPLATE_LIST_TEST_CASE
    // framework of Catch2
    double A[] = {1.0, 2.0, 0.0, -1.0};
    Hypercomplex<double, 4> h1(A);
    Hypercomplex<double, 8> hexpanded = h1.expand<8>();
    REQUIRE( hexpanded[0] == h1[0] );
    REQUIRE( hexpanded[1] == h1[1] );
    REQUIRE( hexpanded[2] == h1[2] );
    REQUIRE( hexpanded[3] == h1[3] );
    REQUIRE( hexpanded[4] == 0.0 );
    REQUIRE( hexpanded[5] == 0.0 );
    REQUIRE( hexpanded[6] == 0.0 );
    REQUIRE( hexpanded[7] == 0.0 );
    REQUIRE_THROWS_AS(
        h1.expand<4>(),
        std::invalid_argument
    );
    const Hypercomplex<double, 4> const_h1(A);
    REQUIRE_NOTHROW(const_h1.expand<8>());
    // MPFR:
    set_mpfr_precision(200);
    std::cout << "Precision: | " << get_mpfr_precision() << std::endl;
    mpfr_t mpfrA[4];
    mpfr_init2(mpfrA[0], MPFR_global_precision);
    mpfr_init2(mpfrA[1], MPFR_global_precision);
    mpfr_init2(mpfrA[2], MPFR_global_precision);
    mpfr_init2(mpfrA[3], MPFR_global_precision);
    mpfr_set_d(mpfrA[0], 1.0, MPFR_RNDN);
    mpfr_set_d(mpfrA[1], 2.0, MPFR_RNDN);
    mpfr_set_d(mpfrA[2], 0.0, MPFR_RNDN);
    mpfr_set_d(mpfrA[3], -1.0, MPFR_RNDN);
    Hypercomplex<mpfr_t, 4> mpfrh1(mpfrA);
    Hypercomplex<mpfr_t, 8> mpfrhexpanded = mpfrh1.expand<8>();
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[0], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[1], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[2], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[3], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[4], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[5], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[6], MPFR_RNDN);
    std::cout << std::endl;
    mpfr_out_str(stdout, 10, 0, mpfrhexpanded[7], MPFR_RNDN);
    std::cout << std::endl;
    REQUIRE_THROWS_AS(
        mpfrh1.expand<4>(),
        std::invalid_argument
    );
    const Hypercomplex<mpfr_t, 4> const_mpfrh1(mpfrA);
    REQUIRE_NOTHROW(const_mpfrh1.expand<8>());
    mpfr_clear(mpfrA[0]);
    mpfr_clear(mpfrA[1]);
    mpfr_clear(mpfrA[2]);
    mpfr_clear(mpfrA[3]);
    clear_mpfr_memory();
}

TEST_CASE( "MPFR lib test", "[unit]" ) {
    //
    SECTION( "Main constructor & functions" ) {
        set_mpfr_precision(200);
        std::cout << "Precision: | " << get_mpfr_precision() << std::endl;
        mpfr_t A[4];
        mpfr_init2(A[0], MPFR_global_precision);
        mpfr_init2(A[1], MPFR_global_precision);
        mpfr_init2(A[2], MPFR_global_precision);
        mpfr_init2(A[3], MPFR_global_precision);
        mpfr_set_d(A[0], 1.0, MPFR_RNDN);
        mpfr_set_d(A[1], 2.0, MPFR_RNDN);
        mpfr_set_d(A[2], 0.0, MPFR_RNDN);
        mpfr_set_d(A[3], -1.0, MPFR_RNDN);
        Hypercomplex<mpfr_t, 4> h1(A);
        REQUIRE_THROWS_AS(
            MPFR_Hypercomplex3(A),
            std::invalid_argument
        );

        SECTION( "Getters" ) {
            REQUIRE( h1._() == 4 );
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
        }

        SECTION( "Norm" ) {
            mpfr_t norm;
            mpfr_init2(norm, MPFR_global_precision);
            std::cout << 2.45 << std::endl;
            h1.norm(norm);
            mpfr_out_str(stdout, 10, 0, norm, MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(norm);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Inverse" ) {
            mpfr_t target;
            mpfr_init2(target, MPFR_global_precision);
            mpfr_set_d(target, 0.166, MPFR_RNDN);
            Hypercomplex<mpfr_t, 4> invh1 = h1.inv();
            std::cout << 0.166 << std::endl;
            mpfr_out_str(stdout, 10, 0, invh1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_t A0[2];
            mpfr_init2(A0[0], MPFR_global_precision);
            mpfr_init2(A0[1], MPFR_global_precision);
            mpfr_set_zero(A0[0], 0);
            mpfr_set_zero(A0[1], 0);
            REQUIRE_THROWS_AS(
                MPFR_Hypercomplex2(A0).inv(),
                std::invalid_argument
            );
            mpfr_clear(target);
            mpfr_clear(A0[0]);
            mpfr_clear(A0[1]);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
        }

        SECTION( "Real part" ) {
            Hypercomplex<mpfr_t, 4> real_h1 = Re(h1);
            mpfr_out_str(stdout, 10, 0, real_h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, real_h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, real_h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, real_h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Imaginary part" ) {
            Hypercomplex<mpfr_t, 4> imaginary_h1 = Im(h1);
            mpfr_out_str(stdout, 10, 0, imaginary_h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, imaginary_h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, imaginary_h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, imaginary_h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Hypercomplex exponentiation" ) {
            mpfr_t target;
            mpfr_init2(target, MPFR_global_precision);
            mpfr_set_d(target, -1.678, MPFR_RNDN);
            Hypercomplex<mpfr_t, 4> exp_h1 = exp(h1);
            std::cout << -1.678 << std::endl;
            mpfr_out_str(stdout, 10, 0, exp_h1[0], MPFR_RNDN);
            std::cout << std::endl;
            std::cout << "-----" << std::endl;
            mpfr_t B[4];
            mpfr_init2(B[0], MPFR_global_precision);
            mpfr_init2(B[1], MPFR_global_precision);
            mpfr_init2(B[2], MPFR_global_precision);
            mpfr_init2(B[3], MPFR_global_precision);
            mpfr_set_d(B[0], 5.0, MPFR_RNDN);
            mpfr_set_d(B[1], 0.0, MPFR_RNDN);
            mpfr_set_d(B[2], 0.0, MPFR_RNDN);
            mpfr_set_d(B[3], 0.0, MPFR_RNDN);
            Hypercomplex<mpfr_t, 4> h2(B);
            Hypercomplex<mpfr_t, 4> exp_h2 = exp(h2);
            mpfr_set_d(target, 148.413, MPFR_RNDN);
            std::cout << 148.413 << std::endl;
            mpfr_out_str(stdout, 10, 0, exp_h2[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, exp_h2[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, exp_h2[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, exp_h2[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(target);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }
    }

    SECTION( "Main constructor: exception" ) {
        set_mpfr_precision(200);
        mpfr_t A0[] = {};
        REQUIRE_THROWS_AS(
            MPFR_Hypercomplex0(A0),
            std::invalid_argument
        );
        clear_mpfr_memory();
    }

    SECTION( "Copy constructor" ) {
        set_mpfr_precision(200);
        mpfr_t A[4];
        mpfr_init2(A[0], MPFR_global_precision);
        mpfr_init2(A[1], MPFR_global_precision);
        mpfr_init2(A[2], MPFR_global_precision);
        mpfr_init2(A[3], MPFR_global_precision);
        mpfr_set_d(A[0], 1.0, MPFR_RNDN);
        mpfr_set_d(A[1], 2.0, MPFR_RNDN);
        mpfr_set_d(A[2], 0.0, MPFR_RNDN);
        mpfr_set_d(A[3], -1.0, MPFR_RNDN);
        Hypercomplex<mpfr_t, 4> h1(A);
        Hypercomplex<mpfr_t, 4> h2(h1);
        Hypercomplex<mpfr_t, 4> h3 = h2;
        REQUIRE( &h1 != &h2 );
        REQUIRE( &h2 != &h3 );
        REQUIRE( &h3 != &h1 );
        REQUIRE( h1._() == h2._() );
        REQUIRE( h2._() == h3._() );
        REQUIRE( h3._() == h1._() );
        REQUIRE( !mpfr_cmp(h1[0], h2[0]) );
        REQUIRE( !mpfr_cmp(h2[0], h3[0]) );
        REQUIRE( !mpfr_cmp(h3[0], h1[0]) );
        mpfr_clear(A[0]);
        mpfr_clear(A[1]);
        mpfr_clear(A[2]);
        mpfr_clear(A[3]);
        clear_mpfr_memory();
    }

    SECTION( "Destructor" ) {
        const unsigned int dim = 4;
        set_mpfr_precision(200);
        mpfr_t A[4];
        mpfr_init2(A[0], MPFR_global_precision);
        mpfr_init2(A[1], MPFR_global_precision);
        mpfr_init2(A[2], MPFR_global_precision);
        mpfr_init2(A[3], MPFR_global_precision);
        mpfr_set_d(A[0], 1.0, MPFR_RNDN);
        mpfr_set_d(A[1], 2.0, MPFR_RNDN);
        mpfr_set_d(A[2], 0.0, MPFR_RNDN);
        mpfr_set_d(A[3], -1.0, MPFR_RNDN);
        Hypercomplex<mpfr_t, dim>* h = new Hypercomplex<mpfr_t, dim>(A);
        delete h;
        mpfr_clear(A[0]);
        mpfr_clear(A[1]);
        mpfr_clear(A[2]);
        mpfr_clear(A[3]);
        clear_mpfr_memory();
        REQUIRE( true );
    }
}

TEST_CASE( "MPFR: const objects", "[unit]" ) {
    const unsigned int dim = 4;
    const unsigned int cui = 2;    
    set_mpfr_precision(200);
    mpfr_t norm;
    mpfr_init2(norm, MPFR_global_precision);
    mpfr_t A[4], B[4];
    mpfr_init2(A[0], MPFR_global_precision);
    mpfr_init2(A[1], MPFR_global_precision);
    mpfr_init2(A[2], MPFR_global_precision);
    mpfr_init2(A[3], MPFR_global_precision);
    mpfr_init2(B[0], MPFR_global_precision);
    mpfr_init2(B[1], MPFR_global_precision);
    mpfr_init2(B[2], MPFR_global_precision);
    mpfr_init2(B[3], MPFR_global_precision);
    mpfr_set_d(A[0], 1.0, MPFR_RNDN);
    mpfr_set_d(A[1], 2.0, MPFR_RNDN);
    mpfr_set_d(A[2], 0.0, MPFR_RNDN);
    mpfr_set_d(A[3], -1.0, MPFR_RNDN);
    mpfr_set_d(B[0], -0.5, MPFR_RNDN);
    mpfr_set_d(B[1], 1.0, MPFR_RNDN);
    mpfr_set_d(B[2], 0.0, MPFR_RNDN);
    mpfr_set_d(B[3], 6.0, MPFR_RNDN);
    const Hypercomplex<mpfr_t, dim> const_h1(A);
    const Hypercomplex<mpfr_t, dim> const_h2(B);
    REQUIRE_NOTHROW(const_h1._());
    REQUIRE_NOTHROW(const_h1.norm(norm));
    REQUIRE_NOTHROW(const_h1.inv());
    REQUIRE_NOTHROW(~const_h1);
    REQUIRE_NOTHROW(-const_h1);
    REQUIRE_NOTHROW(const_h1[0]);
    REQUIRE_NOTHROW(const_h1 == const_h2);
    REQUIRE_NOTHROW(const_h1 != const_h2);
    REQUIRE_NOTHROW(const_h1 + const_h2);
    REQUIRE_NOTHROW(const_h1 - const_h2);
    REQUIRE_NOTHROW(const_h1 * const_h2);
    REQUIRE_NOTHROW(const_h1 / const_h2);
    REQUIRE_NOTHROW(const_h1 ^ cui);
    REQUIRE_NOTHROW(std::cout << const_h1 << std::endl);
    REQUIRE_NOTHROW(Re(const_h1));
    REQUIRE_NOTHROW(Im(const_h1));
    REQUIRE_NOTHROW(exp(const_h1));
    mpfr_clear(norm);
    mpfr_clear(A[0]);
    mpfr_clear(A[1]);
    mpfr_clear(A[2]);
    mpfr_clear(A[3]);
    mpfr_clear(B[0]);
    mpfr_clear(B[1]);
    mpfr_clear(B[2]);
    mpfr_clear(B[3]);
    clear_mpfr_memory();
}

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}
