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
#include "hypercomplex/Polynomial.hpp"
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

using Polynomial4_Hypercomplex0 = Hypercomplex<Polynomial<4>, 0>;
using Polynomial4_Hypercomplex1 = Hypercomplex<Polynomial<4>, 1>;
using Polynomial4_Hypercomplex2 = Hypercomplex<Polynomial<4>, 2>;
using Polynomial4_Hypercomplex3 = Hypercomplex<Polynomial<4>, 3>;

using TestTypes = std::tuple<float, double, long double>;

TEMPLATE_LIST_TEST_CASE(
    "Hypercomplex: Class Structure",
    "[unit]",
    TestTypes
) {
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

TEMPLATE_LIST_TEST_CASE(
    "Hypercomplex: Overloading Operators",
    "[unit]",
    TestTypes
) {
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

TEMPLATE_LIST_TEST_CASE( "Hypercomplex: Special", "[usecase]", TestTypes ) {
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

TEST_CASE( "Hypercomplex: Expansion", "[unit]" ) {
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
    // Polynomial:
    int64_t array1[] = {0, 0, 2, 0, 2};
    int64_t array2[] = {1, 1, 2, 0, 2};
    int64_t array3[] = {3, 0, 2, 1, 2};
    int64_t array4[] = {0, 1, 1, 0, 3};
    Polynomial<4> polynomial1(array1);
    Polynomial<4> polynomial2(array2);
    Polynomial<4> polynomial3(array3);
    Polynomial<4> polynomial4(array4);
    Polynomial<4> coefficients[] = {
        polynomial1, polynomial2, polynomial3, polynomial4
    };
    Polynomial<4> zero;
    Hypercomplex<Polynomial<4>, 4> h(coefficients);
    Hypercomplex<Polynomial<4>, 8> phexpanded = h.expand<8>();
    REQUIRE( phexpanded[0] == h[0] );
    REQUIRE( phexpanded[1] == h[1] );
    REQUIRE( phexpanded[2] == h[2] );
    REQUIRE( phexpanded[3] == h[3] );
    REQUIRE( phexpanded[4] == zero );
    REQUIRE( phexpanded[5] == zero );
    REQUIRE( phexpanded[6] == zero );
    REQUIRE( phexpanded[7] == zero );
    REQUIRE_THROWS_AS(
        h.expand<4>(),
        std::invalid_argument
    );
    const Hypercomplex<Polynomial<4>, 4> const_h(coefficients);
    REQUIRE_NOTHROW(const_h.expand<8>());
}

TEST_CASE( "Hypercomplex: MPFR lib test", "[unit]" ) {
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

    SECTION( "Overloading Operators" ) {
        set_mpfr_precision(200);
        const unsigned int dim2 = 2;
        const unsigned int dim4 = 4;
        mpfr_t A[4], B[4], C[2];
        mpfr_init2(A[0], MPFR_global_precision);
        mpfr_init2(A[1], MPFR_global_precision);
        mpfr_init2(A[2], MPFR_global_precision);
        mpfr_init2(A[3], MPFR_global_precision);
        mpfr_init2(B[0], MPFR_global_precision);
        mpfr_init2(B[1], MPFR_global_precision);
        mpfr_init2(B[2], MPFR_global_precision);
        mpfr_init2(B[3], MPFR_global_precision);
        mpfr_init2(C[0], MPFR_global_precision);
        mpfr_init2(C[1], MPFR_global_precision);
        mpfr_set_d(A[0], 1.0, MPFR_RNDN);
        mpfr_set_d(A[1], 2.0, MPFR_RNDN);
        mpfr_set_d(A[2], 0.0, MPFR_RNDN);
        mpfr_set_d(A[3], -1.0, MPFR_RNDN);
        mpfr_set_d(B[0], -0.5, MPFR_RNDN);
        mpfr_set_d(B[1], 1.0, MPFR_RNDN);
        mpfr_set_d(B[2], 0.0, MPFR_RNDN);
        mpfr_set_d(B[3], 6.0, MPFR_RNDN);
        mpfr_set_d(C[0], 10.0, MPFR_RNDN);
        mpfr_set_d(C[1], -10.0, MPFR_RNDN);
        Hypercomplex<mpfr_t, dim4> h1(A);
        Hypercomplex<mpfr_t, dim4> h2(B);
        Hypercomplex<mpfr_t, dim2> h3(C);

        SECTION( "Conjugate operator" ) {
            Hypercomplex<mpfr_t, dim4> h1_ = ~h1;
            REQUIRE( &h1 != &h1_ );
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            std::cout << "-----" << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[3], MPFR_RNDN);
            std::cout << std::endl;
            unsigned int dim = (~h1)._();
            REQUIRE( dim == dim4 );
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Negation operator" ) {
            Hypercomplex<mpfr_t, dim4> h1_ = -h1;
            REQUIRE( &h1 != &h1_ );
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            std::cout << "-----" << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1_[3], MPFR_RNDN);
            std::cout << std::endl;
            unsigned int dim = (-h1)._();
            REQUIRE( dim == dim4 );
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Access operator" ) {
            mpfr_t v;
            mpfr_init2(v, MPFR_global_precision);
            mpfr_set_d(v, 100.0, MPFR_RNDN);
            mpfr_set(h1[0], v, MPFR_RNDN);
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(v);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Equality operator" ) {
            bool result;
            result = h1 == h2;
            REQUIRE( result == false );
            result = h1 == h1;
            REQUIRE( result == true );
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Inequality operator" ) {
            bool result;
            result = h1 != h2;
            REQUIRE( result == true );
            result = h1 != h1;
            REQUIRE( result == false );
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Assignment operator" ) {
            mpfr_t a[4], b[4], c[4];
            mpfr_init2(a[0], MPFR_global_precision);
            mpfr_init2(a[1], MPFR_global_precision);
            mpfr_init2(a[2], MPFR_global_precision);
            mpfr_init2(a[3], MPFR_global_precision);
            mpfr_init2(b[0], MPFR_global_precision);
            mpfr_init2(b[1], MPFR_global_precision);
            mpfr_init2(b[2], MPFR_global_precision);
            mpfr_init2(b[3], MPFR_global_precision);
            mpfr_init2(c[0], MPFR_global_precision);
            mpfr_init2(c[1], MPFR_global_precision);
            mpfr_init2(c[2], MPFR_global_precision);
            mpfr_init2(c[3], MPFR_global_precision);
            mpfr_set_d(a[0], -3.0, MPFR_RNDN);
            mpfr_set_d(a[1], 5.0, MPFR_RNDN);
            mpfr_set_d(a[2], 2.0, MPFR_RNDN);
            mpfr_set_d(a[3], 1.0, MPFR_RNDN);
            mpfr_set_d(b[0], 9.0, MPFR_RNDN);
            mpfr_set_d(b[1], 0.0, MPFR_RNDN);
            mpfr_set_d(b[2], -4.0, MPFR_RNDN);
            mpfr_set_d(b[3], 1.0, MPFR_RNDN);
            mpfr_set_d(c[0], 5.0, MPFR_RNDN);
            mpfr_set_d(c[1], 8.0, MPFR_RNDN);
            mpfr_set_d(c[2], 0.0, MPFR_RNDN);
            mpfr_set_d(c[3], -8.0, MPFR_RNDN);
            Hypercomplex<mpfr_t, dim4> ha(a);
            Hypercomplex<mpfr_t, dim4> hb(b);
            Hypercomplex<mpfr_t, dim4> hc(c);
            REQUIRE( &h1 != &ha );
            REQUIRE( mpfr_cmp(h1[0], ha[0]) );
            ha = h1;
            REQUIRE( &h1 != &ha );
            REQUIRE( !mpfr_cmp(h1[0], ha[0]) );
            hc = hb = ha;
            REQUIRE( &ha != &hb );
            REQUIRE( &hb != &hc );
            REQUIRE( &hc != &ha );
            REQUIRE( !mpfr_cmp(ha[0], hb[0]) );
            REQUIRE( !mpfr_cmp(hb[0], hc[0]) );
            REQUIRE( !mpfr_cmp(hc[0], ha[0]) );
            h1 = h1;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            mpfr_clear(a[0]);
            mpfr_clear(a[1]);
            mpfr_clear(a[2]);
            mpfr_clear(a[3]);
            mpfr_clear(b[0]);
            mpfr_clear(b[1]);
            mpfr_clear(b[2]);
            mpfr_clear(b[3]);
            mpfr_clear(c[0]);
            mpfr_clear(c[1]);
            mpfr_clear(c[2]);
            mpfr_clear(c[3]);
            clear_mpfr_memory();
        }

        SECTION( "Addition operator" ) {
            Hypercomplex<mpfr_t, dim4> h = h1 + h2;
            mpfr_out_str(stdout, 10, 0, h[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Subtraction operator" ) {
            Hypercomplex<mpfr_t, dim4> h = h1 - h2;
            mpfr_out_str(stdout, 10, 0, h[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Addition-Assignment operator" ) {
            h1 += h2;
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Subtraction-Assignment operator" ) {
            h1 -= h2;
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Multiplication operator" ) {
            Hypercomplex<mpfr_t, dim4> h = h1 * h2;
            mpfr_out_str(stdout, 10, 0, h[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Multiplication-Assignment operator" ) {
            h1 *= h2;
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Power operator" ) {
            REQUIRE_THROWS_AS(h1 ^ 0, std::invalid_argument);
            REQUIRE_NOTHROW(h1 ^ 1);
            REQUIRE_NOTHROW(h1 ^ 2);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Power-Assignment operator" ) {
            REQUIRE_THROWS_AS(h1 ^= 0, std::invalid_argument);
            REQUIRE_NOTHROW(h1 ^= 1);
            REQUIRE_NOTHROW(h1 ^= 2);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }

        SECTION( "Division operator" ) {
            Hypercomplex<mpfr_t, dim4> h = h1 / h2;
            mpfr_out_str(stdout, 10, 0, h[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_t D[4];
            mpfr_init2(D[0], MPFR_global_precision);
            mpfr_init2(D[1], MPFR_global_precision);
            mpfr_init2(D[2], MPFR_global_precision);
            mpfr_init2(D[3], MPFR_global_precision);
            mpfr_set_zero(D[0], 0);
            mpfr_set_zero(D[1], 0);
            mpfr_set_zero(D[2], 0);
            mpfr_set_zero(D[3], 0);
            Hypercomplex<mpfr_t, dim4> h4(D);
            REQUIRE_THROWS_AS(h1 / h4, std::invalid_argument);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            mpfr_clear(D[0]);
            mpfr_clear(D[1]);
            mpfr_clear(D[2]);
            mpfr_clear(D[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Division-Assignment operator" ) {
            h1 /= h2;
            mpfr_out_str(stdout, 10, 0, h1[0], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[1], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[2], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_out_str(stdout, 10, 0, h1[3], MPFR_RNDN);
            std::cout << std::endl;
            mpfr_t D[4];
            mpfr_init2(D[0], MPFR_global_precision);
            mpfr_init2(D[1], MPFR_global_precision);
            mpfr_init2(D[2], MPFR_global_precision);
            mpfr_init2(D[3], MPFR_global_precision);
            mpfr_set_zero(D[0], 0);
            mpfr_set_zero(D[1], 0);
            mpfr_set_zero(D[2], 0);
            mpfr_set_zero(D[3], 0);
            Hypercomplex<mpfr_t, dim4> h4(D);
            REQUIRE_THROWS_AS(h1 / h4, std::invalid_argument);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            mpfr_clear(D[0]);
            mpfr_clear(D[1]);
            mpfr_clear(D[2]);
            mpfr_clear(D[3]);
            clear_mpfr_memory();
            REQUIRE( true );
        }

        SECTION( "Output stream operator" ) {
            mpfr_t X[8];
            mpfr_init2(X[0], MPFR_global_precision);
            mpfr_init2(X[1], MPFR_global_precision);
            mpfr_init2(X[2], MPFR_global_precision);
            mpfr_init2(X[3], MPFR_global_precision);
            mpfr_init2(X[4], MPFR_global_precision);
            mpfr_init2(X[5], MPFR_global_precision);
            mpfr_init2(X[6], MPFR_global_precision);
            mpfr_init2(X[7], MPFR_global_precision);
            mpfr_set_d(X[0], 0.0, MPFR_RNDN);
            mpfr_set_d(X[1], 1.0, MPFR_RNDN);
            mpfr_set_d(X[2], -1.0, MPFR_RNDN);
            mpfr_set_d(X[3], 123.456, MPFR_RNDN);
            mpfr_set_d(X[4], -99.9, MPFR_RNDN);
            mpfr_set_d(X[5], 0.0000000001, MPFR_RNDN);
            mpfr_set_d(X[6], -1.0000000001, MPFR_RNDN);
            mpfr_set_d(X[7], 123456789.123456789, MPFR_RNDN);
            Hypercomplex<mpfr_t, 8> hx(X);
            REQUIRE_NOTHROW(std::cout << hx << std::endl);
            mpfr_clear(X[0]);
            mpfr_clear(X[1]);
            mpfr_clear(X[2]);
            mpfr_clear(X[3]);
            mpfr_clear(X[4]);
            mpfr_clear(X[5]);
            mpfr_clear(X[6]);
            mpfr_clear(X[7]);
            mpfr_clear(A[0]);
            mpfr_clear(A[1]);
            mpfr_clear(A[2]);
            mpfr_clear(A[3]);
            mpfr_clear(B[0]);
            mpfr_clear(B[1]);
            mpfr_clear(B[2]);
            mpfr_clear(B[3]);
            mpfr_clear(C[0]);
            mpfr_clear(C[1]);
            clear_mpfr_memory();
        }
    }
}

TEST_CASE( "Hypercomplex: MPFR: const objects", "[unit]" ) {
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

TEST_CASE( "RingInverse", "[unit]" ) {
    //
    SECTION( "[+]" ) {
        REQUIRE( RingInverse(5, 17) == 7 );
        REQUIRE( RingInverse(1, 17) == 1 );
        REQUIRE( RingInverse(10842, 17) == 4 );
        REQUIRE( RingInverse(5, 2) == 1 );
        REQUIRE( RingInverse(1, 2) == 1 );
        REQUIRE( RingInverse(19573442, 3) == 2 );
        REQUIRE( RingInverse(19573444, 3) == 1 );
        REQUIRE( RingInverse(7, 10) == 3 );
    }

    SECTION( "[-]" ) {
        REQUIRE( RingInverse(-5, 17) == 10 );
        REQUIRE( RingInverse(-1, 17) == 16 );
        REQUIRE( RingInverse(-10842, 17) == 13 );
        REQUIRE( RingInverse(-5, 2) == 1 );
        REQUIRE( RingInverse(-1, 2) == 1 );
        REQUIRE( RingInverse(-19573442, 3) == 1 );
        REQUIRE( RingInverse(-19573444, 3) == 2 );
        REQUIRE( RingInverse(-7, 10) == 7 );
    }

    SECTION( "[exception]" ) {
        REQUIRE_THROWS_AS( RingInverse(0, 17), std::invalid_argument );
        REQUIRE_THROWS_AS( RingInverse(0, 2), std::invalid_argument );
        REQUIRE_THROWS_AS( RingInverse(-9, 3), std::invalid_argument );
        REQUIRE_THROWS_AS( RingInverse(0, 10), std::invalid_argument );
        REQUIRE_THROWS_AS( RingInverse(2, 10), std::invalid_argument );
        REQUIRE_THROWS_AS( RingInverse(-6, 10), std::invalid_argument );
    }
}

TEST_CASE( "Polynomial: Class Structure", "[unit]" ) {
    //
    SECTION( "Main constructor" ) {
        const unsigned int deg = 4;
        int64_t coefficients[] = {100, -1, 2, 0, 0};
        Polynomial<deg> P(coefficients);
        REQUIRE( true );
    }

    SECTION( "Copy constructor" ) {
        const unsigned int deg = 4;
        int64_t coefficients[] = {100, -1, 2, 0, 0};
        Polynomial<deg> P1(coefficients);
        Polynomial<deg> P2(P1);
        Polynomial<deg> P3 = P2;
        REQUIRE( &P1 != &P2 );
        REQUIRE( &P2 != &P3 );
        REQUIRE( &P3 != &P1 );
        REQUIRE( P1[4] == P2[4] );
        REQUIRE( P2[4] == P3[4] );
        REQUIRE( P3[4] == P1[4] );
    }

    SECTION( "Default constructor" ) {
        const unsigned int deg = 2;
        Polynomial<deg> P;
        REQUIRE( P[0] == 0 );
        REQUIRE( P[1] == 0 );
        REQUIRE( P[2] == 0 );
    }

    SECTION( "Destructor" ) {
        const unsigned int deg = 4;
        int64_t coefficients[] = {100, -1, 2, 0, 0};
        // dynamic memory allocation for memory leak test:
        Polynomial<deg>* P = new Polynomial<deg>(coefficients);
        delete P;
        REQUIRE( true );
    }
}

TEST_CASE( "Polynomial: Overloading Operators", "[unit]" ) {
    //
    const unsigned int deg = 4;
    int64_t coefficients1[] = {100, -1, 2, 0, 0};
    int64_t coefficients2[] = {-2, 0, 6, 9, -4};
    int64_t coefficients3[] = {-20, 4, 1, -4, 10};
    Polynomial<deg> P1(coefficients1);
    Polynomial<deg> P2(coefficients2);
    Polynomial<deg> P3(coefficients3);

    SECTION( "Access operator" ) {
        REQUIRE( P1[0] == coefficients1[0] );
        REQUIRE( P1[1] == coefficients1[1] );
        REQUIRE( P1[2] == coefficients1[2] );
        REQUIRE( P1[3] == coefficients1[3] );
        REQUIRE( P1[4] == coefficients1[4] );
        P1[0] = 10000;
        REQUIRE( P1[0] == 10000 );
    }

    SECTION( "Equality operator" ) {
        bool result;
        result = P1 == P2;
        REQUIRE( result == false );
        result = P1 == P1;
        REQUIRE( result == true );
    }

    SECTION( "Inequality operator" ) {
        bool result;
        result = P1 != P2;
        REQUIRE( result == true );
        result = P1 != P1;
        REQUIRE( result == false );
    }

    SECTION( "Negation operator" ) {
        Polynomial<deg> P1_ = -P1;
        REQUIRE( &P1 != &P1_ );
        REQUIRE( P1_[0] == -coefficients1[0] );
        REQUIRE( P1_[1] == -coefficients1[1] );
        REQUIRE( P1_[2] == -coefficients1[2] );
        REQUIRE( P1_[3] == -coefficients1[3] );
        REQUIRE( P1_[4] == -coefficients1[4] );
    }

    SECTION( "Assignment operator" ) {
        int64_t coefficientsA[] = {1, 2, 0, 1, -4};
        int64_t coefficientsB[] = {2, 20, 200, 2000, 20000};
        int64_t coefficientsC[] = {-5, -23, -43, -662, -2934};
        Polynomial<deg> PA(coefficientsA);
        Polynomial<deg> PB(coefficientsB);
        Polynomial<deg> PC(coefficientsC);
        REQUIRE( &P1 != &PA );
        REQUIRE( P1[0] != PA[0] );
        PA = P1;
        REQUIRE( &P1 != &PA );
        REQUIRE( P1[0] == PA[0] );
        // chain assignment:
        PC = PB = PA;
        REQUIRE( &PA != &PB );
        REQUIRE( &PB != &PC );
        REQUIRE( &PC != &PA );
        REQUIRE( PA[0] == PB[0] );
        REQUIRE( PB[0] == PC[0] );
        REQUIRE( PC[0] == PA[0] );
        // test self-assignment:
        P1 = P1;
    }

    SECTION( "Addition operator" ) {
        Polynomial<deg> P = P1 + P2;
        REQUIRE( P[0] == 98 );
        REQUIRE( P[1] == -1 );
        REQUIRE( P[2] == 8 );
        REQUIRE( P[3] == 9 );
        REQUIRE( P[4] == -4 );
    }

    SECTION( "Subtraction operator" ) {
        Polynomial<deg> P = P1 - P2;
        REQUIRE( P[0] == 102 );
        REQUIRE( P[1] == -1 );
        REQUIRE( P[2] == -4 );
        REQUIRE( P[3] == -9 );
        REQUIRE( P[4] == 4 );
    }

    SECTION( "Multiplication-by-scalar operator" ) {
        Polynomial<deg> P = 3 * P1;
        REQUIRE( P[0] == 300 );
        REQUIRE( P[1] == -3 );
        REQUIRE( P[2] == 6 );
        REQUIRE( P[3] == 0 );
        REQUIRE( P[4] == 0 );
    }

    SECTION( "Convolution-Multiplication operator" ) {
        Polynomial<deg> P = P1 * P2;
        REQUIRE( P[0] == -178 );
        REQUIRE( P[1] == -6 );
        REQUIRE( P[2] == 596 );
        REQUIRE( P[3] == 894 );
        REQUIRE( P[4] == -397 );
    }

    SECTION( "Output stream operator" ) {
        REQUIRE_NOTHROW(std::cout << P1 << std::endl);
    }

    SECTION( "Modulo operator" ) {
        int64_t coefficients_x[] = {-2, 43, 1, 0, -110, 4125375, -258731};
        Polynomial<6> Px(coefficients_x);
        //
        Polynomial P_87 = Px % 87;
        REQUIRE( P_87[0] == 85 );
        REQUIRE( P_87[1] == 43 );
        REQUIRE( P_87[2] == 1 );
        REQUIRE( P_87[3] == 0 );
        REQUIRE( P_87[4] == 64 );
        REQUIRE( P_87[5] == 9 );
        REQUIRE( P_87[6] == 7 );
        //
        Polynomial P_10 = Px % 10;
        REQUIRE( P_10[0] == 8 );
        REQUIRE( P_10[1] == 3 );
        REQUIRE( P_10[2] == 1 );
        REQUIRE( P_10[3] == 0 );
        REQUIRE( P_10[4] == 0 );
        REQUIRE( P_10[5] == 5 );
        REQUIRE( P_10[6] == 9 );
        //
        Polynomial P_2 = Px % 2;
        REQUIRE( P_2[0] == 0 );
        REQUIRE( P_2[1] == 1 );
        REQUIRE( P_2[2] == 1 );
        REQUIRE( P_2[3] == 0 );
        REQUIRE( P_2[4] == 0 );
        REQUIRE( P_2[5] == 1 );
        REQUIRE( P_2[6] == 1 );
    }
}

TEST_CASE( "Polynomial: const objects", "[unit]" ) {
    const int64_t coefficientsA[] = {1, 2, 0, 1, -4};
    const int64_t coefficientsB[] = {2, 20, 200, 2000, 20000};
    const Polynomial<4> const_PA(coefficientsA);
    const Polynomial<4> const_PB(coefficientsB);
    REQUIRE_NOTHROW(-const_PA);
    REQUIRE_NOTHROW(const_PA[0]);
    REQUIRE_NOTHROW(const_PA == const_PB);
    REQUIRE_NOTHROW(const_PA != const_PB);
    REQUIRE_NOTHROW(const_PA + const_PB);
    REQUIRE_NOTHROW(const_PA - const_PB);
    REQUIRE_NOTHROW(5 * const_PA);
    REQUIRE_NOTHROW(const_PA * const_PB);
    REQUIRE_NOTHROW(const_PA % 3);
    REQUIRE_NOTHROW(std::cout << const_PA << std::endl);
}

TEST_CASE( "Polynomial: CenteredLift function", "[unit]" ) {
    //
    int64_t coefficients_1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    Polynomial<10> P1(coefficients_1);
    CenteredLift(P1, 13);
    REQUIRE( P1[0] == 0 );
    REQUIRE( P1[1] == 1 );
    REQUIRE( P1[2] == 2 );
    REQUIRE( P1[3] == 3 );
    REQUIRE( P1[4] == 4 );
    REQUIRE( P1[5] == 5 );
    REQUIRE( P1[6] == 6 );
    REQUIRE( P1[7] == -6 );
    REQUIRE( P1[8] == -5 );
    REQUIRE( P1[9] == -4 );
    //
    int64_t coefficients_2[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    const Polynomial<10> P2(coefficients_2);
    CenteredLift(P2, 12);
    REQUIRE( P2[0] == 0 );
    REQUIRE( P2[1] == 1 );
    REQUIRE( P2[2] == 2 );
    REQUIRE( P2[3] == 3 );
    REQUIRE( P2[4] == 4 );
    REQUIRE( P2[5] == 5 );
    REQUIRE( P2[6] == 6 );
    REQUIRE( P2[7] == -5 );
    REQUIRE( P2[8] == -4 );
    REQUIRE( P2[9] == -3 );
}

TEST_CASE( "Polynomial: RingInverse function", "[unit]" ) {
    //
    SECTION( "Test 0: x^5-1 | mod5" ) {
        int64_t coefficients[] = {0, 3, 0, 0, 2};
        Polynomial<4> P(coefficients);
        REQUIRE_THROWS_AS(
            RingInverse(P, 5),
            std::invalid_argument
        );
    }

    SECTION( "Test 1: x^11-1 | mod3" ) {
        int64_t coefficients[] = {2, 1, 1, 0, 2, 0, 1, 0, 0, 1, 2};
        Polynomial<10> P(coefficients);
        int64_t coefficients_inv[] = {1, 2, 0, 2, 2, 1, 0, 2, 1, 2, 0};
        Polynomial<10> Pinv(coefficients_inv);
        REQUIRE( RingInverse(P, 3) == Pinv );
    }

    SECTION( "Test 2: x^11-1 | mod37" ) {
        int64_t coefficients[] = {36, 1, 1, 0, 36, 0, 1, 0, 0, 1, 36};
        Polynomial<10> P(coefficients);
        int64_t coefficients_inv[] = {13, 0, 8, 34, 36, 14, 9, 5, 33, 12, 22};
        Polynomial<10> Pinv(coefficients_inv);
        REQUIRE( RingInverse(P, 37) == Pinv );
    }

    SECTION( "Test 3: x^5-1 | mod5" ) {
        int64_t coefficients[] = {0, 0, 2, 0, 2};
        const Polynomial<4> P(coefficients);
        int64_t coefficients_inv[] = {1, 4, 4, 4, 1};
        const Polynomial<4> Pinv(coefficients_inv);
        REQUIRE( RingInverse(P, 5) == Pinv );
    }

    SECTION( "Test 4: x^5-1 | mod5" ) {
        int64_t coefficients[] = {1, 0, 0, 2, 0};
        const Polynomial<4> P(coefficients);
        int64_t coefficients_inv[] = {2, 3, 2, 1, 4};
        const Polynomial<4> Pinv(coefficients_inv);
        REQUIRE( RingInverse(P, 5) == Pinv );
    }

    SECTION( "Test 5: x^5-1 | mod5" ) {
        int64_t coefficients[] = {1, 1, 1, 0, 0};
        const Polynomial<4> P(coefficients);
        int64_t coefficients_inv[] = {4, 3, 3, 4, 3};
        const Polynomial<4> Pinv(coefficients_inv);
        REQUIRE( RingInverse(P, 5) == Pinv );
    }

    SECTION( "Test 6: x^5-1 | mod7" ) {
        int64_t coefficients[] = {0, 0, 0, 0, 0};
        Polynomial<4> P(coefficients);
        REQUIRE_THROWS_AS(
            RingInverse(P, 7),
            std::invalid_argument
        );
    }
}

TEST_CASE( "Hypercomplex: Polynomial lib test", "[unit]" ) {
    //
    SECTION( "Main constructor & functions" ) {
        const unsigned int dim = 4;
        const unsigned int MaxDeg = 4;
        int64_t array1[] = {0, 0, 2, 0, 2};
        int64_t array2[] = {1, 1, 2, 0, 2};
        int64_t array3[] = {3, 0, 2, 1, 2};
        int64_t array4[] = {0, 1, 1, 0, 3};
        Polynomial<MaxDeg> polynomial1(array1);
        Polynomial<MaxDeg> polynomial2(array2);
        Polynomial<MaxDeg> polynomial3(array3);
        Polynomial<MaxDeg> polynomial4(array4);
        Polynomial<MaxDeg> coefficients[] = {
            polynomial1, polynomial2, polynomial3, polynomial4
        };
        Polynomial<MaxDeg> invalid[] = {
            polynomial1, polynomial2, polynomial3
        };
        Hypercomplex<Polynomial<MaxDeg>, dim> h1(coefficients);
        REQUIRE_THROWS_AS(
            Polynomial4_Hypercomplex3(invalid),
            std::invalid_argument
        );

        SECTION( "Getters" ) {
            REQUIRE( h1._() == dim );
        }

        SECTION( "Norm" ) {
            int64_t array[] = {24, 33, 22, 33, 29};
            Polynomial<MaxDeg> polynomial(array);
            REQUIRE( h1.norm() == polynomial );
        }

        SECTION( "Real part" ) {
            Polynomial<MaxDeg> zero;
            Hypercomplex<Polynomial<MaxDeg>, dim> real_h1 = Re(h1);
            REQUIRE( real_h1[0] == h1[0] );
            REQUIRE( real_h1[1] == zero );
            REQUIRE( real_h1[2] == zero );
            REQUIRE( real_h1[3] == zero );
        }

        SECTION( "Imaginary part" ) {
            Polynomial<MaxDeg> zero;
            Hypercomplex<Polynomial<MaxDeg>, dim> imaginary_h1 = Im(h1);
            REQUIRE( imaginary_h1[0] == zero );
            REQUIRE( imaginary_h1[1] == h1[1] );
            REQUIRE( imaginary_h1[2] == h1[2] );
            REQUIRE( imaginary_h1[3] == h1[3] );
        }

        SECTION( "CenteredLift function" ) {
            int64_t coefficients_1[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
            int64_t coefficients_2[] = {0, 1, 2, 3, 4, 10, 9, 8, 7, 6, 5};
            Polynomial<10> P1(coefficients_1);
            Polynomial<10> P2(coefficients_2);
            Polynomial<10> coefficients[] = { P1, P2 };
            Hypercomplex<Polynomial<10>, 2> H(coefficients);
            CenteredLift(H, 13);
            int64_t coefficients_1x[] = {0, 1, 2, 3, 4, 5, 6, -6, -5, -4, -3};
            int64_t coefficients_2x[] = {0, 1, 2, 3, 4, -3, -4, -5, -6, 6, 5};
            Polynomial<10> P1x(coefficients_1x);
            Polynomial<10> P2x(coefficients_2x);
            REQUIRE( H[0] == P1x );
            REQUIRE( H[1] == P2x );
        }
    }

    SECTION( "Main constructor: exception" ) {
        const unsigned int MaxDeg = 4;
        int64_t array1[] = {0, 0, 2, 0, 2};
        Polynomial<MaxDeg> polynomial1(array1);
        Polynomial<MaxDeg> coefficients1[] = { polynomial1 };
        Polynomial<MaxDeg> coefficients0[] = {};
        REQUIRE_NOTHROW(Polynomial4_Hypercomplex1(coefficients1));
        REQUIRE_THROWS_AS(
            Polynomial4_Hypercomplex0(coefficients0),
            std::invalid_argument
        );
    }

    SECTION( "Copy constructor" ) {
        const unsigned int dim = 2;
        const unsigned int MaxDeg = 4;
        int64_t array1[] = {0, 0, 2, 0, 2};
        int64_t array2[] = {3, 3, 2, 0, 2};
        Polynomial<MaxDeg> polynomial1(array1);
        Polynomial<MaxDeg> polynomial2(array2);
        Polynomial<MaxDeg> coefficients[] = { polynomial1, polynomial2 };
        Hypercomplex<Polynomial<MaxDeg>, dim> h1(coefficients);
        Hypercomplex<Polynomial<MaxDeg>, dim> h2(h1);
        Hypercomplex<Polynomial<MaxDeg>, dim> h3 = h2;
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
        const unsigned int dim = 2;
        const unsigned int MaxDeg = 4;
        int64_t array1[] = {0, 0, 2, 0, 2};
        int64_t array2[] = {3, 3, 2, 0, 2};
        Polynomial<MaxDeg> polynomial1(array1);
        Polynomial<MaxDeg> polynomial2(array2);
        Polynomial<MaxDeg> coefficients[] = { polynomial1, polynomial2 };
        // dynamic memory allocation for memory leak test:
        Hypercomplex<Polynomial<MaxDeg>, dim>* h = 
            new Hypercomplex<Polynomial<MaxDeg>, dim>(coefficients);
        delete h;
        REQUIRE( true );
    }

    SECTION( "Overloading Operators" ) {
        const unsigned int dim = 4;
        const unsigned int MaxDeg = 4;
        //
        int64_t array1[] = {1, 0, 2, 0, 2};
        int64_t array2[] = {2, 1, 0, 0, 2};
        int64_t array3[] = {3, 0, 2, 2, 2};
        int64_t array4[] = {0, 3, 1, 1, 3};
        int64_t array5[] = {0, 4, 2, 0, 2};
        int64_t array6[] = {1, 4, 2, 0, 2};
        int64_t array7[] = {3, 0, 2, 4, 4};
        int64_t array8[] = {1, 1, 2, 4, 3};
        int64_t array9[] = {1, 2, 3, 0, 2};
        int64_t array10[] = {3, 4, 1, 1, 1};
        int64_t array11[] = {4, 4, 4, 4, 4};
        int64_t array12[] = {4, 1, 4, 1, 3};
        //
        Polynomial<MaxDeg> polynomial1(array1);
        Polynomial<MaxDeg> polynomial2(array2);
        Polynomial<MaxDeg> polynomial3(array3);
        Polynomial<MaxDeg> polynomial4(array4);
        Polynomial<MaxDeg> polynomial5(array5);
        Polynomial<MaxDeg> polynomial6(array6);
        Polynomial<MaxDeg> polynomial7(array7);
        Polynomial<MaxDeg> polynomial8(array8);
        Polynomial<MaxDeg> polynomial9(array9);
        Polynomial<MaxDeg> polynomial10(array10);
        Polynomial<MaxDeg> polynomial11(array11);
        Polynomial<MaxDeg> polynomial12(array12);
        //
        Polynomial<MaxDeg> coefficientsA[] = {
            polynomial1, polynomial2, polynomial3, polynomial4
        };
        Polynomial<MaxDeg> coefficientsB[] = {
            polynomial5, polynomial6, polynomial7, polynomial8
        };
        Polynomial<MaxDeg> coefficientsC[] = {
            polynomial9, polynomial10, polynomial11, polynomial12
        };
        //
        Hypercomplex<Polynomial<MaxDeg>, dim> hA(coefficientsA);
        Hypercomplex<Polynomial<MaxDeg>, dim> hB(coefficientsB);
        Hypercomplex<Polynomial<MaxDeg>, dim> hC(coefficientsC);

        SECTION( "Conjugate operator" ) {
            Hypercomplex<Polynomial<MaxDeg>, dim> hA_ = ~hA;
            REQUIRE( &hA != &hA_ );
            REQUIRE( hA_[0] == coefficientsA[0] );
            REQUIRE( hA_[1] == -coefficientsA[1] );
            REQUIRE( hA_[2] == -coefficientsA[2] );
            REQUIRE( hA_[3] == -coefficientsA[3] );
            unsigned int dimX = (~hA)._();
            REQUIRE( dimX == dim );
        }

        SECTION( "Access operator" ) {
            REQUIRE( hA[0] == coefficientsA[0] );
            REQUIRE( hA[1] == coefficientsA[1] );
            REQUIRE( hA[2] == coefficientsA[2] );
            REQUIRE( hA[3] == coefficientsA[3] );
            int64_t array[] = {1, 0, 0, 0, 1};
            Polynomial<MaxDeg> polynomial(array);
            hA[0] = polynomial;
            REQUIRE( hA[0] == polynomial );
        }

        SECTION( "Equality operator" ) {
            bool result;
            result = hA == hB;
            REQUIRE( result == false );
            result = hA == hA;
            REQUIRE( result == true );
        }

        SECTION( "Inequality operator" ) {
            bool result;
            result = hA != hB;
            REQUIRE( result == true );
            result = hA != hA;
            REQUIRE( result == false );
        }

        SECTION( "Negation operator" ) {
            Hypercomplex<Polynomial<MaxDeg>, dim> hA_ = -hA;
            REQUIRE( &hA != &hA_ );
            REQUIRE( hA_[0] == -coefficientsA[0] );
            REQUIRE( hA_[1] == -coefficientsA[1] );
            REQUIRE( hA_[2] == -coefficientsA[2] );
            REQUIRE( hA_[3] == -coefficientsA[3] );
            unsigned int dimX = (-hA)._();
            REQUIRE( dimX == dim );
        }

        SECTION( "Assignment operator" ) {
            int64_t array_1[] = {1, 1, 1, 0, 1};
            int64_t array_2[] = {2, 1, 0, 0, 1};
            int64_t array_3[] = {3, 0, 1, 1, 2};
            int64_t array_4[] = {0, 0, 0, 0, 3};
            Polynomial<MaxDeg> polynomial_1(array_1);
            Polynomial<MaxDeg> polynomial_2(array_2);
            Polynomial<MaxDeg> polynomial_3(array_3);
            Polynomial<MaxDeg> polynomial_4(array_4);
            Polynomial<MaxDeg> coefficients_A[] = {
                polynomial_1, polynomial_2, polynomial_3, polynomial_4
            };
            Hypercomplex<Polynomial<MaxDeg>, dim> h_A(coefficients_A);
            //
            REQUIRE( &hA != &h_A );
            REQUIRE( hA[0] != h_A[0] );
            h_A = hA;
            REQUIRE( &hA != &h_A );
            REQUIRE( hA[0] == h_A[0] );
            // chain assignment:
            hC = hB = hA;
            REQUIRE( &hA != &hB );
            REQUIRE( &hB != &hC );
            REQUIRE( &hC != &hA );    
            REQUIRE( hA[0] == hB[0] );
            REQUIRE( hB[0] == hC[0] );
            REQUIRE( hC[0] == hA[0] );    
            // test self-assignment:
            hA = hA;
        }

        SECTION( "Output stream operator" ) {
            REQUIRE_NOTHROW(std::cout << hA << std::endl);
        }
    }
}

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}
