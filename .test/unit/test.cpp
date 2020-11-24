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
#include "hypercomplex/Hypercomplex.h"
#include <stdexcept>
#include <iostream>

TEST_CASE( "Class Structure", "[unit]" ) {
    //
    SECTION( "Main constructor" ) {
        unsigned int dim = 4;
        float A[] = {1.0, 2.0, 0.0, -1.0};
        float invalidA[] = {1.0, 2.0, 0.0};
        Hypercomplex h1 = Hypercomplex(dim, A);
        REQUIRE_THROWS_AS(Hypercomplex(3, invalidA), std::invalid_argument);

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
            float target3 = 0.0;
            Approx target4 = Approx(0.166).epsilon(0.01);
            Hypercomplex invh1 = h1.inv();
            REQUIRE( invh1[0] == target1 );
            REQUIRE( invh1[1] == target2 );
            REQUIRE( invh1[2] == target3 );
            REQUIRE( invh1[3] == target4 );
            float A0[] = {0.0,0.0};
            REQUIRE_THROWS_AS(Hypercomplex(2, A0).inv(), std::invalid_argument);
        }

        SECTION( "Expansion" ) {
            Hypercomplex hexpanded = h1.expand(8);
            REQUIRE( hexpanded[0] == h1[0] );
            REQUIRE( hexpanded[1] == h1[1] );
            REQUIRE( hexpanded[2] == h1[2] );
            REQUIRE( hexpanded[3] == h1[3] );
            REQUIRE( hexpanded[4] == 0.0 );
            REQUIRE( hexpanded[5] == 0.0 );
            REQUIRE( hexpanded[6] == 0.0 );
            REQUIRE( hexpanded[7] == 0.0 );
            REQUIRE_THROWS_AS(h1.expand(4), std::invalid_argument);
        }

        SECTION( "Real part" ) {
            Hypercomplex real_h1 = Re(h1);
            REQUIRE( real_h1[0] == h1[0] );
            REQUIRE( real_h1[1] == 0 );
            REQUIRE( real_h1[2] == 0 );
            REQUIRE( real_h1[3] == 0 );
        }

        SECTION( "Imaginary part" ) {
            Hypercomplex imaginary_h1 = Im(h1);
            REQUIRE( imaginary_h1[0] == 0 );
            REQUIRE( imaginary_h1[1] == h1[1] );
            REQUIRE( imaginary_h1[2] == h1[2] );
            REQUIRE( imaginary_h1[3] == h1[3] );
        }

        SECTION( "Hypercomplex exponentiation" ) {
            Approx target1 = Approx(-1.678).epsilon(0.01);
            Approx target2 = Approx(1.913).epsilon(0.01);
            float target3 = 0.0;
            Approx target4 = Approx(-0.956).epsilon(0.01);
            Hypercomplex exp_h1 = exp(h1);
            REQUIRE( exp_h1[0] == target1 );
            REQUIRE( exp_h1[1] == target2 );
            REQUIRE( exp_h1[2] == target3 );
            REQUIRE( exp_h1[3] == target4 );
        }
    }

    SECTION( "Main constructor: exception" ) {
        float A1[] = {10.10};
        float A0[] = {};
        REQUIRE_NOTHROW(Hypercomplex(1, A1));
        REQUIRE_THROWS_AS(Hypercomplex(0, A0), std::invalid_argument);
    }

    SECTION( "Copy constructor" ) {
        unsigned int dim = 4;
        float A[] = {1.0, 2.0, 0.0, -1.0};
        Hypercomplex h1 = Hypercomplex(dim, A);
        Hypercomplex h2 = Hypercomplex(h1);
        Hypercomplex h3 = h2;
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
        unsigned int dim = 4;
        float A[] = {1.0, 2.0, 0.0, -1.0};
        // dynamic memory allocation for memory leak test:
        Hypercomplex* h = new Hypercomplex(dim, A);
        delete h;
        REQUIRE( true == true );
    }
}

TEST_CASE( "Overloading Operators", "[unit]" ) {
    //
    unsigned int dim2 = 2;
    unsigned int dim4 = 4;
    float A[] = {1.0, 2.0, 0.0, -1.0};
    float B[] = {-0.5, 1.0, 0.0, 6.0};
    float C[] = {10.0, -10.0};

    Hypercomplex h1 = Hypercomplex(dim4, A);
    Hypercomplex h2 = Hypercomplex(dim4, B);
    Hypercomplex h3 = Hypercomplex(dim2, C);

    SECTION( "Conjugate operator" ) {
        Hypercomplex h1_ = ~h1;
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
        result = h1 == h3;
        REQUIRE( result == false );
        result = h1 == h1;
        REQUIRE( result == true );
    }

    SECTION( "Inequality operator" ) {
        bool result;
        result = h1 != h2;
        REQUIRE( result == true );
        result = h1 != h3;
        REQUIRE( result == true );
        result = h1 != h1;
        REQUIRE( result == false );
    }

    SECTION( "Negation operator" ) {
        Hypercomplex h1_ = -h1;
        REQUIRE( &h1 != &h1_ );
        REQUIRE( h1_[0] == -A[0] );
        REQUIRE( h1_[1] == -A[1] );
        REQUIRE( h1_[2] == -A[2] );
        REQUIRE( h1_[3] == -A[3] );
        unsigned int dim = (-h1)._();
        REQUIRE( dim == dim4 );
    }

    SECTION( "Assignment operator" ) {
        float a[] = {-3.0, 5.0, 2.0, 1.0};
        float b[] = {9.0, 0.0, -4.0, 1.0};
        float c[] = {5.0, 8.0, 0.0, -8.0};
        Hypercomplex ha = Hypercomplex(dim4, a);
        Hypercomplex hb = Hypercomplex(dim4, b);
        Hypercomplex hc = Hypercomplex(dim4, c);
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
        Hypercomplex h = h1 + h2;
        REQUIRE( h[0] == 0.5 );
        REQUIRE( h[1] == 3.0 );
        REQUIRE( h[2] == 0.0 );
        REQUIRE( h[3] == 5.0 );
        REQUIRE_THROWS_AS(h1 + h3, std::invalid_argument);
    }

    SECTION( "Subtraction operator" ) {
        Hypercomplex h = h1 - h2;
        REQUIRE( h[0] == 1.5 );
        REQUIRE( h[1] == 1.0 );
        REQUIRE( h[2] == 0.0 );
        REQUIRE( h[3] == -7.0 );
        REQUIRE_THROWS_AS(h1 - h3, std::invalid_argument);
    }

    SECTION( "Addition-Assignment operator" ) {
        h1 += h2;
        REQUIRE( h1[0] == 0.5 );
        REQUIRE( h1[1] == 3.0 );
        REQUIRE( h1[2] == 0.0 );
        REQUIRE( h1[3] == 5.0 );
        REQUIRE_THROWS_AS(h1 += h3, std::invalid_argument);
    }

    SECTION( "Subtraction-Assignment operator" ) {
        h1 -= h2;
        REQUIRE( h1[0] == 1.5 );
        REQUIRE( h1[1] == 1.0 );
        REQUIRE( h1[2] == 0.0 );
        REQUIRE( h1[3] == -7.0 );
        REQUIRE_THROWS_AS(h1 -= h3, std::invalid_argument);
    }

    SECTION( "Multiplication operator" ) {
        Hypercomplex h = h1 * h2;
        REQUIRE( h[0] == 3.5 );
        REQUIRE( h[1] == 0.0 );
        REQUIRE( h[2] == -13.0 );
        REQUIRE( h[3] == 6.5 );
        REQUIRE_THROWS_AS(h1 * h3, std::invalid_argument);
    }

    SECTION( "Multiplication-Assignment operator" ) {
        h1 *= h2;
        REQUIRE( h1[0] == 3.5 );
        REQUIRE( h1[1] == 0.0 );
        REQUIRE( h1[2] == -13.0 );
        REQUIRE( h1[3] == 6.5 );
        REQUIRE_THROWS_AS(h1 *= h3, std::invalid_argument);
    }

    SECTION( "Power operator" ) {
        REQUIRE_THROWS_AS(h1 ^ 0, std::invalid_argument);
        REQUIRE_NOTHROW(h1 ^ 1);
        Hypercomplex h = h1 ^ 2;
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
        Hypercomplex h = h1 / h2;
        REQUIRE( h[0] == target1 );
        REQUIRE( h[1] == target2 );
        REQUIRE( h[2] == target3 );
        REQUIRE( h[3] == target4 );
        REQUIRE_THROWS_AS(h1 / h3, std::invalid_argument);
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
        REQUIRE_THROWS_AS(h1 /= h3, std::invalid_argument);
    }

    SECTION( "Output stream operator" ) {
        REQUIRE_NOTHROW(std::cout << h1 << std::endl);
    }
}

TEST_CASE( "Special", "[usecase]" ) {
    //
    SECTION( "Multiplication optimization" ) {
        float A[] = {1.51, -1.13, 2.28, -10.77, -2.63, -9.11, 0.01, 4.02};
        float B[] = {-7.32, -0.70, 0.91, 99.32, 8.09, -9.33, 0.84, -5.32};
        float C[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        Hypercomplex h1 = Hypercomplex(8, A);
        Hypercomplex h2 = Hypercomplex(8, B);
        Hypercomplex result = Hypercomplex(8, C);
        for (unsigned int i=0; i < 10000; i++) result = h1 * h2;
        REQUIRE( true == true );
    }

    SECTION( "Const objects" ) {
        const unsigned int dim = 4;
        const unsigned int newdim = 8;
        const unsigned int cui = 2;
        const float A[] = {1.0, 2.0, 0.0, -1.0};
        const float B[] = {-0.5, 1.0, 0.0, 6.0};
        const Hypercomplex const_h1 = Hypercomplex(dim, A);
        const Hypercomplex const_h2 = Hypercomplex(dim, B);
        REQUIRE_NOTHROW(const_h1._());
        REQUIRE_NOTHROW(const_h1.norm());
        REQUIRE_NOTHROW(const_h1.inv());
        REQUIRE_NOTHROW(~const_h1);
        REQUIRE_NOTHROW(-const_h1);
        REQUIRE_NOTHROW(const_h1.expand(newdim));
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

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}
