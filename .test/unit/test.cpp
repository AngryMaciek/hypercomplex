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
#include <iostream>

TEST_CASE( "Class Structure", "[class]" ) {
    //
    SECTION( "Main constructor" ) {
        unsigned int dim = 4;
        float A[] = {1.0, 2.0, 0.0, -1.0};
        Hypercomplex h1 = Hypercomplex(dim, A);

        SECTION( "Getters" ) {
            REQUIRE( h1._() == dim );
        }
    }

    // default constructor

    // copy constructor

    // BOOST: octonion & university lib - methods?

    SECTION( "Destructor" ) {
        unsigned int dim = 4;
        float A[] = {1.0, 2.0, 0.0, -1.0};
        // dynamic memory allocation for memory leak test:
        Hypercomplex* h = new Hypercomplex(dim, A);
        delete h;
        REQUIRE( true == true );
    }
}

TEST_CASE( "Overloading Operators", "[operators]" ) {
    //
    unsigned int dim = 4;
    float A[] = {1.0, 2.0, 0.0, -1.0};
    float B[] = {-0.5, 1.0, 0.0, 6.0};

    Hypercomplex h1 = Hypercomplex(dim, A);
    Hypercomplex h2 = Hypercomplex(dim, B);

    SECTION( "Conjugate operator" ) {
        Hypercomplex h1_ = ~h1;
        REQUIRE( h1_.arr[0] == A[0] );
        REQUIRE( h1_.arr[1] == -A[1] );
        REQUIRE( h1_.arr[2] == -A[2] );
        REQUIRE( h1_.arr[3] == -A[3] );
    }

    SECTION( "Access operator" ) {
        REQUIRE( h1[0] == A[0] );
        REQUIRE( h1[1] == A[1] );
        REQUIRE( h1[2] == A[2] );
        REQUIRE( h1[3] == A[3] );
    }

    SECTION( "Equality operator" ) {
        bool result = h1 == h2;
        REQUIRE( result == false );
    }

    SECTION( "Inequality operator" ) {
        bool result = h1 != h2;
        REQUIRE( result == true );
    }

    SECTION( "Negation operator" ) {
        Hypercomplex h1_ = -h1;
        REQUIRE( h1_.arr[0] == -A[0] );
        REQUIRE( h1_.arr[1] == -A[1] );
        REQUIRE( h1_.arr[2] == -A[2] );
        REQUIRE( h1_.arr[3] == -A[3] );
    }

    SECTION( "Assignment operator" ) {
        float a[] = {-3.0, 5.0, 2.0, 1.0};
        float b[] = {9.0, 0.0, -4.0, 1.0};
        float c[] = {5.0, 8.0, 0.0, -8.0};
        Hypercomplex ha = Hypercomplex(dim, a);
        Hypercomplex hb = Hypercomplex(dim, b);
        Hypercomplex hc = Hypercomplex(dim, c);
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

    // + - += -=

    // * / *= /=

    // ^ ^=
}

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}
