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

TEST_CASE( "Overloading Operators", "[operators]" ) {
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
        REQUIRE( h1[0] == A[0] );
        REQUIRE( h1[1] == -A[1] );
        REQUIRE( h1[2] == -A[2] );
        REQUIRE( h1[3] == -A[3] );
    }

    SECTION( "Access operator" ) {
        REQUIRE( h1[0] == A[0] );
        REQUIRE( h1[1] == A[1] );
        REQUIRE( h1[2] == A[2] );
        REQUIRE( h1[3] == A[3] );
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
        REQUIRE( h1[0] == -A[0] );
        REQUIRE( h1[1] == -A[1] );
        REQUIRE( h1[2] == -A[2] );
        REQUIRE( h1[3] == -A[3] );
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

    // + - += -=

    // * / *= /=

    // ^ ^=
}

int main(int argc, char* const argv[]) {
    return Catch::Session().run(argc, argv);
}
