// Copyright 2020 <Maciej Bak>
/*
###############################################################################
#
#   Hypercomplex library: header file
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 22-10-2020
#   LICENSE: MIT
#
###############################################################################
*/

// mark that we included this header
#ifndef HYPERCOMPLEX_HYPERCOMPLEX_H_
#define HYPERCOMPLEX_HYPERCOMPLEX_H_

/*
Main class of the library
*/
class Hypercomplex {
 private:
    unsigned int d;  // dimension
    float* arr;  // elements
 public:
    explicit Hypercomplex(unsigned int d, float* arr);  // main constructor
    Hypercomplex(const Hypercomplex& H);  // copy constructor
    Hypercomplex() = delete;  // forbid default constructor | c++11
    ~Hypercomplex();  // destructor
    float _() { return d; }  // space dimension
    Hypercomplex operator~ ();  // conjugate operator
    Hypercomplex operator- ();  // negation operator
    Hypercomplex& operator= (const Hypercomplex &H);  // assignment operator
    float& operator[] (unsigned int i);  // access operator
    bool operator== (const Hypercomplex& H);  // equality operator
    bool operator!= (const Hypercomplex& H);  // inequality operator
};

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_H_
