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
    unsigned int d;
 public:
    float* arr;  // move to private once [] is implemented
    explicit Hypercomplex(unsigned int d, float* arr);  // constructor
    ~Hypercomplex();  // destructor
    float _() { return d; }  // space dimension
    Hypercomplex operator~ ();  // conjugate operator
    Hypercomplex operator- ();
    Hypercomplex& operator= (const Hypercomplex &H)
    float& operator[] (unsigned int i);
    bool operator== (const Hypercomplex& H);
    bool operator!= (const Hypercomplex& H);
};

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_H_
