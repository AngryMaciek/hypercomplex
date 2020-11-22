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

#include <iosfwd>

/*
Main class of the library
*/
class Hypercomplex {

 private:
    unsigned int d;
    float *arr;

 public:
    explicit
        Hypercomplex(const unsigned int arg_d, const float* arg_arr);
    Hypercomplex(const Hypercomplex& H);
    Hypercomplex() = delete;  // forbid default constructor | c++11
    ~Hypercomplex();
    float _() const { return d; }
    float norm() const;
    Hypercomplex inv() const;
    Hypercomplex expand(const unsigned int arg_d) const;
    Hypercomplex Re() const;
    Hypercomplex Im() const;
    Hypercomplex operator~ () const;
    Hypercomplex operator- () const;
    Hypercomplex& operator= (const Hypercomplex &H);
    float& operator[] (const unsigned int i) const;
    Hypercomplex& operator+= (const Hypercomplex &H);
    Hypercomplex& operator-= (const Hypercomplex &H);
    Hypercomplex& operator*= (const Hypercomplex &H);
    Hypercomplex& operator^= (const unsigned int x);
    Hypercomplex& operator/= (const Hypercomplex &H);
};

/*
Operators
*/
bool operator== (const Hypercomplex &H1, const Hypercomplex &H2);
bool operator!= (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator+ (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator- (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator* (const Hypercomplex &H1, const Hypercomplex &H2);
Hypercomplex operator^ (const Hypercomplex &H, const unsigned int x);
Hypercomplex operator/ (const Hypercomplex &H1, const Hypercomplex &H2);
std::ostream& operator<< (std::ostream &os, const Hypercomplex &H);

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_H_
