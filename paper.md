---
title: 'Precise calculations on hypercomplex numbers in C++'
tags:
- C++
- Algebra
- Arbitrary-precision arithmetic
authors:
- name: Maciej Bak
  orcid: 0000-0003-1361-7301
  affiliation: "1, 2"
affiliations:
 - name: Biozentrum, University of Basel
   index: 1
 - name: Swiss Institute of Bioinformatics
   index: 2
date: 25 August 2020
bibliography: paper.bib
---

# Summary

The following repository contains a C++ library which is dedicated to performing arbitrary-precise computations on hypercomplex numbers from the Cayley-Dickson algebras. Basic arithmetical operations as well as a few miscellaneous functions are implemented. It aims to aid other developers in computational research.

# Statement of need 

This is a reference [@Tadelis:2012]. This is another one [@axelrodproject]. And another one [@Macnamara:2020].

# Key features

- template mechanism for built-in types, open for specialisation with custom classes, umbrella to gather the rest
- Header-only C++ library, easy to use
- arbitrary high precision of calculations (mpfr)
- high-quality software engineering

The following library aims to deliver a simple method to construct hypercomplex numbers from any of the Cayley-Dickson algebras and later perform calculations in a arbitrary-precise arithmetic. It is dedicated mostly to computational mathematicians and computational scientists who use high-dimensional numbers and/or need to carry out highly-precise calculations. As a header-only C++ template code it's greatest advantage is the combination of speed, generic programming and convenience for the end user.

# State of the field

- https://www.boost.org/doc/libs/1_34_1/doc/html/boost_math/
- https://ece.uwaterloo.ca/~dwharder/C++/CQOST/Sedenion/
- https://github.com/ferd36/quaternions

# References

