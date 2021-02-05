---
title: 'Precise calculations on hypercomplex numbers in C++'
tags:
- Algebra
- Arbitrary-precision arithmetic
- C++
authors:
- name: Maciej Bak
  orcid: 0000-0003-1361-7301
  affiliation: "1, 2"
affiliations:
 - name: Biozentrum, University of Basel
   index: 1
 - name: Swiss Institute of Bioinformatics
   index: 2
date: 4 February 2021
bibliography: paper.bib
---

# Summary

The following repository contains a *C++* library which is dedicated to performing arbitrary-precise calculations on hypercomplex numbers from the Cayley-Dickson algebras. Basic arithmetical operations as well as a few miscellaneous functions are implemented. Its focus is to aid other developers in computational research.

# Statement of need

This is a highly specialised software aimed mostly for computational mathematicians and computational scientists who use high-dimensional numbers and/or need to carry out highly-precise calculations.
(generality, fills in gap in the open source community?)

# Key features

- As a header-only *C++* template code it's greatest advantage is the combination of speed, generic programming and convenience for the end user. Open Source license together with template specialisation mechanism allows contributors to add-in support for custom objects, define specific functions and extend the scope of the library.
- One of such specialisation is already included in the library itself - a support for arbitrary high precision of calculations via GNU MPFR library(c)...
- State of the art technology for software engineering:
  - CI/CD mechanism set up with GitHub Actions: automatic tests for library installation, source code inclusion, compilation and execution,
  - extensive unit testing with Catch2 framework(c) and code coverage measurement uploaded to Codecov(c); current coverage: 100%,
  - source code linting with cpplint(c),
  - automatic documentation generation and hosting on GitHub Pages (generation via Doxygen(c), publishing via Actions).

# State of the field

abc [@boost] [@quaternions] [@cd][@fousse:inria-00070266] [@schafer2017introduction]

Neither of the available repositories provide both 

# References

