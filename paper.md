---
title: 'Hypercomplex: abstract & fast header-only C++ template library for lattice-based cryptosystems in high-dimensional algebras'
tags:
- Cryptography
- Algebra
- Arbitrary-precision arithmetic
- C++
authors:
- name: Maciek Bak
  orcid: 0000-0003-1361-7301
date: 6 February 2023
bibliography: paper.bib
---

# Summary

The following work presents a *C++* library which is dedicated to performing arbitrary-precise calculations on hypercomplex numbers from the Cayley-Dickson algebras [@schafer2017introduction]. Basic arithmetical operations as well as a few miscellaneous functions are implemented. 
Its most important feature is the support for encryption/decryption procedures for public-key lattice-based cryptosystems in high-dimensional algebras over truncated polynomial rings.

# Statement of need

This is a highly specialised software aimed mostly for computational mathematicians and computational scientists who operate on high-dimensional numbers, need to carry out arbitrary-precise calculations or whose focus is the study of lattice-based, post-quantum cryptography. The library is well suited for wide range of computationally-challenging projects: from investigating general algebraic properties _per se_ to applied research where hypercomplex framework serves merely as a mean to an end (as in previously mentioned cryptosystems).

# Key features

- As a header-only *C++* template code it's greatest advantage is the combination of speed, generic programming and convenience for the end user. Open Source license together with template specialisation mechanism allows contributors to add-in support for custom objects, define specific functions and extend the scope of the library.
- The most important specialisation, already included in the library itself, is the introduction of operations in hypercomplex algebras over truncated polynomial rings. These allow for many cryptographic applications as described in a dedicated section below. 
- Another template class specialisation introduces the support for arbitrary high precision of calculations via GNU MPFR library [@fousse:inria-00070266], for which the operators have been overloaded such that all the instructions are carried out on specific data structures.
- State of the art technology for software engineering:
  - CI/CD mechanism set up with GitHub Actions: automatic tests for library installation, source code inclusion, compilation and execution,
  - extensive unit testing with Catch2 framework [@catch2] alongside code coverage measurement uploaded to Codecov; current coverage: 100%,
  - source code linting with cpplint [@cpplint] - Google code style enforced,
  - automatic documentation generation and hosting on GitHub Pages: build via Doxygen [@doxygen], publishing via Actions.

# Cryptographic applications


Let $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations

\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}

and refer to \autoref{eq:fourier} from text.

![
  Examples of _Hypercomplex_ applications for cryptography.
  **(a)** Message (M) composed of seven, relatively-big
  secret numbers is encoded in binary in a [64x7] matrix.
  This is later encrypted (E) and decrypted (D) with a public-key
  cryptosystem, allowing to transfer the numbers in a secure manner.
  **(b)** Graphical representation of an encrypted/decrypted QR code,
  encoded in a [32x29] matrix (padded image).
  **(c)** Encrypted/Decrypted [128x127] meme; credits: [www.nyan.cat](www.nyan.cat).
](img/Fig1.png)

# State of the field

When it comes to a general hypercomplex framework the well-known _boost C++_ libraries deserve the most notable mention here [@boost]. Unfortunately their scope is limitted as they only provide quaterions and octonions classes (however as an upside - all the operations are well optimised). Moreover, these libraries do not support operations on MPFR types natively. It may also be worth to mention the existence of smaller repositories like: [@quaternions] or [@cd], but, unlike our work, they often lack proper test suites, code coverage reports, documentation and are also significantly restricted in functionality which is a major drawback.

However, (most importantly) to our best knowledge there is currently no high-quality open-source library which natively supports cryptosystems based on truncated polynomial rings.
Previous research described distinct versions of NTRU [c]: 4-dimensional QTRU [c], 8-dimensional OTRU [c]; some proposed 16-dimenisional STRU, which correctness has not yet been verified [c].
Despite these efforts no generalization has been provided yet.
Our work is a first to: present that these procedures are vaild in arbitrary-high-dimensional Cayley-Dickson algebras (provided a careful choice of parameters of the system)
and to provide reproducible examples of a successful encryption/decryption procedures.

# Acknowledgments

We would like to express our wholehearted gratitidue towards: the members of
a facebook group _>implying we can discuss mathematics_, who aided us
with clarifications and suggestions related to the topic of research
as well as a _Cryptography Stack Exchange_ user: _DanielS_, who helped us
analyse and understand specifics of lattice-based cryptosystems.

# References

