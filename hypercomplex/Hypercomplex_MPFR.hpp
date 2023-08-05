// Copyright 2020 <Maciek Bak>
/*! \file */
/*
###############################################################################
#
#   Hypercomplex header-only library.
#   Explicit template specialisation & function overloading for mpfr_t type
#
#   AUTHOR: Maciek_Bak
#   AFFILIATION: Department_of_Mathematics_City_University_of_London
#   CONTACT: wsciekly.maciek@gmail.com
#   CREATED: 05-08-2023
#   LICENSE: Apache 2.0
#
###############################################################################
*/

// mark that we included this library
#ifndef HYPERCOMPLEX_HYPERCOMPLEX_MPFR_HPP_
#define HYPERCOMPLEX_HYPERCOMPLEX_MPFR_HPP_

#include "./Hypercomplex.hpp"

static uint64_t MPFR_global_precision;

/** \brief Getter for the global precision of the MPFR variables
  * \return precision in bits
  */
uint64_t get_mpfr_precision() {
    return MPFR_global_precision;
}

/** \brief Setter for the global precision of the MPFR variables
  * \param [in] n positive integer (precision in bits)
  */
void set_mpfr_precision(uint64_t n) {
    MPFR_global_precision = n;
}

/** \brief Wrapper for MPFR memory cleanup
  */
void clear_mpfr_memory() {
    mpfr_free_cache();
    assert(!mpfr_mp_memory_cleanup());
}

/** Partial specialisation of the main class for high precision
  */
template <const uint64_t dim>
class Hypercomplex<mpfr_t, dim> {
 private:
    mpfr_t* arr = new mpfr_t[dim]; // NOLINT
    static inline uint64_t** baseprodabs;
    static inline bool** baseprodpos;

 public:
    /** \brief Basis multiplication table initialiser
    */
    static void init() {
        int64_t** M = new int64_t*[dim];
        for (uint64_t i = 0; i < dim; i++) M[i] = new int64_t[dim];
        M[0][0] = 1;
        uint64_t n = 1;
        while (n != dim) {
            for (uint64_t i=0; i < n; i++) {
                for (uint64_t j=0; j < n; j++) {
                    M[i][n+j] = M[j][i] > 0 ? M[j][i] + n : M[j][i] - n;
                    M[i+n][j] = M[i][j] > 0 ? M[i][j] + n : M[i][j] - n;
                    M[i+n][j] = M[i+n][j] * (j ? -1 : 1);
                    M[i+n][j+n] = -M[j][i] * (j ? -1 : 1);
                }
            }
            n *= 2;
        }
        baseprodabs = new uint64_t*[dim];
        baseprodpos = new bool*[dim];
        for (uint64_t i = 0; i < dim; i++) {
            baseprodabs[i] = new uint64_t[dim];
            baseprodpos[i] = new bool[dim];
        }
        for (uint64_t i=0; i < dim; i++) {
            for (uint64_t j=0; j < dim; j++) {
                baseprodabs[i][j] = std::abs(M[i][j]) - 1;
                baseprodpos[i][j] = (0 < M[i][j]);
            }
        }
        for (uint64_t i = 0; i < dim; i++) {
            delete[] M[i];
        }
        delete[] M;
    }

    /** \brief Cleanup function: free all memory
    */
    static void clear() {
        for (uint64_t i = 0; i < dim; i++) {
            delete[] baseprodabs[i];
            delete[] baseprodpos[i];
        }
        delete[] baseprodabs;
        delete[] baseprodpos;
    }

    /** \brief Optimised multiplication function
    * \param [in] H1 LHS operand
    * \param [in] H2 RHS operand
    * \return new class instance
    */
    static Hypercomplex MUL(
        const Hypercomplex &H1,
        const Hypercomplex &H2
    ) {
        mpfr_t prod;
        mpfr_init2(prod, MPFR_global_precision);
        mpfr_t temparr[dim]; // NOLINT
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(temparr[i], MPFR_global_precision);
        for (uint64_t i=0; i < dim; i++) mpfr_set_zero(temparr[i], 0);
        for (uint64_t i=0; i < dim; i++) {
            for (uint64_t j=0; j < dim; j++) {
                mpfr_mul(prod, H1[i], H2[j], MPFR_RNDN);
                if (Hypercomplex::baseprodpos[i][j]) {
                    mpfr_add(
                        temparr[Hypercomplex::baseprodabs[i][j]],
                        temparr[Hypercomplex::baseprodabs[i][j]],
                        prod,
                        MPFR_RNDN);
                } else {
                    mpfr_sub(
                        temparr[Hypercomplex::baseprodabs[i][j]],
                        temparr[Hypercomplex::baseprodabs[i][j]],
                        prod,
                        MPFR_RNDN);
                }
            }
        }
        Hypercomplex H(temparr);
        for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
        mpfr_clear(prod);
        return H;
    }

    /** \brief This is the main constructor
      * \param [in] ARR array of MPFR numbers
      *
      * Template parameters are:
      * * dimensionality of the algebra
      */
    explicit Hypercomplex(const mpfr_t* ARR) {
        if (dim == 0) {
            delete[] arr;
            throw std::invalid_argument("invalid dimension");
        }
        if ((dim & (dim - 1)) != 0) {
            delete[] arr;
            throw std::invalid_argument("invalid dimension");
        }
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(arr[i], MPFR_global_precision);
        for (uint64_t i=0; i < dim; i++)
            mpfr_set(arr[i], ARR[i], MPFR_RNDN);
    }

    /** \brief This is the copy constructor
      * \param [in] H existing class instance
      *
      * Template parameters are:
      * * dimensionality of the algebra
      */
    Hypercomplex(const Hypercomplex &H) {
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(arr[i], MPFR_global_precision);
        for (uint64_t i=0; i < dim; i++)
            mpfr_set(arr[i], H[i], MPFR_RNDN);
    }

    Hypercomplex() = delete;

    ~Hypercomplex() {
        for (uint64_t i=0; i < dim; i++) mpfr_clear(arr[i]);
        delete[] arr;
    }

    /** \brief Dimensionality getter
      * \return algebraic dimension of the underlying object
      */
    uint64_t _() const { return dim; }

    /** \brief Calculate Euclidean norm of a number
      * \param [in,out] norm MPFR variable for the calculated norm
      * \return exit status
      *
      * Note that the interface of this method is different
      * than for the fully template class.
      * Following the MPFR logic: function takes as an
      * argument a variable to store the output in.
      */
    int32_t norm(mpfr_t norm) const {
        mpfr_t temp;
        mpfr_init2(temp, MPFR_global_precision);
        mpfr_set_zero(norm, 0);
        for (uint64_t i=0; i < dim; i++) {
            mpfr_mul(temp, arr[i], arr[i], MPFR_RNDN);
            mpfr_add(norm, norm, temp, MPFR_RNDN);
        }
        mpfr_sqrt(norm, norm, MPFR_RNDN);
        mpfr_clear(temp);
        return 0;
    }

    /** \brief Calculate inverse of a given number
      * \return new class instance
      */
    Hypercomplex inv() const {
        mpfr_t zero, norm;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_init2(norm, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        (*this).norm(norm);
        if (mpfr_equal_p(norm, zero)) {
            mpfr_clear(zero);
            mpfr_clear(norm);
            throw std::invalid_argument("division by zero");
        } else {
            mpfr_t temparr[dim]; // NOLINT
            for (uint64_t i=0; i < dim; i++)
                mpfr_init2(temparr[i], MPFR_global_precision);
            mpfr_mul(norm, norm, norm, MPFR_RNDN);
            mpfr_div(temparr[0], arr[0], norm, MPFR_RNDN);
            for (uint64_t i=1; i < dim; i++) {
                mpfr_div(temparr[i], arr[i], norm, MPFR_RNDN);
                mpfr_sub(temparr[i], zero, temparr[i], MPFR_RNDN);
            }
            Hypercomplex<mpfr_t, dim> H(temparr);
            mpfr_clear(zero);
            mpfr_clear(norm);
            for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
            return H;
        }
    }

    /** \brief Cast a number into a higher dimension
      * \return new class instance
      *
      * New dimension is passed as a function template parameter,
      * as the return class is not the same as the caller's class.
      */
    template <const uint64_t newdim>
    Hypercomplex<mpfr_t, newdim> expand() const {
        if (newdim <= dim) throw std::invalid_argument("invalid dimension");
        mpfr_t temparr[newdim]; // NOLINT
        for (uint64_t i=0; i < newdim; i++)
            mpfr_init2(temparr[i], MPFR_global_precision);
        for (uint64_t i=0; i < dim; i++)
            mpfr_set(temparr[i], arr[i], MPFR_RNDN);
        for (uint64_t i=dim; i < newdim; i++) mpfr_set_zero(temparr[i], 0);
        Hypercomplex<mpfr_t, newdim> H(temparr);
        for (uint64_t i=0; i < newdim; i++) mpfr_clear(temparr[i]);
        return H;
    }

    /** \brief Create a complex conjugate
      * \return new class instance
      */
    Hypercomplex operator~ () const {
        mpfr_t zero;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        mpfr_t temparr[dim]; // NOLINT
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(temparr[i], MPFR_global_precision);
        mpfr_set(temparr[0], arr[0], MPFR_RNDN);
        for (uint64_t i=1; i < dim; i++)
            mpfr_sub(temparr[i], zero, arr[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
        mpfr_clear(zero);
        return H;
    }

    /** \brief Create an additive inverse of a given number
      * \return new class instance
      */
    Hypercomplex operator- () const {
        mpfr_t zero;
        mpfr_init2(zero, MPFR_global_precision);
        mpfr_set_zero(zero, 0);
        mpfr_t temparr[dim]; // NOLINT
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(temparr[i], MPFR_global_precision);
        for (uint64_t i=0; i < dim; i++)
            mpfr_sub(temparr[i], zero, arr[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
        mpfr_clear(zero);
        return H;
    }

    /** \brief Assignment operator
      * \param [in] H existing class instance
      * \return Reference to the caller (for chained assignments)
      */
    Hypercomplex& operator= (const Hypercomplex &H) {
        if (this == &H) return *this;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set(arr[i], H[i], MPFR_RNDN);
        return *this;
    }

    /** \brief Access operator (const)
      * \param [in] i index for the element to access
      * \return i-th element of the number
      */
    mpfr_t const & operator[] (const uint64_t i) const {
        assert(0 <= i && i < dim);
        return arr[i];
    }

    /** \brief Access operator (non-const)
      * \param [in] i index for the element to access
      * \return i-th element of the number
      */
    mpfr_t & operator[] (const uint64_t i) {
        assert(0 <= i && i < dim);
        return arr[i];
    }

    /** \brief Addition-Assignment operator
      * \param [in] H existing class instance
      * \return Reference to the caller
      */
    Hypercomplex& operator+= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) + H;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set((*this)[i], result[i], MPFR_RNDN);
        return *this;
    }

    /** \brief Subtraction-Assignment operator
      * \param [in] H existing class instance
      * \return Reference to the caller
      */
    Hypercomplex& operator-= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) - H;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set((*this)[i], result[i], MPFR_RNDN);
        return *this;
    }

    /** \brief Multiplication-Assignment operator
      * \param [in] H existing class instance
      * \return Reference to the caller
      */
    Hypercomplex& operator*= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) * H;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set((*this)[i], result[i], MPFR_RNDN);
        return *this;
    }

    /** \brief Power-Assignment operator
      * \param [in] x positive integer
      * \return Reference to the caller
      */
    Hypercomplex& operator^= (const uint64_t x) {
        Hypercomplex<mpfr_t, dim> result = (*this) ^ x;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set((*this)[i], result[i], MPFR_RNDN);
        return *this;
    }

    /** \brief Division-Assignment operator
      * \param [in] H existing class instance
      * \return Reference to the caller
      */
    Hypercomplex& operator/= (const Hypercomplex &H) {
        Hypercomplex<mpfr_t, dim> result = (*this) / H;
        for (uint64_t i=0; i < dim; i++)
            mpfr_set((*this)[i], result[i], MPFR_RNDN);
        return *this;
    }
};

/** \brief Equality operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return boolean value after the comparison
  */
template <const uint64_t dim>
bool operator==(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    for (uint64_t i=0; i < dim; i++) {
        if (!mpfr_equal_p(H1[i], H2[i])) return false;
    }
    return true;
}

/** \brief Inequality operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return boolean value after the comparison
  */
template <const uint64_t dim>
bool operator!=(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    return !(H1 == H2);
}

/** \brief Addition operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> operator+(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    mpfr_t temparr[dim]; // NOLINT
    for (uint64_t i=0; i < dim; i++)
        mpfr_init2(temparr[i], MPFR_global_precision);
    for (uint64_t i=0; i < dim; i++)
        mpfr_add(temparr[i], H1[i], H2[i], MPFR_RNDN);
    Hypercomplex<mpfr_t, dim> H(temparr);
    for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
    return H;
}

/** \brief Subtraction operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> operator-(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    mpfr_t temparr[dim]; // NOLINT
    for (uint64_t i=0; i < dim; i++)
        mpfr_init2(temparr[i], MPFR_global_precision);
    for (uint64_t i=0; i < dim; i++)
        mpfr_sub(temparr[i], H1[i], H2[i], MPFR_RNDN);
    Hypercomplex<mpfr_t, dim> H(temparr);
    for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
    return H;
}

/** \brief Multiplication operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> operator*(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    // recursion base:
    if constexpr (dim == 1) {
        mpfr_t result;
        mpfr_init2(result, MPFR_global_precision);
        mpfr_mul(result, H1[0], H2[0], MPFR_RNDN);
        mpfr_t temparr[1];
        mpfr_init2(temparr[0], MPFR_global_precision);
        mpfr_set(temparr[0], result, MPFR_RNDN);
        Hypercomplex<mpfr_t, 1> H_(temparr);
        mpfr_clear(result);
        mpfr_clear(temparr[0]);
        return H_;
    // recursion step:
    } else {
        // shared objects:
        const uint64_t halfd = dim / 2;
        mpfr_t temparr[dim]; // NOLINT
        for (uint64_t i=0; i < dim; i++)
            mpfr_init2(temparr[i], MPFR_global_precision);
        // construct helper objects:
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i], H1[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, halfd> H1a(temparr);
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i], H1[i+halfd], MPFR_RNDN);
        Hypercomplex<mpfr_t, halfd> H1b(temparr);
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i], H2[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, halfd> H2a(temparr);
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i], H2[i+halfd], MPFR_RNDN);
        Hypercomplex<mpfr_t, halfd> H2b(temparr);
        // multiply recursively:
        Hypercomplex<mpfr_t, halfd> H1a2a = H1a * H2a;
        Hypercomplex<mpfr_t, halfd> H2b_1b = ~H2b * H1b;
        Hypercomplex<mpfr_t, halfd> H2b1a = H2b * H1a;
        Hypercomplex<mpfr_t, halfd> H1b2a_ = H1b * ~H2a;
        // construct the final object
        Hypercomplex<mpfr_t, halfd> Ha = H1a2a - H2b_1b;
        Hypercomplex<mpfr_t, halfd> Hb = H2b1a + H1b2a_;
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i], Ha[i], MPFR_RNDN);
        for (uint64_t i=0; i < halfd; i++)
            mpfr_set(temparr[i+halfd], Hb[i], MPFR_RNDN);
        Hypercomplex<mpfr_t, dim> H(temparr);
        for (uint64_t i=0; i < dim; i++) mpfr_clear(temparr[i]);
        return H;
    }
}

/** \brief Power operator
  * \param [in] H LHS operand
  * \param [in] x RHS operand
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> operator^(
    const Hypercomplex<mpfr_t, dim> &H,
    const uint64_t x
) {
    if (!(x)) {
        throw std::invalid_argument("zero is not a valid argument");
    } else {
        Hypercomplex<mpfr_t, dim> Hx(H);
        for (uint64_t i=0; i < x-1; i++) Hx *= H;
        return Hx;
    }
}

/** \brief Division operator
  * \param [in] H1 LHS operand
  * \param [in] H2 RHS operand
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> operator/(
    const Hypercomplex<mpfr_t, dim> &H1,
    const Hypercomplex<mpfr_t, dim> &H2
) {
    Hypercomplex<mpfr_t, dim> H = H1 * H2.inv();
    return(H);
}

/** \brief Print operator
  * \param [in,out] os output stream
  * \param [in] H existing class instance
  * \return output stream
  */
template <const uint64_t dim>
std::ostream& operator<<(
    std::ostream &os,
    const Hypercomplex<mpfr_t, dim> &H
) {
    mpfr_exp_t exponent = 0;
    char* outstr;
    for (uint64_t i=0; i < dim - 1; i++) {
        outstr = mpfr_get_str(NULL, &exponent, 10, 0, H[i], MPFR_RNDN);
        os << outstr << "E" << exponent << " ";
        mpfr_free_str(outstr);
    }
    outstr = mpfr_get_str(NULL, &exponent, 10, 0, H[dim - 1], MPFR_RNDN);
    os << outstr << "E" << exponent;
    mpfr_free_str(outstr);
    return os;
}

/** \brief Real part of a hypercomplex number
  * \param [in] H existing class instance
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> Re(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = H;
    for (uint64_t i=1; i < dim; i++) mpfr_set_zero(result[i], 0);
    return result;
}

/** \brief Imaginary part of a hypercomplex number
  * \param [in] H existing class instance
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> Im(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = H;
    mpfr_set_zero(result[0], 0);
    return result;
}

/** \brief Exponentiation operation on a hypercomplex number
  * \param [in] H existing class instance
  * \return new class instance
  */
template <const uint64_t dim>
Hypercomplex<mpfr_t, dim> exp(const Hypercomplex<mpfr_t, dim> &H) {
    Hypercomplex<mpfr_t, dim> result = Im(H);
    mpfr_t zero, norm, expreal;
    mpfr_init2(zero, MPFR_global_precision);
    mpfr_init2(norm, MPFR_global_precision);
    mpfr_init2(expreal, MPFR_global_precision);
    mpfr_set_zero(zero, 0);
    result.norm(norm);
    mpfr_exp(expreal, H[0], MPFR_RNDN);

    if (mpfr_equal_p(norm, zero)) {
        mpfr_set(result[0], expreal, MPFR_RNDN);
        for (uint64_t i=1; i < dim; i++) mpfr_set_zero(result[i], 0);
    } else {
        mpfr_t sinv_v;
        mpfr_init2(sinv_v, MPFR_global_precision);
        mpfr_sin(sinv_v, norm, MPFR_RNDN);
        mpfr_div(sinv_v, sinv_v, norm, MPFR_RNDN);
        for (uint64_t i=0; i < dim; i++) {
            mpfr_mul(result[i], result[i], sinv_v, MPFR_RNDN);
        }
        mpfr_cos(norm, norm, MPFR_RNDN);
        mpfr_add(result[0], result[0], norm, MPFR_RNDN);
        for (uint64_t i=0; i < dim; i++) {
            mpfr_mul(result[i], result[i], expreal, MPFR_RNDN);
        }
        mpfr_clear(sinv_v);
    }
    mpfr_clear(zero);
    mpfr_clear(norm);
    mpfr_clear(expreal);
    return result;
}

#endif  // HYPERCOMPLEX_HYPERCOMPLEX_MPFR_HPP_
