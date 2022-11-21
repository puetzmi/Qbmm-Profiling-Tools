/**
 * @file core_inversion.hpp
 * @author Michele Puetz
 * @brief Algorithms to compute tridiagonal Jacobi matrix from moments.
 * @version 0.1
 * @date 2022-09-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef CORE_INVERSION_HPP
#define CORE_INVERSION_HPP

#include "global_defs.hpp"
#include <stdexcept>
#include <cmath>
#include <mkl.h>
#include <factory.hpp>
#include "cholesky.hpp"


class CoreInversionAlgorithm;       // forward declaration

/// Factory for `CoreInversionAlgorithm` objects
using CoreInversionAlgorithmFactory = Factory<CoreInversionAlgorithm, int>;


/**
 * @brief Base class for core inversion algorithms computing the tridiagonal Jacobi matrix from moments.
 * 
 * Derived classes provide functions that take a realizable set of the first @f$ 2N @f$ integer moments 
 * as input and compute the tridiagonal Jacobi matrix
 * @f\[
 * \mathbf{J} = 
 *  \begin{pmatrix}
 *  \alpha_0 & \gamma_1 &  \\
 *  \gamma_1 & \alpha_1 & \ddots  \\
 *   & \ddots & \ddots & \gamma_{N-1} \\
 *   & & \gamma_{N-1} & \alpha_{N-1}
 *  \end{pmatrix},
 * @f\]
 * where @f$ \gamma_k =  \sqrt{\beta_k} @f$ and @f$ \alpha_k @f$ and @f$ \beta_k @f$ are the recurrence coefficients
 * of the orthogonal polynomials with respect to @f$ n(\xi) @f$ (see e.g. @cite Gautschi2004).
 * 
 */
class CoreInversionAlgorithm {

protected:

    int nMoments_;          ///< Number of moments.
    double *work_;      ///< Workspace (size depends on specific algorithm)

public:

    /**
     * @brief Construct a new `CoreInversionAlgorithm` object.
     * 
     * @param nMoments Number of moments.
     * @param nWorkAlloc Number of doubles for which workspace must be allocated.
     */
    CoreInversionAlgorithm(int nMoments, int nWorkAlloc);

    /**
     * @brief Destroy the `CoreInversionAlgorithm` object.
     * 
     */
    virtual ~CoreInversionAlgorithm();

    /**
     * @brief Compute tridiagonal Jacobi matrix from moments, see class description of `CoreInversionAlgorithm`.
     * 
     * @param[in] moments Set of moments.
     * @param[out] alpha Main diagonal elements of the Jacobi matrix (allocated memory for at least N doubles).
     * @param[out] gamma First upper/lower diagonal elements of the Jacobi matrix (allocated memory for at least N-1 doubles).
     * @return int Error flag.
     */
    virtual int compute(double *moments, double *alpha, double *gamma) = 0;

    /**
     * @brief Alternative method to compute tridiagonal Jacobi matrix from moments (alternative method), see `compute` method.
     * 
     */
    int operator()(double *moments, double *alpha, double *gamma);


    /**
     * @brief Create new `CoreInversionAlgorithm` object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `CoreInversionAlgorithmFactory`.
     * @param nMoments Number of moments.
     * @return std::unique_ptr<CoreInversionAlgorithm> Unique pointer to new `CoreInversionAlgorithm` object.
     */
    static std::unique_ptr<CoreInversionAlgorithm> makeUnique(const std::string &typeName, int nMoments);
    

    /**
     * @brief Create new `CoreInversionAlgorithm` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `CoreInversionAlgorithmFactory`.
     * @param nMoments Number of moments.
     * @return std::shared_ptr<CoreInversionAlgorithm> Shared pointer to new `CoreInversionAlgorithm` object.
     */
    static std::shared_ptr<CoreInversionAlgorithm> makeShared(const std::string &typeName, int nMoments);
};


/**
 * @brief Class for the computation of the tridiagonal Jacobi matrix from moments using the 
 * long quotient-modified difference algorithm @cite Sack1971 as described in @cite John2012.
 * 
 */
class LqmdAlgorithm : public CoreInversionAlgorithm {

public:

    /**
     * @brief Construct a new `LqmdAlgorithm` object.
     * 
     * @param nMoments Number of moments.
     */
    LqmdAlgorithm(int nMoments);

    /**
     * @brief Destroy the `LqmdAlgorithm` object.
     * 
     */
    virtual ~LqmdAlgorithm();


    /**
     * @brief Compute the tridiagonal Jacobi matrix from moments using the 
     * long quotient-modified difference algorithm @cite Sack1971 as described in @cite John2012.
     * 
     * @param[in] moments Set of moments.
     * @param[out] alpha Main diagonal elements of the Jacobi matrix (allocated memory for at least N doubles).
     * @param[out] gamma First upper/lower diagonal elements of the Jacobi matrix (allocated memory for at least N-1 doubles).
     * @return int Error flag.
     */
    virtual int compute(double *moments, double *alpha, double *gamma);

    /**
     * @brief Alternative method to compute tridiagonal Jacobi matrix from moments (alternative method), see `compute` method.
     * 
     */
    int operator()(double *moments, double *alpha, double *gamma);

private:

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, LqmdAlgorithm);

};


/**
 * @brief Class for the computation of the tridiagonal Jacobi matrix from moments using the 
 * algorithm proposed by Golub & Welsch @cite Golub1969.
 * 
 * @tparam CholeskyDecompositionType The type of Cholesky decomposition used in the algorithm.
 */
template <class CholeskyDecompositionType>
class GolubWelschAlgorithm : public CoreInversionAlgorithm {

public:

    /**
     * @brief Construct a new `GolubWelschAlgorithm` object.
     * 
     * @param nMoments Number of moments.
     */
    GolubWelschAlgorithm(int nMoments);

    /**
     * @brief Destroy the `GolubWelschAlgorithm` object.
     * 
     */
    virtual ~GolubWelschAlgorithm();


    /**
     * @brief Compute the tridiagonal Jacobi matrix from moments using the algorithm due to Golub and Welsch @cite Golub1969.
     * 
     * @param[in] moments Set of moments.
     * @param[out] alpha Main diagonal elements of the Jacobi matrix (allocated memory for at least N doubles).
     * @param[out] gamma First upper/lower diagonal elements of the Jacobi matrix (allocated memory for at least N-1 doubles).
     * @return int Error flag.
     */
    virtual int compute(double *mom, double *alpha, double *gamma);

    /**
     * @brief Alternative method to compute tridiagonal Jacobi matrix from moments (alternative method), see `compute` method.
     * 
     */
    int operator()(double *moments, double *alpha, double *gamma);

private:

    std::unique_ptr<CholeskyDecomposition> choleskyDecomposition_;    ///< Pointer to CholeskyDecomposition object

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, GolubWelschAlgorithm);

};


class GolubWelschAlgorithmLapackPotrf 
    : 
    public GolubWelschAlgorithm<CholeskyDecompositionLapackPotrf>
{

public:

    /**
     * @brief Construct a new `GolubWelschAlgorithmLapackPotrf` object
     * 
     * @param nMoments Number of moments.
     */
    GolubWelschAlgorithmLapackPotrf(int nMoments);


    /**
     * @brief Destroy the `GolubWelschAlgorithmLapackPotrf` object
     * 
     */
    ~GolubWelschAlgorithmLapackPotrf();

private:

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, GolubWelschAlgorithmLapackPotrf);

};


class GolubWelschAlgorithmLapackPotrf2 
    : 
    public GolubWelschAlgorithm<CholeskyDecompositionLapackPotrf2>
{

public:

    /**
     * @brief Construct a new `GolubWelschAlgorithmLapackPotrf2` object
     * 
     * @param nMoments Number of moments.
     */
    GolubWelschAlgorithmLapackPotrf2(int nMoments);


    /**
     * @brief Destroy the `GolubWelschAlgorithmLapackPotrf2` object
     * 
     */
    ~GolubWelschAlgorithmLapackPotrf2();

private:

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, GolubWelschAlgorithmLapackPotrf2);

};


class GolubWelschAlgorithmEigenlib 
    : 
    public GolubWelschAlgorithm<CholeskyDecompositionEigenlib>
{

public:

    /**
     * @brief Construct a new `GolubWelschAlgorithmEigenlib` object
     * 
     * @param nMoments Number of moments.
     */
    GolubWelschAlgorithmEigenlib(int nMoments);


    /**
     * @brief Destroy the `GolubWelschAlgorithmEigenlib` object
     * 
     */
    ~GolubWelschAlgorithmEigenlib();

private:

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, GolubWelschAlgorithmEigenlib);

};


class GolubWelschAlgorithmPlainCxx 
    : 
    public GolubWelschAlgorithm<CholeskyDecompositionPlainCxx>
{

public:

    /**
     * @brief Construct a new `GolubWelschAlgorithmPlainCxx` object
     * 
     * @param nMoments Number of moments.
     */
    GolubWelschAlgorithmPlainCxx(int nMoments);


    /**
     * @brief Destroy the `GolubWelschAlgorithmPlainCxx` object
     * 
     */
    ~GolubWelschAlgorithmPlainCxx();

private:

    static inline bool registered_ = REGISTER_TYPE(CoreInversionAlgorithmFactory, GolubWelschAlgorithmPlainCxx);

};

#endif // CORE_INVERSION_HPP