/**
 * @file cholesky.hpp
 * @author M. Puetz
 * @brief Classes for Cholesky decompositions.
 * @date 2022-09-20
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef CHOLESKY_HPP
#define CHOLESKY_HPP

#include "global_defs.hpp"
#include "eigen_defs.hpp"
#include "factory.hpp"
#include <mkl.h>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Cholesky>
#include <functional>
#include <memory>


class CholeskyDecomposition;    // forward declaration

/// Factory for `CholeskyDecomposition` objects
using CholeskyDecompositionFactory = Factory<CholeskyDecomposition, int>;


/**
 * @brief Base class for Cholesky decompositions.
 * 
 * Classes of this type compute a Cholesky decomposition, i.e. the decomposition
 * @f\[
 *  \mathbf{A} = \mathbf{L} \mathbf{L}^{\mathrm{T}}
 * @f\]
 * where @f$ \mathbf{A} @f$ is a symmetric positive-definite matrix and @f$ \mathbf{L} @f$ 
 * is a lower triangular matrix.
 * 
 */
class CholeskyDecomposition
{

protected:

    const int maxSize_;                 ///< The maximum size of input matrices
    int size_;                          ///< The actual size of input matrices

public:

    /**
     * @brief Construct a new `CholeskyDecomposition` object from given maximum size.
     * 
     * @param maxSize The maximum size of input matrices.
     */
    CholeskyDecomposition(int maxSize);


    /**
     * @brief Destroy the `CholeskyDecomposition` object.
     * 
     */
    virtual ~CholeskyDecomposition();


    /**
     * @brief Compute the Cholesky decomposition of a given matrix.
     * 
     * The input matrix, given as a C-style (row-major) array, must be symmetric and positive-definite. 
     * Symmetry is not checked and only the lower triangular part of the matrix is used.
     * 
     * @param matrix Input matrix in row-major format.
     * @return int Error flag
     */
    virtual int compute(double *matrix) = 0;


    /**
     * @brief Alternative way to call the `compute` method.
     * 
     */
    int operator()(double *matrix)
    {
        return this->compute(matrix);
    }


    /**
     * @brief Create new `CholeskyDecomposition` object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `CholeskyDecompositionFactory`.
     * @param maxSize Maximum size of input matrices.
     * @return std::unique_ptr<CholeskyDecomposition> Unique pointer to new `CholeskyDecomposition` object.
     * 
     */
    static std::unique_ptr<CholeskyDecomposition> makeUnique(
            const std::string& typeName, int maxSize);


    /**
     * @brief Create new `CholeskyDecomposition` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `CholeskyDecompositionFactory`.
     * @param maxSize Maximum size of input matrices.
     * @return std::shared_ptr<CholeskyDecomposition> Shared pointer to new `CholeskyDecomposition` object.
     * 
     */
    static std::shared_ptr<CholeskyDecomposition> makeShared(
            const std::string& typeName, int maxSize);

};


/**
 * @brief Base class serving as interface to MKL/LAPACK functions for Cholesky decompositions.
 * 
 */
class CholeskyDecompositionLapack : public CholeskyDecomposition
{

    /// @brief Type for MKL/LAPACK functions performing Cholesky decompositions
    using MklLapackCholeskyFunction = std::function<int (int, char, lapack_int, double*, lapack_int)>;

public:

    /**
     * @brief Construct a new `CholeskyDecompositionLapack` object given the maximum size of input matrices
     * and the MKL/LAPACK function.
     * 
     * @param maxSize Maximum size of input matrices.
     * @param lapackFunction MKL/LAPACK function to be called to perform a Cholesky decomposition.
     */
    CholeskyDecompositionLapack(int maxSize, MklLapackCholeskyFunction lapackFunction);

    /**
     * @brief Destroy the `CholeskyDecompositionLapackobject`.
     * 
     */
    virtual ~CholeskyDecompositionLapack();

    virtual int compute(double *matrix);


protected:

    MklLapackCholeskyFunction lapackFunction_;      ///< MKL/LAPACK function to be called for computation

};


/**
 * @brief Class serving as interface to the MKL/LAPACK function 'POTRF' for a Cholesky decomposition.
 * 
 */
class CholeskyDecompositionLapackPotrf : public CholeskyDecompositionLapack
{

public:

    /**
     * @brief Construct a new `CholeskyDecompositionLapackPotrf` object given the maximum size of input matrices.
     * 
     * @param maxSize Maximum size of input matrices.
     */
    CholeskyDecompositionLapackPotrf(int maxSize);

    /**
     * @brief Destroy the `CholeskyDecompositionLapackPotrf` object.
     * 
     */
    virtual ~CholeskyDecompositionLapackPotrf();


private:

    static bool inline registered_ = 
        REGISTER_TYPE(CholeskyDecompositionFactory, CholeskyDecompositionLapackPotrf);

};


/**
 * @brief Class serving as interface to the MKL/LAPACK function 'POTRF2' for a Cholesky decomposition.
 * 
 */
class CholeskyDecompositionLapackPotrf2 : public CholeskyDecompositionLapack
{

public:

    /**
     * @brief Construct a new `CholeskyDecompositionLapackPotrf2` object given the maximum size of input matrices.
     * 
     * @param maxSize Maximum size of input matrices.
     */
    CholeskyDecompositionLapackPotrf2(int maxSize);

    /**
     * @brief Destroy the `CholeskyDecompositionLapackPotrf2` object.
     * 
     */
    virtual ~CholeskyDecompositionLapackPotrf2();


private:

    static bool inline registered_ = 
        REGISTER_TYPE(CholeskyDecompositionFactory, CholeskyDecompositionLapackPotrf2);

};


/**
 * @brief Class serving as interface for Cholesky decompositions using the Eigen library.
 * 
 */
class CholeskyDecompositionEigenlib : public CholeskyDecomposition
{

public:

    /**
     * @brief Construct a new `CholeskyDecompositionEigenlib` object.
     * 
     * @param maxSize Maximum size of input matrices.
     */
    CholeskyDecompositionEigenlib(int maxSize);

    /**
     * @brief Destroy the `CholeskyDecompositionEigenlib` object.
     * 
     */
    virtual ~CholeskyDecompositionEigenlib();

    virtual int compute(double *matrix);

private:

    MatrixMap matrixMap_;                                               ///< Eigen::Map mapping raw pointer to Matrix and vice versa.

    static bool inline registered_ = 
        REGISTER_TYPE(CholeskyDecompositionFactory, CholeskyDecompositionEigenlib);

};


/**
 * @brief Class implementing the Cholesky decomposition in plain C++.
 * 
 */
class CholeskyDecompositionPlainCxx : public CholeskyDecomposition
{

public:

    /**
     * @brief Construct a new `CholeskyDecompositionPlainCxx` object.
     * 
     * @param maxSize Maximum size of input matrices.
     */
    CholeskyDecompositionPlainCxx(int maxSize);

    /**
     * @brief Destroy the `CholeskyDecompositionPlainCxx` object.
     * 
     */
    virtual ~CholeskyDecompositionPlainCxx();

    virtual int compute(double *matrix);

private:

    static bool inline registered_ = 
        REGISTER_TYPE(CholeskyDecompositionFactory, CholeskyDecompositionPlainCxx);

};

#endif // CHOLESKY_HPP