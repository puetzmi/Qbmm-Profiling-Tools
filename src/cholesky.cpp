/**
 * @file cholesky.cpp
 * @author M. Puetz
 * @brief Implementation of Cholesky decomposition methods.
 * @date 2022-09-20
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "cholesky.hpp"
#include <mkl.h>
#include <eigen3/Eigen/Cholesky>


int CholeskyDecompositionPlainCxx::compute(double *matrix)
{

    double *ithRowPtr = matrix;

    for (int i=0; i<size_; i++) {

        double *jthRowPtr = matrix;
        for (int j=0; j<=i; j++) {

            double sum = 0;
            for (int k=0; k<j; k++) {
                sum += ithRowPtr[k] * jthRowPtr[k];
            }

            if (i == j) {
                double d = ithRowPtr[i] - sum;
                if (d <= 0) return i;
                ithRowPtr[i] = std::sqrt(d);
            }
            else {
                ithRowPtr[j] = 1.0 / jthRowPtr[j] * (ithRowPtr[j] - sum);
            }
            
            jthRowPtr += size_;
        }
        ithRowPtr += size_;
    }

    return 0;
}


int CholeskyDecompositionLapack::compute(double *matrix)
{

    int info = lapackFunction_(LAPACK_ROW_MAJOR, 'L', size_, matrix, size_);

    return info;
}


int CholeskyDecompositionEigenlib::compute(double *matrix)
{

    if (matrix != matrixMap_.data())
    {
        // This is the recommended way to reassign data (see Eigen documentation)
        new (&matrixMap_) MatrixMap(matrix, size_, size_);
    }

    // This is not the nicest imaginable way but it seems to be the only way since
    // "external" Eigen functions and classes like `LLT` do not return any information
    // on non-positive minors
    int info = Eigen::internal::llt_inplace<double,Eigen::UpLoType::Lower>::blocked(matrixMap_);
    return (info == -1) ? 0 : info;
}


CholeskyDecomposition::CholeskyDecomposition(int maxSize)
    :
    maxSize_(maxSize),
    size_(maxSize)
{
}


CholeskyDecomposition::~CholeskyDecomposition()
{
}


std::unique_ptr<CholeskyDecomposition> CholeskyDecomposition::makeUnique(
        const std::string& typeName, int maxSize)
{
    return CholeskyDecompositionFactory::makeUnique(typeName, maxSize);
}


std::shared_ptr<CholeskyDecomposition> CholeskyDecomposition::makeShared(
        const std::string& typeName, int maxSize)
{
    return CholeskyDecompositionFactory::makeShared(typeName, maxSize);
}


CholeskyDecompositionLapack::CholeskyDecompositionLapack(
    int maxSize, CholeskyDecompositionLapack::MklLapackCholeskyFunction lapackFunction)
    :
    CholeskyDecomposition(maxSize),
    lapackFunction_(lapackFunction)
{
}


CholeskyDecompositionLapack::~CholeskyDecompositionLapack()
{
}


CholeskyDecompositionLapackPotrf::CholeskyDecompositionLapackPotrf(int maxSize)
    :
    CholeskyDecompositionLapack(maxSize, LAPACKE_dpotrf)
{
}


CholeskyDecompositionLapackPotrf2::~CholeskyDecompositionLapackPotrf2()
{
}


CholeskyDecompositionLapackPotrf2::CholeskyDecompositionLapackPotrf2(int maxSize)
    :
    CholeskyDecompositionLapack(maxSize, LAPACKE_dpotrf2)
{
}


CholeskyDecompositionLapackPotrf::~CholeskyDecompositionLapackPotrf()
{
}


CholeskyDecompositionEigenlib::CholeskyDecompositionEigenlib(int maxSize)
    :
    CholeskyDecomposition(maxSize),
    matrixMap_(nullptr, size_, size_)
{
}


CholeskyDecompositionEigenlib::~CholeskyDecompositionEigenlib()
{
}


CholeskyDecompositionPlainCxx::CholeskyDecompositionPlainCxx(int maxSize)
    :
    CholeskyDecomposition(maxSize)
{
}

CholeskyDecompositionPlainCxx::~CholeskyDecompositionPlainCxx()
{
}