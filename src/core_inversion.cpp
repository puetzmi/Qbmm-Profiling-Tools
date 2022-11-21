/**
 * @file core_inversion.cpp
 * @author M. Puetz
 * @brief Implementation of algorithms computing the tridiagonal Jacobi matrix from moments.
 * @date 2022-10-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "core_inversion.hpp"
#include <mkl.h>
#include <stdexcept>
#include <string>
#include <limits>


CoreInversionAlgorithm::CoreInversionAlgorithm(int nMoments, int nWorkAlloc)
    :
    nMoments_ (nMoments),
    work_ (static_cast<double*>(mkl_malloc(nWorkAlloc*sizeof(double), MALLOC_ALIGN)))
{
    if (nMoments_ < 2 || (nMoments % 2) != 0) {
        std::string msg = "Found invalid argument: `nMom` must be an even positive integer.";
        throw std::invalid_argument(msg);
    }
}


CoreInversionAlgorithm::~CoreInversionAlgorithm() {
    mkl_free(work_);
}


std::unique_ptr<CoreInversionAlgorithm> 
    CoreInversionAlgorithm::makeUnique(const std::string &typeName, int nMoments) 
{
    return CoreInversionAlgorithmFactory::makeUnique(typeName, nMoments);
}


std::shared_ptr<CoreInversionAlgorithm> 
    CoreInversionAlgorithm::makeShared(const std::string &typeName, int nMoments) 
{
    return CoreInversionAlgorithmFactory::makeUnique(typeName, nMoments);
}


int CoreInversionAlgorithm::operator()(double *moments, double *alpha, double *gamma) {
    return this->compute(moments, alpha, gamma);
}


LqmdAlgorithm::LqmdAlgorithm(int nMoments) 
    :
    CoreInversionAlgorithm(nMoments, nMoments*nMoments/2)
{
}


LqmdAlgorithm::~LqmdAlgorithm() {}


int LqmdAlgorithm::compute(double *moments, double *alpha, double *gamma) {

    double *sigma = work_;
    int n = nMoments_/2;
    int nMinus1 = n - 1;
    int twoNMinus1 = nMoments_ - 1;

    // Set first row of auxiliary matrix sigma to normalized moments
    for (int j=0; j<nMoments_; j++) {
        sigma[j] = moments[j]/moments[0];
    }

    // Start of algorithm, i = 0
    alpha[0] = sigma[1];
    gamma[0] =  sigma[2] - alpha[0]*sigma[1];
    for (int j=2; j<twoNMinus1; j++) {
        sigma[nMoments_ + j] = (sigma[j+1] - alpha[0]*sigma[j]) / gamma[0];
    }
    // Square roots of recurrence coefficients yield sub-/superdiagonal elements of Jacobi matrix
    gamma[0] = std::sqrt(gamma[0]);

    // Pointers to (i-1)th, ith and (i+1)th row of auxiliary matrix sigma, respectively
    double *sigmaIMinus1;
    double *sigmaI = sigma;
    double *sigmaIPlus1 = &(sigma[nMoments_]);

    // Recursion for i=1..n-2
    for (int i=1; i<nMinus1; i++) {
        sigmaIMinus1 = sigmaI;
        sigmaI = sigmaIPlus1;
        sigmaIPlus1 += nMoments_;
        alpha[i] = sigmaI[i+1] - sigmaIMinus1[i];
        gamma[i] =  sigmaI[i+2] - (alpha[i]*sigmaI[i+1] + sigmaIMinus1[i+1]);
        for (int j=2; j<twoNMinus1; j++) {
            if (gamma[i] == 0) break;   // see explanation below
            sigmaIPlus1[j] = (sigmaI[j+1] - (alpha[i]*sigmaI[j] + sigmaIMinus1[j])) / gamma[i];
        }

        // Square roots of recurrence coefficients yield sub-/superdiagonal elements of Jacobi matrix;
        // if a non-positive value is encountered the Jacobi matrix corresponds to a unrealizable moment set,
        // in that case, all following coefficents are set to zero, the algorithm is terminated and the current i
        // is returned.
        if (gamma[i] <= 0) {
            for (int k=i; i<nMinus1; i++) {
                gamma[k] = 0.;
                alpha[k] = 0.;
            }
            alpha[nMinus1] = 0.;
            return i;
        }
        gamma[i] = std::sqrt(gamma[i]);
    }

    // Last recurrence coefficient / diagonal element in Jacobi matrix
    alpha[nMinus1] = sigma[n*twoNMinus1] - sigma[n*(nMoments_ - 3) - 1];

    return 0;
}


template <class CholeskyDecompositionType>
GolubWelschAlgorithm<CholeskyDecompositionType>::GolubWelschAlgorithm(int nMoments) 
    :
    CoreInversionAlgorithm(nMoments, nMoments*nMoments/2),
    choleskyDecomposition_(new CholeskyDecompositionType(nMoments/2+1))
{
}


template <class CholeskyDecompositionType>
GolubWelschAlgorithm<CholeskyDecompositionType>::~GolubWelschAlgorithm() {}


template <class CholeskyDecompositionType>
int GolubWelschAlgorithm<CholeskyDecompositionType>::compute(double *moments, double *alpha, double *gamma) {

    double *matrixL = work_;
    int n = nMoments_/2;
    int nPlus1 = n + 1;
    constexpr double maxDouble = std::numeric_limits<double>::max();

    // Assemble Hankel matrix (only upper triangular part is required for the Cholesky decomposition)
    for (int i=0; i<nPlus1; i++) {
        for (int j=0; j<i; j++) {
            matrixL[i*nPlus1+j] = moments[i+j];
        }
    }
    for (int i=0; i<n; i++) {
        matrixL[i*(nPlus1 + 1)] = moments[2*i];
    }
    // This should suffice to make Hankel matrix positive definite if nothing is wrong with the moment set
    // (only needed to enable Cholesky decomposition by third-party libraries, not required to assemble the Jacobi matrix)
    matrixL[nPlus1*nPlus1 - 1] = 0.1*maxDouble;

    // Cholesky decomposition
    int info = choleskyDecomposition_->compute(matrixL);

    // Calculation of Jacobi matrix
    // for j = 0
    double ljj = matrixL[0];
    double ljjInv = 1/ljj;
    double c = matrixL[nPlus1]*ljjInv;
    alpha[0] = c;
    // ... and for j={1,..,N-1}
    int nSpd = (info == 0) ? n : info - 1;
    for (int j=1; j<nSpd; j++) {
        alpha[j] = -c;
        ljj = matrixL[nPlus1*j + j];
        gamma[j-1] = ljj*ljjInv;      // Here ljjI has still the value from previous iteration j - 1
        ljjInv = 1/ljj;
        c = matrixL[nPlus1*(j + 1) + j]*ljjInv;
        alpha[j] += c;
    }
    for (int j=nSpd; j<n; j++) {
        alpha[j] = 0;
        gamma[j-1] = 0;
    }

    return info;                    // Should be zero if nothing went wrong in Cholesky decomposition
}


GolubWelschAlgorithmLapackPotrf::GolubWelschAlgorithmLapackPotrf(int nMoments)
    :
    GolubWelschAlgorithm<CholeskyDecompositionLapackPotrf>(nMoments)
{
}

GolubWelschAlgorithmLapackPotrf::~GolubWelschAlgorithmLapackPotrf()
{
}


GolubWelschAlgorithmLapackPotrf2::GolubWelschAlgorithmLapackPotrf2(int nMoments)
    :
    GolubWelschAlgorithm<CholeskyDecompositionLapackPotrf2>(nMoments)
{
}

GolubWelschAlgorithmLapackPotrf2::~GolubWelschAlgorithmLapackPotrf2()
{
}


GolubWelschAlgorithmEigenlib::GolubWelschAlgorithmEigenlib(int nMoments)
    :
    GolubWelschAlgorithm<CholeskyDecompositionEigenlib>(nMoments)
{
}

GolubWelschAlgorithmEigenlib::~GolubWelschAlgorithmEigenlib()
{
}


GolubWelschAlgorithmPlainCxx::GolubWelschAlgorithmPlainCxx(int nMoments)
    :
    GolubWelschAlgorithm<CholeskyDecompositionPlainCxx>(nMoments)
{
}

GolubWelschAlgorithmPlainCxx::~GolubWelschAlgorithmPlainCxx()
{
}