/**
 * @file test_cholesky.cpp
 * @author M. Puetz
 * @brief Test functionality of all core inversion classes registered to the object factory.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "global_defs.hpp"
#include "eigen_defs.hpp"
#include "core_inversion.hpp"
#include <cstdio>
#include <cmath>
#include <random>
#include <functional>
#include <eigen3/Eigen/Eigen>

#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>


int main() {

    // Minumum and maximum size of Jacobi matrix
    int nMin = 2;
    int nMax = 10;

    double relTol = 1e-5;
    const int randomSeed = DEFAULT_RANDOM_SEED;
    std::mt19937_64 randomGen;
    randomGen.seed(randomSeed);

    // Normal distribution with mean 0 and standard deviation 1
    double mu = 0.;
    double sigma = 1.;
    std::normal_distribution<double> normalDistribution(mu, sigma);

    // Gamma distribution with shape parameter 1
    double k = 1;
    std::gamma_distribution<double> gammaDistribution(k);

    // Allocate and initialize arrays for tridiagonal symmetric Jacobi matrix
    double *jacobiMainDiagonal = static_cast<double*>(mkl_malloc(nMax*sizeof(double), MALLOC_ALIGN));
    double *jacobiSubDiagonal = static_cast<double*>(mkl_malloc((nMax-1)*sizeof(double), MALLOC_ALIGN));
    jacobiMainDiagonal[0] = normalDistribution(randomGen);
    for (int i=1; i<=nMax; i++) {
        jacobiMainDiagonal[i] = normalDistribution(randomGen);
        jacobiSubDiagonal[i-1] = std::sqrt(gammaDistribution(randomGen));
    }

    // Allocate space for moments and computed coefficients
    double *moments = static_cast<double*>(mkl_malloc(nMax*2*sizeof(double), MALLOC_ALIGN));
    double *alpha = static_cast<double*>(mkl_malloc(nMax*sizeof(double), MALLOC_ALIGN));
    double *gamma = static_cast<double*>(mkl_malloc((nMax - 1)*sizeof(double), MALLOC_ALIGN));

    // Loop over all available types
    const std::vector<std::string>& allCoreInversionTypeNames = CoreInversionAlgorithmFactory::keys();
    for (auto& coreInversionTypeName : allCoreInversionTypeNames) {

        for (int n=nMin; n<=nMax; n++) {

            int nMoments = n*2;

            // Compute moments from Jacobi matrix
            VectorMap jacobiMainDiagonalMap(jacobiMainDiagonal, n);
            VectorMap jacobiSubDiagonalMap(jacobiSubDiagonal, n - 1);
            Eigen::MatrixXd jacobiMatrix(Eigen::MatrixXd::Zero(n, n));
            jacobiMatrix.diagonal() = jacobiMainDiagonalMap;
            jacobiMatrix.diagonal(-1) = jacobiSubDiagonalMap;
            jacobiMatrix.diagonal(1) = jacobiSubDiagonalMap;
            Eigen::MatrixXd jacobiMatrixPower(Eigen::MatrixXd::Identity(n, n));
            for (int k=0; k<nMoments; k++) {
                moments[k] = jacobiMatrixPower(0,0);
                jacobiMatrixPower *= jacobiMatrix;
            }

            // Create object and compute
            auto coreInversionPtr = CoreInversionAlgorithm::makeShared(coreInversionTypeName, nMoments);
            coreInversionPtr->compute(moments, alpha, gamma);

            // Check result
            for (int k=0; k<n-1; k++) {
                assert(std::abs((alpha[k] - jacobiMainDiagonal[k])/jacobiMainDiagonal[k]) < relTol);
                assert(std::abs((gamma[k] - jacobiSubDiagonal[k])/jacobiSubDiagonal[k]) < relTol);
            }
            assert(std::abs((alpha[n-1] - jacobiMainDiagonal[n-1])/jacobiMainDiagonal[n-1]) < relTol);
        } 
    }

    mkl_free(jacobiMainDiagonal);
    mkl_free(jacobiSubDiagonal);
    mkl_free(moments);
    mkl_free(alpha);
    mkl_free(gamma);

    return 0;
}
