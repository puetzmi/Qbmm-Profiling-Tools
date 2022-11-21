/**
 * @file test_cholesky.cpp
 * @author M. Puetz
 * @brief Test functionality of all implemented Cholesky decomposition methods of which the respective classes are registered to the object factory.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "global_defs.hpp"
#include "cholesky.hpp"
#include <cstdio>
#include <cmath>
#include <random>
#include <functional>

#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>


void testCholesky(std::shared_ptr<CholeskyDecomposition> choleskyDecompPtr, 
                    const double *matrix, int n, double *resultMatrix, double *workMatrix, double absTol)
{

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            resultMatrix[i*n + j] = matrix[i*n + j];
            workMatrix[i*n + j] = matrix[i*n + j];
        }
    }

    int info = choleskyDecompPtr->compute(resultMatrix);
    info += 0;
    assert(info == 0);

    // Make lower triangular
    for (int i=0; i<n; i++) {
        for (int j=0; j<i; j++) {
            resultMatrix[j*n + i] = 0;
        }
    }

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1, 
                resultMatrix, n, resultMatrix, n, -1, workMatrix, n);
    
    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            assert(std::fabs(workMatrix[i*n + j]) < absTol);
        }
    }

}


void prepareMatrices(int n, double *workMatrix, double *resultMatrix, double *spdMatrix, 
                        const std::function<double()>& randomSample) {

    for (int i=0; i<n; i++) {
        for (int j=0; j<n; j++) {
            workMatrix[i*n + j] = randomSample();
            spdMatrix[i*n + j] = 0.;
            resultMatrix[i*n + j] = 0.;
        }
    }

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n, n, n, 1, 
                workMatrix, n, workMatrix, n, 0, spdMatrix, n);
    
}

int main() {

    double absTol = 1e-12;
    int randomSeed = 1120012;
    std::mt19937_64 randomGen;
    randomGen.seed(randomSeed);

    // Normal distribution with mean 0 and standard deviation 1
    double mu = 0.;
    double sigma = 1.;
    std::normal_distribution<double> normal(mu, sigma);
    std::function<double()> sample(std::bind(normal, randomGen));

    int nMin = 2;
    int nMax = 40;

    double *workMatrix = static_cast<double*>(mkl_malloc(nMax*nMax*sizeof(double), MALLOC_ALIGN));
    double *resultMatrix = static_cast<double*>(mkl_malloc(nMax*nMax*sizeof(double), MALLOC_ALIGN));
    double *spdMatrix = static_cast<double*>(mkl_malloc(nMax*nMax*sizeof(double), MALLOC_ALIGN));


    const std::vector<std::string>& allCholeskyDecompositionTypeNames = CholeskyDecompositionFactory::keys();
    for (auto& choleskyDecompositionTypeName : allCholeskyDecompositionTypeNames) {

        for (int n=nMin; n<=nMax; n++) {

            prepareMatrices(n, workMatrix, resultMatrix, spdMatrix, sample);
            auto choleskyDecompPtr = CholeskyDecomposition::makeShared(choleskyDecompositionTypeName, n);
            testCholesky(choleskyDecompPtr, spdMatrix, n, resultMatrix, workMatrix, absTol);
        } 
    }

    mkl_free(workMatrix);
    mkl_free(resultMatrix);
    mkl_free(spdMatrix);

    return 0;
}
