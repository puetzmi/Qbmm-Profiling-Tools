/**
 * @file test_eigensolve.cpp
 * @author M. Puetz
 * @brief Test of all implemented eigen solver classes that have been registered to the RealEigenSolver object factory.
 * @date 2022-09-17
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "global_defs.hpp"
#include "eigensolve.hpp"
#include <random>
#include <mkl.h>
#include <memory>
#include <cassert>

// Make sure assertions are checked
#ifdef NDEBUG
    #undef NDEBUG
#endif


/**
 * @brief Test correctness of eigenvalues and -vectors for tridiagonal symmetric matrix.
 * 
 * @param n Matrix size.
 * @param eigenSolverTypeName Name of the EigenSolver.
 * @param mainDiag Main diagonal.
 * @param firstDiag First upper/lower diagonal.
 * @param absTol Absolute tolerance for 'numerical equality' checks.
 */
void testAllEigenValues(int n, const std::string& eigenSolverTypeName,
        double *mainDiag, double *firstDiag, [[maybe_unused]] double absTol)
{

    [[maybe_unused]] int info;
    // Test computation of eigenvalues and eigenvectors
    auto eigenSolver = RealEigenSolver::makeUnique(eigenSolverTypeName, n, EigenProblemType::EigenPairs);
    auto eigenVals = eigenSolver->eigenValues();
    auto eigenVecs = eigenSolver->eigenVectors();


    // If size is > 2 test twice, once with decreased size to make sure that works
    int end = 1;
    if (n > 2) {
        n--;
        eigenSolver->setSize(n);
        end = 2;
    }
    for (int j=0; j<end; j++) {

        info = eigenSolver->compute(mainDiag, firstDiag);
        assert(info == 0);

        eigenVals = eigenSolver->eigenValues();
        eigenVecs = eigenSolver->eigenVectors();

        // As the matrix, say A, is symmetric, the equation 'A = V * lambda * V^T' must be satisfied,
        // where 'V' is the matrix of eigenvectors and lambda is the diagonal matrix of eigenvalues
        const double *V = eigenVecs.get();      // eigenvectors
        double *lambda = new double[n*n]();     // diagonal matrix of eigenvalues
        for (int i=0; i<n; i++) {
            lambda[i*n + i] = eigenVals[i];
        }
        double *A = new double[n*n]();          // result

        // A <- V*lambda
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1, V, n, lambda, n, 0, A, n);

        // lambda <- A*V^T
        for (int i=0; i<n; i++) {
            lambda[i*n + i] = 0.;
        }
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, n, 1, A, n, V, n, 0, lambda, n);

        // A <- lambda
        delete[] A;
        A = lambda;

        // Ensure approximate equality with given tolerance
        for (int i=0; i<n; i++) {
            for (int j=0; j<n; j++) {
                [[maybe_unused]] double aijRef = 0.;
                switch (i-j) {
                    case 0:     //main diagonal
                        aijRef = mainDiag[i];
                        break;
                    case 1:     // first lower diagonal
                        aijRef = firstDiag[j];
                        break;
                    case -1:    // first upper diagonal
                        aijRef = firstDiag[i];
                        break;
                }

                assert(std::fabs(aijRef - A[i*n + j]) < absTol);
            }
        }
        delete[] lambda;

        eigenSolver->setSize(n);
    }

    eigenSolver->setSize(n);

    // Now compute only eigenvalues to make sure that there are no memory-related issues
    auto eigenSolver1 = RealEigenSolver::makeUnique(eigenSolverTypeName, n, EigenProblemType::EigenValsOnly);

    info = eigenSolver1->compute(mainDiag, firstDiag);
    assert(info == 0);

    auto eigenVals1 = eigenSolver1->eigenValues();

    // Check if that yields the equal eigenvalues
    for (int i=0; i<n; i++) {
        assert(std::fabs(eigenVals[i] - eigenVals1[i]) < absTol);
    }

}

// TODO: Test computation of subsets as well

int main () {

    double absTol = 1e-12;
    const int randomSeed = DEFAULT_RANDOM_SEED;
    std::mt19937_64 randomGen;
    randomGen.seed(randomSeed);

    // Normal distribution with mean 0 and standard deviation 1
    double mu = 0.;
    double sigma = 1.;
    std::normal_distribution<double> normal(mu, sigma);

    int nMin = 2;
    int nMax = 20;      // Matrices are small when using QMOM

    // Tridiagonal symmetric matrix with independently normally distributed elements
    double *mainDiag = static_cast<double*>(mkl_malloc(nMax*sizeof(double), MALLOC_ALIGN));
    for (int i=0; i<nMax; i++) {
        mainDiag[i] = normal(randomGen);
    }
    double *firstDiag = static_cast<double*>(mkl_malloc((nMax-1)*sizeof(double), MALLOC_ALIGN));
    for (int i=0; i<nMax-1; i++) {
        firstDiag[i] = normal(randomGen);
    }

    // Test every selectable (non-abstract) eigen solver type
    const std::vector<std::string>& allEigenSolverTypeNames = RealEigenSolverFactory::keys();
    for (auto& eigenSolverTypeName : allEigenSolverTypeNames) {
        // Test for all n
        for (int n=nMin; n<nMax+1; n++) {
                testAllEigenValues(n, eigenSolverTypeName, mainDiag, firstDiag, absTol);
        }
    }

    mkl_free(mainDiag);
    mkl_free(firstDiag);

    return 0;
}
