/**
 * @file test_linsolve.cpp
 * @author M. Puetz
 * @brief Test of all implemented (and registered) linear solvers using a random Vandermonde system.
 * @date 2022-09-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "global_defs.hpp"
#include "linsolve.hpp"
#include <random>
#include <mkl.h>
#include <vector>

// Make sure assertions are checked even in non-debug mode
#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>


/**
 * @brief Test linear solver of given type using a random Vandermonde system of given size.
 * 
 * @param linearSolverPtr Shared pointer to linear solver object.
 * @param n Size of the linear system.
 * @param randomGen Pseudo-random number generator.
 * @param absTol Absolute tolerance used for comparison between results and reference.
 */
void test(std::shared_ptr<LinearSolver> linearSolverPtr, int n, std::mt19937& randomGen, double absTol) {

    std::normal_distribution<double> normal{};

    double *aPtr = static_cast<double*>(mkl_malloc(n*n*sizeof(double), MALLOC_ALIGN));
    double *xPtr = static_cast<double*>(mkl_malloc(n*sizeof(double), MALLOC_ALIGN));
    double *xRefPtr = static_cast<double*>(mkl_malloc(n*sizeof(double), MALLOC_ALIGN));
    double *bPtr = static_cast<double*>(mkl_malloc(n*sizeof(double), MALLOC_ALIGN));
     
    // Maps mapping Eigen data types to raw pointers and vice versa
    MatrixMap A(aPtr, n, n);
    VectorMap x(xPtr, n);
    VectorMap xRef(xRefPtr, n);
    VectorMap b(bPtr, n);

    // Successive solution of two systems to make sure the solver is updated correctly
    for (int k=0; k<2; k++) {

        // Generate random Vandermonde matrix and solution vector
        for (int j=0; j<n; j++) {
            xRefPtr[j] = normal(randomGen);
            aPtr[j] = 1.;
            double aj = normal(randomGen);
            for (int i=1; i<n; i++) {
                aPtr[i*n + j] = aPtr[(i-1)*n + j]*aj;
            }
        }

        // R.h.s vector
        b = A*xRef;

        // Solve system
        linearSolverPtr->solve(aPtr, bPtr, xPtr);

        // Check results with given tolerance
        for (int i=0; i<n; i++) {
            assert(std::fabs(xPtr[i] - xRefPtr[i]) < absTol);
        }
    }

    mkl_free(aPtr);
    mkl_free(bPtr);
    mkl_free(xPtr);
    mkl_free(xRefPtr);

}


int main () 
{

    // Parameters
    int nMin = 2;
    int nMax = 10;
    double absTol = 1e-7;


    const int seed = DEFAULT_RANDOM_SEED;
    std::mt19937 randomGen;
    randomGen.seed(seed);

    // Test for all registered types and given sizes
    const std::vector<std::string> &allLinearSolverTypeNames = LinearSolverFactory::keys();
    for (auto linearSolverTypeName : allLinearSolverTypeNames) {
        for (int n=nMin; n<=nMax; n++) {

            auto linearSolverPtr = LinearSolver::makeShared(linearSolverTypeName, n);
            test(linearSolverPtr, n, randomGen, absTol);

            // test if decreasing the size works
            linearSolverPtr->setSize(n-1);
            test(linearSolverPtr, n-1, randomGen, absTol);
        }
    }

    return 0;
}
