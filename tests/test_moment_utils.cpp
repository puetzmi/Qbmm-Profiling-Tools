/**
 * @file test_moment_utils.cpp
 * @author M. Puetz
 * @brief Test of functions in moment utils.
 * @date 2022-12-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <stdexcept>
#include <cmath>
#include <mkl.h>

#include "constants.hpp"
#include "global_defs.hpp"
#include "moment_utils.hpp"


void testMomentUtils()
{

    // Parameters
    int nMoments = 12;
    double absTol = 1e-10;
    double relTol = 1e-8;
    double transformShift = -1.2;
    double transformScale = 1.23;

    // Chebyshev nodes (first kind) and weights
    int nNodes = (nMoments + 1)/2;
    double *quadratureNodes = static_cast<double*>
    (
        mkl_malloc(nNodes*sizeof(double), MALLOC_ALIGN)
    );
    double *quadratureWeights = static_cast<double*>
    (
        mkl_malloc(nNodes*sizeof(double), MALLOC_ALIGN)
    );
    for (int j=1; j<=nNodes; j++) {
        quadratureNodes[j-1] = std::cos(0.5*(2*j - 1)/nNodes*constants::pi);
        quadratureWeights[j-1] = constants::pi/nNodes;
    }

    // Test `getChebyshevMoment` and `computeMomentsFromQuadrature`
    // using standard Chebyshev moments of the first kind.
    double *moments = static_cast<double*>
    (
        mkl_malloc(nMoments*sizeof(double), MALLOC_ALIGN)
    );
    computeMomentsFromQuadrature(quadratureNodes, quadratureWeights, 
        nNodes, nMoments, moments);
    for (int k=0; k<nMoments; k++) {
        assert(std::abs(moments[k] - getChebyshevMoment(k)) 
            < absTol + relTol*std::abs(moments[k]));
    }

    // Compute transformed moments and verify results with moments
    // from shifted nodes and weights
    double *transformedMoments = static_cast<double*>
    (
        mkl_malloc(nMoments*sizeof(double), MALLOC_ALIGN)
    );
    linearMomentTransform(transformShift, transformScale,
        nMoments, moments, transformedMoments);
    for (int j=0; j<nNodes; j++) {
        quadratureNodes[j] *= transformScale;
        quadratureNodes[j] += transformShift;
    }
    computeMomentsFromQuadrature(quadratureNodes, quadratureWeights, 
        nNodes, nMoments, moments);
    for (int k=0; k<nMoments; k++) {
        assert(std::abs(moments[k] - transformedMoments[k]) 
            < absTol + relTol*std::abs(moments[k]));
    }

    // Assemble Jacobi matrix for Chebyshev polynomials of the first kind
    // (see Gautschi2004, Table 1.1) and verify result from 
    // `computeMomentsFromJacobiMatrix`.
    double *mainDiagonal = static_cast<double*>
    (
        mkl_malloc(nNodes*sizeof(double), MALLOC_ALIGN)
    );
    double *superDiagonal = static_cast<double*>
    (
        mkl_malloc((nNodes-1)*sizeof(double), MALLOC_ALIGN)
    );
    for (int j=0; j<nNodes; j++) {
        mainDiagonal[j] = 0.;
    }
    superDiagonal[0] = std::sqrt(0.5);
    for (int j=1; j<nNodes-1; j++) {
        superDiagonal[j] = 0.5;
    }

    // Modify Jacobi matrix by shifting the main diagonal elements and
    // scaling the off-diagonal elements; verify moments computed by
    // `computeMomentsFromJacobiMatrix`
    mainDiagonal[nNodes - 1] += transformShift;
    for (int j=0; j<nNodes-1; j++) {
        mainDiagonal[j] += transformShift;
        superDiagonal[j] *= transformScale;
    }
    computeMomentsFromJacobiMatrix(mainDiagonal, superDiagonal,
        nMoments, moments, getChebyshevMoment(0));
    for (int k=0; k<nMoments; k++) {
        assert(std::abs(moments[k] - transformedMoments[k]) 
            < absTol + relTol*std::abs(moments[k]));
    }


    mkl_free(quadratureNodes);
    mkl_free(quadratureWeights);
    mkl_free(moments);
    mkl_free(transformedMoments);
    mkl_free(mainDiagonal);
    mkl_free(superDiagonal);
}


int main()
{

    testMomentUtils();
    return 0;
}