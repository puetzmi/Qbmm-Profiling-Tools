/**
 * @file moment_utils.hpp
 * @author Michele Puetz
 * @brief This file contains basic functions related to moments, moment realizability, random moments etc.
 * @date 2022-06-28
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifndef MOMENT_UTILS_HPP
#define MOMENT_UTILS_HPP

#include "global_defs.hpp"
#include <random>
#include <mkl.h>
#include <eigen3/Eigen/Eigen>
#include "eigen_defs.hpp"


/**
 * @brief Compute moment sequence up to given order from quadrature.
 * 
 * @param[in] nodes The quadrature nodes.
 * @param[in] weights The quadrature weights.
 * @param[in] numberOfNodes Number of quadrature nodes.
 * @param[in] numberOfMoments Number of moments, which is the maximum moment order + 1.
 * @param[out] moments The computed moment sequence.
 */
void computeMomentsFromQuadrature(double *nodes, double *weights, int numberOfNodes, 
                                    int numberOfMoments, double *moments)
{

    for (int k=0; k<numberOfMoments; k++) {
        moments[k] = 0.;
    }

    for (int j=0; j<numberOfNodes; j++) {
        double prod = weights[j];
        for (int k=0; k<numberOfMoments; k++) {
            moments[k] += prod;
            prod *= nodes[j];
        }
    }

}


/**
 * @brief Compute a moment sequence of given even length from tridiagonal Jacobi matrix.
 * 
 * The @f$ k @f$th moment @f$ m_k @f$ is computed from the @f$ k @f$th matrix power of the 
 * Jacobi matrix @f$ \mathbf{J} \in \mathbb{R}^{n \times n} @f$ as follows:
 *   
 * @f\[
 *  m_k = \mathbf{e}_1^{T} \mathbf{J}^k \mathbf{e}_1,
 * @f\]
 * where @f$ \mathbf{e}_1 @f$ is the first @f$ n @f$-dimensional standard unit vector.
 * 
 * @param[in] mainDiagonal Main diagonal of the Jacobi matrix.
 * @param[in] firstDiagonal Sub-/superdiagonal of the Jacobi matrix.
 * @param[in] nMoments Length of the moment sequence.
 * @param[out] moments The computed moment sequence (allocated space must be at least `nMoments`) .
 * @param[in] m0 Zeroth moment (scaling factor of the moment sequence).
 * @return int Dummy error flag.
 */
int computeMomentsFromJacobiMatrix(double *mainDiagonal, double *firstDiagonal, int nMoments, double *moments, double m0=1.)
{

    // Due to this, the function works only for even `nMoments` at the moment
    int jacobiMatrixSize = nMoments/2;

    // Eigen::Map objects and matrices for the simple reconstruction of moments
    VectorMap mainDiagonalMap(mainDiagonal, jacobiMatrixSize);
    VectorMap firstDiagonalMap(firstDiagonal, jacobiMatrixSize - 1);
    Eigen::MatrixXd jacobiMatrix(Eigen::MatrixXd::Zero(jacobiMatrixSize, jacobiMatrixSize));
    Eigen::MatrixXd jacobiMatrixPower(m0 * Eigen::MatrixXd::Identity(jacobiMatrixSize, jacobiMatrixSize));

    jacobiMatrix.diagonal() = mainDiagonalMap;
    jacobiMatrix.diagonal(-1) = firstDiagonalMap;
    jacobiMatrix.diagonal(1) = firstDiagonalMap;

    // The kth moment is the very first upper left element of the kth power of the Jacobi matrix
    for (int k=0; k<nMoments; k++) {
        moments[k] = jacobiMatrixPower(0,0);
        jacobiMatrixPower = jacobiMatrixPower*jacobiMatrix;
    }

    return 0;
}

#endif // MOMENT_UTILS_HPP