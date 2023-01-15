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

#include <random>
#include <mkl.h>
#include <eigen3/Eigen/Eigen>
#include "global_defs.hpp"
#include "eigen_defs.hpp"
#include "constants.hpp"


/**
 * @brief Compute moment sequence up to given order from quadrature.
 * 
 * @param[in] nodes The quadrature nodes.
 * @param[in] weights The quadrature weights.
 * @param[in] numberOfNodes Number of quadrature nodes.
 * @param[in] numberOfMoments Number of moments, which is the maximum moment order + 1.
 * @param[out] moments The computed moment sequence.
 */
void computeMomentsFromQuadrature(const double *nodes, const double *weights, int numberOfNodes, 
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


/**
 * @brief Compute Chebyshev moment of the first kind.
 * 
 * @param order Moment order.
 * @return double Chebyshev moment.
 */
double getChebyshevMoment(int order)
{
    if (order == 0) {
        return constants::pi;
    }

    if (order % 2 == 1)
        return 0;

    double sqrtPi = std::sqrt(constants::pi);

    return 2*sqrtPi*std::tgamma(0.5*(1 + order))/(order*std::tgamma(0.5*order));
}


/**
 * @brief Linear moment transformation given shifting distance and scaling factor.
 * 
 * @param[in] shift Shifting distance.
 * @param[in] scale Scaling factor.
 * @param[in] nMoments Number of moments.
 * @param[in] moments Original moment sequence.
 * @param[out] transformedMoments Transformed moment sequence.
 * @return int Dummy error flag.
 */
int linearMomentTransform(double shift, double scale, int nMoments, double *moments,
    double *transformedMoments)
{

    for (int k=0; k<nMoments; k++) {

        double sum = 0;
        for (int j=0; j<k+1; j++) {

            // This is a quick awfully slow solution, but ok for mere moment
            // initializations
            int binomialCoefficient = static_cast<int>(
                std::tgamma(k+1)/(std::tgamma(j+1)*std::tgamma(k-j+1)) + 0.1
            );
            
            sum += binomialCoefficient
                   *std::pow(scale,  k - j)
                   *std::pow(shift, j)
                   *moments[k-j];
        }
        transformedMoments[k] = sum;
    }

    return 0;
}


#endif // MOMENT_UTILS_HPP
