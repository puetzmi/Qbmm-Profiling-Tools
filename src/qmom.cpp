/**
 * @file qmom.cpp
 * @author Michele Puetz
 * @brief Implementation of QMOM routines.
 * @version 0.1
 * @date 2022-05-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include "qmom.hpp"
#include "core_inversion.hpp"
#include <cmath>
#include <mkl.h>
#include "constants.hpp"
#include <stdexcept>
#include <functional>


Qmom::Qmom(int nMoments, int nNodes, std::shared_ptr<CoreInversionAlgorithm> coreInversion,
            std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver)
    :
    nMoments_ (nMoments),
    N_ (nMoments/2),
    twoN_ (2*N_),
    nNodes_ (nNodes),
    coreInversion_ (coreInversion),
    eigenSolver_ (eigenSolver),
    linearSolver_ (linearSolver),
    computeWeightsFromEigenVectors_ (!(static_cast<bool>(linearSolver_))),
    quadratureNodes_(static_cast<double*>(mkl_malloc(nNodes*sizeof(double), MALLOC_ALIGN)), &mkl_free),
    quadratureWeights_(static_cast<double*>(mkl_malloc(nNodes*sizeof(double), MALLOC_ALIGN)), &mkl_free),
    vandermondeMatrix_(static_cast<double*>(mkl_malloc(nMoments*nNodes*sizeof(double), MALLOC_ALIGN))),
    alpha_ (static_cast<double*>(mkl_malloc(N_*sizeof(double), MALLOC_ALIGN))),
    gamma_ (static_cast<double*>(mkl_malloc((N_-1)*sizeof(double), MALLOC_ALIGN))),
    eigenVecs_ (eigenSolver_->eigenVectors()),
    momentsRateOfChange_(static_cast<double*>(mkl_malloc(2*N_*sizeof(double), MALLOC_ALIGN)), &mkl_free)
{

    // Do not compute eigenvectors if linear solver is provided
    if (linearSolver_) {
        eigenSolver_->setProblemType(EigenProblemType::EigenValsOnly);
    }
    else {
        eigenSolver_->setProblemType(EigenProblemType::EigenPairs);
    }

    // First row of Vandermonde matrix is always 1
    for (int j=0; j<N_; j++) {
        vandermondeMatrix_[j] = 1.;
    }

}


Qmom::~Qmom() {

    mkl_free(vandermondeMatrix_); 
    mkl_free(alpha_);
    mkl_free(gamma_); 
}


void Qmom::updateVandermondeMatrix(int end, int begin) {

    double *Vim1;                   // (i-1)th row
    double *Vi = &(vandermondeMatrix_[begin-1]);    // ith row

    for (int i=begin; i<end; i++) {
        Vim1 = Vi;
        Vi = &(vandermondeMatrix_[i*N_]);
        for (int j=0; j<N_; j++) {
            Vi[j] = Vim1[j]*quadratureNodes_[j];
        }
    }
}


int Qmom::computeMomentsRateOfChange(int nMoments, const std::function<double(double)> &xiDot) {

    for (int k=0; k<nMoments; k++) {
        momentsRateOfChange_[k] = 0.;
    }

    for (int j=0; j<nNodes_; j++) {
        double prod = quadratureWeights_[j] * xiDot(quadratureNodes_[j]);
        for (int k=0; k<nMoments; k++) {
            momentsRateOfChange_[k] += prod;
            prod *= quadratureNodes_[j];
        }
    }

    return 0;
}

int Qmom::computeJacobiMatrix(double *moments)
{
    return coreInversion_->compute(moments, alpha_, gamma_);
}


int Qmom::compute(double *moments, const std::function<double(double)> &gFunc) {

    // Single-node case
    if (N_ == 1) {
        quadratureNodes_[0] = moments[1]/moments[0];
        quadratureWeights_[0] = moments[0];
        return 0;
    }

    int info = this->computeJacobiMatrix(moments);

    info += this->computeQuadrature(moments);

    info += this->computeMomentsRateOfChange(nMoments_, gFunc);

    return info;
}

int Qmom::computeMoments(int kmin, int kmax, double *mom) {

    for (int k=0; k<twoN_; k++) {
        mom[k] = 0.;
    }

    for (int j=0; j<nNodes_; j++) {
        double prod = quadratureWeights_[j];
        for (int k=kmin; k<kmax; k++) {
            mom[k] += prod;
            prod *= quadratureNodes_[j];
        }
    }

    return 0;

}

void Qmom::setJacobiMatrix(double *mainDiagonal, double *firstDiagonal)
{
    for (int i=0; i<N_-1; i++) {
        alpha_[i] = mainDiagonal[i];
        gamma_[i] = firstDiagonal[i];
    }
    alpha_[N_-1] = mainDiagonal[N_-1];
}

QmomStd::QmomStd(int nMoments, std::shared_ptr<CoreInversionAlgorithm> coreInversion, 
    std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver)
    :
    Qmom(nMoments, numberOfNodes(nMoments), coreInversion, eigenSolver, linearSolver)
{

    if (nMoments_ % 2 != 0) {
        std::string msg = "QmomStd::Qmomstd got invalid argument:"
            " The number of moments `nMoments` must be a positive even integer";
        throw std::invalid_argument(msg);
    }
}

QmomStd::~QmomStd(){}


int QmomStd::numberOfNodes(int nMoments) const {
    return nMoments/2;
}


int QmomStd::computeQuadrature(double *moments) 
{
    // Jacobi matrix -> quadrature
    int info = eigenSolver_->compute(alpha_, gamma_);

    // The quadrature nodes equal the eigenvalues
    for (int i=0; i<N_; i++) {
        quadratureNodes_[i] = eigenSolver_->eigenValues()[i];
    }

    // Determine weights from eigenvectors if selected
    if (computeWeightsFromEigenVectors_) {
        for (int i=0; i<N_; i++) {
            quadratureWeights_[i] = eigenVecs_[i]*eigenVecs_[i];
        }
    }

    // ... or solve Vandermonde system otherwise.
    else {

        // First row of Vandermonde matrix that is not required to solve the linear system
        int end = linearSolver_->endOfElementsNeeded()/N_;

        updateVandermondeMatrix(end);

        linearSolver_->solve(vandermondeMatrix_, moments, quadratureWeights_.get());

    }

    return info;
}


QmomGaG::QmomGaG(int nMoments, std::shared_ptr<CoreInversionAlgorithm> coreInversion, 
    std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver)
    :
    Qmom(nMoments, numberOfNodes(nMoments), coreInversion, eigenSolver, linearSolver),
    work_(static_cast<double*>(mkl_malloc(N_*sizeof(double), MALLOC_ALIGN)))
{

    if (nMoments_ % 2 != 0) {
        std::string msg = "QmomGaG::QmomGaG got invalid argument:"
            " The number of moments `nMoments` must be a positive even integer";
        throw std::invalid_argument(msg);
    }
}


QmomGaG::~QmomGaG() {

    mkl_free(work_);
}


int QmomGaG::numberOfNodes(int nMoments) const {
    return nMoments - 1;
}


int QmomGaG::computeQuadrature(double *moments) 
{

    // Modify Jacobi matrix for anti-Gaussian quadrature
    gamma_[N_-2] *= constants::sqrt2;

    int info = 0;

    // Compute Gaussian quadrature in the first iteration and then the
    // anti-Gaussian quadrature using the modified Jacobi matrix
    for (int i=1; i>-1; i--) {

        // Set size of eigenvalue problem
        eigenSolver_->setSize(N_ - i);

        // Jacobi matrix -> quadrature
        info += eigenSolver_->compute(alpha_, gamma_);

        // The quadrature nodes equal the eigenvalues
        for (int j=0; j<N_-i; j++) {
            quadratureNodes_[2*j + i] = eigenSolver_->eigenValues()[j];
        }

        // Determine weights from eigenvectors if selected
        if (computeWeightsFromEigenVectors_) {
            for (int j=0; j<N_-i; j++) {
                quadratureWeights_[2*j + i] = eigenVecs_[j]*eigenVecs_[j];
            }
        }

        // ... or solve Vandermonde system otherwise.
        else {

            // First row of Vandermonde matrix that is not required to solve the linear system
            int end = linearSolver_->endOfElementsNeeded()/N_;

            updateVandermondeMatrix(end);

            linearSolver_->solve(vandermondeMatrix_, moments, work_);

        }

    }

    // Averaged quadrature
    for (int j=0; j<nNodes_; j++) {
        quadratureWeights_[j] *= 0.5;
    }

    return info;
}