/**
 * @file linsolve.cpp
 * @author M. Puetz
 * @brief Implementation of linear solver functions.
 * @date 2022-09-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "linsolve.hpp"
#include <mkl.h>

int solveVandermondeSystemBjorckPereyra(int n, double *V, double *b, double *x) {

    // The Vandermonde matrix is completely determined by the first row (assuming indexing beginning at zero)
    const double *v1 = &V[n];

    // START: Procedure `pvand` in the publication by Bjorck and Pereyra
    for (int i=0; i<n; i++) {
        x[i] = b[i];
    }
    for (int k=0; k<n; k++) {
        for (int j=n-1; j>k; j--) {
            x[j] -= v1[k]*x[j-1];
        }
    }
    for (int k=n-1; k>0; k--) {
        for (int j=k; j<n; j++) {
            x[j] /= v1[j] - v1[j-k];
            x[j-1] -= x[j];
        }
    }
    // END

    return 0;
}


LinearSolver::LinearSolver(unsigned int maxSize)
    :
    maxSize_(maxSize),
    size_(maxSize)
{
}


LinearSolver::~LinearSolver()
{
}


std::unique_ptr<LinearSolver> LinearSolver::makeUnique(const std::string &typeName, unsigned int maxSize)
{
    return LinearSolverFactory::makeUnique(typeName, maxSize);
}


std::shared_ptr<LinearSolver> LinearSolver::makeShared(const std::string &typeName, unsigned int maxSize)
{
    return LinearSolverFactory::makeShared(typeName, maxSize);
}


int LinearSolver::beginOfElementsNeeded() const {
    return 0;
}


int LinearSolver::endOfElementsNeeded() const {
    return size_*size_;
}


int LinearSolver::operator()(double* A, double* b, double* x) {
    return this->solve(A, b, x);
}


LinearVandermondeSolver::LinearVandermondeSolver(unsigned int maxSize)
    :
    LinearSolver(maxSize)
{
}


LinearVandermondeSolver::~LinearVandermondeSolver() 
{
}


int LinearVandermondeSolver::solve(double *V, double *b, double *x) 
{
    return solveVandermondeSystemBjorckPereyra(size_, V, b, x);
}


int LinearVandermondeSolver::beginOfElementsNeeded() const {
    return size_;
}


int LinearVandermondeSolver::endOfElementsNeeded() const {
    return 2*size_;
}


LinearLapackGesvSolver::LinearLapackGesvSolver(unsigned int maxSize)
    :
    LinearSolver(maxSize),
    work_(static_cast<double*>(mkl_malloc(maxSize_*maxSize_*sizeof(double), MALLOC_ALIGN))),
    pivotIndices_(static_cast<lapack_int*>(mkl_malloc(maxSize_*sizeof(lapack_int), MALLOC_ALIGN))) 
{
}


LinearLapackGesvSolver::~LinearLapackGesvSolver() {

    mkl_free(work_);
    mkl_free(pivotIndices_);
}


int LinearLapackGesvSolver::solve(double *A, double *b, double *x)
{
    lapack_int nRhs = 1;    // number of RHS vectors

    // Passed RHS vector is overwritten in LAPACK function
    for (unsigned int i=0; i<size_; i++) {
        x[i] = b[i];
    }
    // Copy matrix into workspace (the matrix passed to `LAPACKE_dgesv` is overwritten during the computation)
    for (unsigned int i=0; i<size_*size_; i++) {
        work_[i] = A[i];
    }

    lapack_int info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, size_, nRhs, work_, size_, pivotIndices_, x, nRhs);

    return info; 
}


LinearEigenlibPartialPivLuSolver::LinearEigenlibPartialPivLuSolver(unsigned int maxSize)
    :
    LinearSolver(maxSize),
    matrixMap_(nullptr, maxSize_, maxSize_),
    xVectorMap_(nullptr, maxSize_),
    bVectorMap_(nullptr, maxSize_),
    partialPivLuSolver_(nullptr)
{
}


LinearEigenlibPartialPivLuSolver::~LinearEigenlibPartialPivLuSolver() 
{
}


int LinearEigenlibPartialPivLuSolver::solve(double *A, double *b, double *x)
{

    // If size has changed maps must be adjusted to the new one
    if (bVectorMap_.size() != size_) {
        new (&matrixMap_) MatrixMap(A, size_, size_);
        new (&bVectorMap_) VectorMap(b, size_);
        new (&xVectorMap_) VectorMap(x, size_);
    }
    else
    {
        // Assign pointers to matrix, r.h.s vector and solution vector if
        // location has changed
        if (matrixMap_.data() != A) {
            new (&matrixMap_) MatrixMap(A, size_, size_);
        }
        if (bVectorMap_.data() != b) {
            new (&bVectorMap_) VectorMap(b, size_);
        }
        if (xVectorMap_.data() != x) {
            new (&xVectorMap_) VectorMap(x, size_);
        }
    }

    if (!partialPivLuSolver_) {
        partialPivLuSolver_.reset(
            new Eigen::PartialPivLU<Matrix>(matrixMap_) // LU decomposition is computed on construction
        );
    }
    else {
        partialPivLuSolver_->compute(matrixMap_);
    }

    xVectorMap_ = partialPivLuSolver_->solve(bVectorMap_);

    return 0;
}