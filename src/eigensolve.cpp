/**
 * @file eigensolve.cpp
 * @author M. Puetz
 * @brief Implementation of solvers for eigenvalue problems.
 * @date 2022-10-08
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "eigensolve.hpp"
#include <mkl.h>


int EigenTridiagonalSymmLapackQR::compute(double *mainDiag, double *firstDiag) {

    eigenVals_[size_-1] = mainDiag[size_-1];
    for (int i=0; i<size_-1; i++) {
        eigenVals_[i] = mainDiag[i];      // array is overwritten with eigenvalues in dstev
        firstDiag_[i] = firstDiag[i];
    }

    lapack_int ldz = static_cast<lapack_int>(size_);   // Leading dimensions, see `ldz` variable in Lapack function `dstev`
    int info = LAPACKE_dstev(LAPACK_ROW_MAJOR, static_cast<char>(problemType_), size_, eigenVals_.get(), firstDiag_, eigenVecs_.get(), ldz);

    return info;
}


int EigenTridiagonalSymmLapackDivideAndConquer::compute(double *mainDiag, double *firstDiag) {

    eigenVals_[size_-1] = mainDiag[size_-1];
    for (int i=0; i<size_-1; i++) {
        eigenVals_[i] = mainDiag[i];      // array is overwritten with eigenvalues in dstevd
        firstDiag_[i] = firstDiag[i];
    }

    lapack_int ldz = static_cast<lapack_int>(size_);   // Leading dimensions, see `ldz` variable in Lapack function `dstevd`
    int info = LAPACKE_dstevd(LAPACK_ROW_MAJOR, static_cast<char>(problemType_), size_, 
                                eigenVals_.get(), firstDiag_, eigenVecs_.get(), ldz);

    return info;
}


int EigenTridiagonalSymmLapackBisection::compute(double *mainDiag, double *firstDiag) {

    mainDiag_[size_-1] = mainDiag[size_-1];
    for (int i=0; i<size_-1; i++) {
        mainDiag_[i] = mainDiag[i];
        firstDiag_[i] = firstDiag[i];
    }

    lapack_int eigenValsFound = -1;
    lapack_int ldz = static_cast<lapack_int>(size_);   // Leading dimensions, see `ldz` variable in Lapack function `dstevx`

    int info = LAPACKE_dstevx(LAPACK_ROW_MAJOR, static_cast<char>(problemType_), static_cast<char>(eigenValSet_), size_, mainDiag_, firstDiag_, 
                                lowerBound_, upperBound_, lowerIndex_, upperIndex_, absTol_, &eigenValsFound, 
                                eigenVals_.get(), eigenVecs_.get(), ldz, lapackIntWorkspace_);

    return info;
}


int EigenTridiagonalSymmLapackRRR::compute(double *mainDiag, double *firstDiag) {

    mainDiag_[size_-1] = mainDiag[size_-1];
    for (int i=0; i<size_-1; i++) {
        mainDiag_[i] = mainDiag[i];
        firstDiag_[i] = firstDiag[i];
    }

    lapack_int eigenValsFound = size_;
    lapack_int ldz = static_cast<lapack_int>(size_);   // Leading dimensions, see `ldz` variable in Lapack function `dstevr`

    int info = LAPACKE_dstevr(LAPACK_ROW_MAJOR, static_cast<char>(problemType_), static_cast<char>(eigenValSet_), size_, mainDiag_, firstDiag_, 
                                lowerBound_, upperBound_, lowerIndex_, upperIndex_, absTol_, &eigenValsFound, 
                                eigenVals_.get(), eigenVecs_.get(), ldz, lapackIntWorkspace_);

    return info;
}


int EigenTridiagonalSymmEigenlibQR::compute(double *mainDiag, double *firstDiag)
{

    for (int i=0; i<size_-1; i++) {
        mainDiag_[i] = mainDiag[i];
        firstDiag_[i] = firstDiag[i];
    }
    mainDiag_[size_-1] = mainDiag[size_-1];

    selfAdjointEigenSolver_.computeFromTridiagonal(
        mainDiagMap_.head(size_), firstDiagMap_.head(size_ - 1), eigenProblemType_
    );

    eigenValsMap_.head(size_) = selfAdjointEigenSolver_.eigenvalues();
    
    // This condition is necessary for the program not to crash in debug mode
    if (eigenProblemType_ == Eigen::DecompositionOptions::ComputeEigenvectors) {
        for (int i=0; i<size_; i++) {
            for (int j=0; j<size_; j++) {
                eigenVecs_[i*size_ + j] = selfAdjointEigenSolver_.eigenvectors()(i,j);
            }
        }
        //eigenVecs_ = selfAdjointEigenSolver_.eigenvectors().data();
    }

    return 0;
}


RealEigenSolver::RealEigenSolver(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, 
                    double absTol, double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    maxSize_(maxSize),
    size_(maxSize),
    problemType_(problemType),
    eigenValSet_(eigenValSet),
    absTol_(absTol),
    lowerBound_(lowerBound),
    upperBound_(upperBound),
    lowerIndex_(lowerIndex),
    upperIndex_(upperIndex),
    eigenVals_(static_cast<double*>(mkl_malloc(maxSize*sizeof(double), MALLOC_ALIGN)), &mkl_free),
    eigenVecs_(static_cast<double*>(mkl_malloc(maxSize*maxSize*sizeof(double), MALLOC_ALIGN)), &mkl_free)
{
}

RealEigenSolver::~RealEigenSolver()
{
}

std::unique_ptr<RealEigenSolver> RealEigenSolver::makeUnique(
                    const std::string& typeName,
                    int size, 
                    EigenProblemType problemType, 
                    EigenValSetType eigenValSet, 
                    double absTol,
                    double lowerBound,
                    double upperBound, 
                    int lowerIndex, 
                    int upperIndex
                    )
{
    return RealEigenSolverFactory::makeUnique(typeName, size, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex);
}

std::shared_ptr<RealEigenSolver> RealEigenSolver::makeShared(
                    const std::string& typeName,
                    int size, 
                    EigenProblemType problemType, 
                    EigenValSetType eigenValSet, 
                    double absTol,
                    double lowerBound,
                    double upperBound, 
                    int lowerIndex, 
                    int upperIndex
                    )
{
    return RealEigenSolverFactory::makeShared(typeName, size, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex);
}

void RealEigenSolver::setProblemType(EigenProblemType problemType) {
    problemType_ = problemType;
}

TridiagonalSymmetricEigenSolver::TridiagonalSymmetricEigenSolver(int maxSize, EigenProblemType problemType,
                        EigenValSetType eigenValSet, double lowerBound, double upperBound, double absTol,
                        int lowerIndex, int upperIndex)
    :
    RealEigenSolver(maxSize, problemType, eigenValSet, lowerBound, upperBound, absTol, lowerIndex, upperIndex),
    mainDiag_(static_cast<double*>(mkl_malloc(maxSize*sizeof(double), MALLOC_ALIGN))),
    firstDiag_(static_cast<double*>(mkl_malloc((maxSize - 1)*sizeof(double), MALLOC_ALIGN)))
{
}

TridiagonalSymmetricEigenSolver::~TridiagonalSymmetricEigenSolver() {
    mkl_free(mainDiag_);
    mkl_free(firstDiag_);
}


EigenTridiagonalSymmLapackBisection::EigenTridiagonalSymmLapackBisection(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, double absTol, 
                    double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    TridiagonalSymmetricEigenSolver(maxSize, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex),
    lapackIntWorkspace_(static_cast<lapack_int*>(mkl_malloc(maxSize*sizeof(lapack_int), MALLOC_ALIGN)))
{
}

EigenTridiagonalSymmLapackBisection::~EigenTridiagonalSymmLapackBisection()
{
    mkl_free(lapackIntWorkspace_);
}


EigenTridiagonalSymmLapackQR::EigenTridiagonalSymmLapackQR(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, double absTol, 
                    double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    TridiagonalSymmetricEigenSolver(maxSize, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex)
{
}

EigenTridiagonalSymmLapackQR::~EigenTridiagonalSymmLapackQR(){}


EigenTridiagonalSymmLapackDivideAndConquer::EigenTridiagonalSymmLapackDivideAndConquer(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, 
                    double absTol, double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    TridiagonalSymmetricEigenSolver(maxSize, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex)
{
}

EigenTridiagonalSymmLapackDivideAndConquer::~EigenTridiagonalSymmLapackDivideAndConquer(){}


EigenTridiagonalSymmLapackRRR::EigenTridiagonalSymmLapackRRR(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, double absTol, 
                    double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    TridiagonalSymmetricEigenSolver(maxSize, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex),
    lapackIntWorkspace_(static_cast<lapack_int*>(mkl_malloc(2*maxSize*sizeof(lapack_int), MALLOC_ALIGN)))
{
}

EigenTridiagonalSymmLapackRRR::~EigenTridiagonalSymmLapackRRR()
{
    mkl_free(lapackIntWorkspace_);
}


EigenTridiagonalSymmEigenlibQR::EigenTridiagonalSymmEigenlibQR(int maxSize, EigenProblemType problemType, EigenValSetType eigenValSet, double absTol, 
                    double lowerBound, double upperBound, int lowerIndex, int upperIndex)
    :
    TridiagonalSymmetricEigenSolver(maxSize, problemType, eigenValSet, absTol, lowerBound, upperBound, lowerIndex, upperIndex),
    selfAdjointEigenSolver_(size_),
    eigenProblemType_((problemType == EigenProblemType::EigenPairs) 
                        ? Eigen::DecompositionOptions::ComputeEigenvectors 
                        : Eigen::DecompositionOptions::EigenvaluesOnly),
    mainDiagMap_(mainDiag_, size_),
    firstDiagMap_(firstDiag_, size_ - 1),
    eigenValsMap_(eigenVals_.get(), size_)
{
}

EigenTridiagonalSymmEigenlibQR::~EigenTridiagonalSymmEigenlibQR(){}

void EigenTridiagonalSymmEigenlibQR::setProblemType(EigenProblemType problemType) {

    RealEigenSolver::setProblemType(problemType);

    if (problemType_ == EigenProblemType::EigenPairs) {
        eigenProblemType_ = Eigen::DecompositionOptions::ComputeEigenvectors;
    }
    else if (problemType_ == EigenProblemType::EigenValsOnly) {
        eigenProblemType_ = Eigen::DecompositionOptions::EigenvaluesOnly;
    }
    else {
        // This should never happen
        throw std::invalid_argument("Invalid argument for eigen problemType.");
    }
}