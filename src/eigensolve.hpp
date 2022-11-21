/**
 * @file eigensolve.hpp
 * @author M. Puetz
 * @brief Solvers for eigenvalue problems, particularly those associated with real symmetric tridiagonal matrices.
 * @date 2022-05-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef EIGENSOLVE_HPP
#define EIGENSOLVE_HPP

#include "global_defs.hpp"
#include "eigen_defs.hpp"
#include <limits>
#include <mkl.h>
#include <memory>
#include <map>
#include <functional>
#include <vector>
#include "factory.hpp"
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/Eigenvalues>


/// @brief Problem type in terms of target quantities.
enum class EigenProblemType : char {
    EigenValsOnly = 'N',   ///< Compute eigenvalues only
    EigenPairs    = 'V',   ///< Compute eigenvalues and eigenvectors
};


/// @brief Specification of which eigenvalues shall be computed.
enum class EigenValSetType : char {
    All      = 'A',   ///< Compute all eigenvalues (and optionally eigenvectors)
    Interval = 'V',   ///< Compute eigenvalues (and optionally eigenvectors) within a given interval
    Indices  = 'I'    ///< Compute eigenvalues (and optionally eigenvectors) with given indices
};

class RealEigenSolver;      // forward declaration


/// @brief Factory for `RealEigenSolver` objects
using RealEigenSolverFactory = Factory< RealEigenSolver, 
                                    int,
                                    EigenProblemType,
                                    EigenValSetType,
                                    double,
                                    double,
                                    double,
                                    int,
                                    int >;

/// @brief Default parameters used in several constructors.
namespace EigenDefaultParams {
    constexpr EigenValSetType eigenValSet = EigenValSetType::All;

    // This should correspond to `2*LAPACKE_dlamch('S')` as suggested in the MLK/LAPACKE documentation for max. accuracy
    constexpr double absTol = 2*std::numeric_limits<double>::min();

    constexpr double lowerBound = std::numeric_limits<double>::max();
    constexpr double upperBound = std::numeric_limits<double>::min();
    constexpr int lowerIndex = std::numeric_limits<int>::max();
    constexpr int upperIndex = std::numeric_limits<int>::min();
}


/**
 * @brief Base class for all eigenvalue problem solvers assuming real eigenvalues and eigenvectors.
 * 
 * Derived classes provide interfaces to different algorithms to solve the eigenvalue problem
 * @f\[ 
 *      \mathbf{A} \mathbf{v}_j = \lambda_j \mathbf{v}_j,
 * @f\]
 * where @f$ \mathbf{A} @f$ is a real square matrix, and @f$ \lambda_j @f$ and @f$ \mathbf{v}_j @f$
 * are the @f$ j @f$-th eigenvalue and eigenvector, respectively.
 * 
 */
class RealEigenSolver {

public:

    /**
     * @brief Construct a new `RealEigenSolver` object
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     */
   RealEigenSolver(
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
                );

    /**
     * @brief Destroy the `RealEigenSolver` object.
     * 
     */
    virtual ~RealEigenSolver();


    /**
     * @brief Create new eigen solver object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `RealEigenSolverFactory`.
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @return std::unique_ptr<RealEigenSolver> Unique pointer to new `RealEigenSolver` object.
     */
    static std::unique_ptr<RealEigenSolver> makeUnique(
                        const std::string& typeName,
                        int maxSize, 
                        EigenProblemType problemType, 
                        EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                        double absTol=EigenDefaultParams::absTol,
                        double lowerBound=EigenDefaultParams::lowerBound,
                        double upperBound=EigenDefaultParams::upperBound, 
                        int lowerIndex=EigenDefaultParams::lowerIndex, 
                        int upperIndex=EigenDefaultParams::upperIndex
                        );

    /**
     * @brief Create new eigen solver object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `RealEigenSolverFactory`.
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @return std::shared_ptr<RealEigenSolver> Shared pointer to new `RealEigenSolver` object.
     */
    static std::shared_ptr<RealEigenSolver> makeShared(
                        const std::string& typeName,
                        int maxSize, 
                        EigenProblemType problemType, 
                        EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                        double absTol=EigenDefaultParams::absTol,
                        double lowerBound=EigenDefaultParams::lowerBound,
                        double upperBound=EigenDefaultParams::upperBound, 
                        int lowerIndex=EigenDefaultParams::lowerIndex, 
                        int upperIndex=EigenDefaultParams::upperIndex
                        );

    /**
     * @brief Get pointer to eigenvalues.
     * 
     * @return std::shared_ptr<const double[]> Shared pointer to eigenvalues.
     */
    inline std::shared_ptr<const double[]> eigenValues() {
        return eigenVals_;
    }

    /**
     * @brief Get pointer to eigenvectors.
     * 
     * @return std::shared_ptr<const double[]> Shared pointer to eigenvectors.
     */
    inline std::shared_ptr<const double[]> eigenVectors() {
        return eigenVecs_;
    }
    

    /**
     * @brief Compute eigenvalues/-vectors.
     * 
     * @return int Error flag.
     */
    virtual int compute(double*, double*) = 0;


    /**
     * @brief Set the size of input matrices.
     * 
     * @param size New size of input matrices, must be less than or equal to `maxSize`.
     */
    inline void setSize(int size) { // TODO: This needs to be tested
        
        // check validity
        if (size > maxSize_) {
            std::string msg = "RealEigenSolver::setSize got invalid parameter:"
                "The parameter `size` must not be greater than `maxSize_`";
            throw std::invalid_argument(msg);
        }

        size_ = size;
    }


    /**
     * @brief Set the `problemType` indicating whether or not eigenvectors are to be computed.
     * 
     * @param problemType The problem type, see enum class `EigenProblemType`.
     */
    virtual void setProblemType(EigenProblemType problemType);


protected:
    const int maxSize_;                     ///< Maximum matrix size allowed based on allocated memory.
    int size_;                              ///< Matrix size.
    EigenProblemType problemType_;          ///< Problem type in terms of target quantities.
    EigenValSetType eigenValSet_;           ///< Specification of eigenValSet of eigenvalues to be computed, see enum EigenValSetType.
    double absTol_;                         ///< Absolute tolerance (needed for some iterative algorithms).
    double lowerBound_;                     ///< Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
    double upperBound_;                     ///< Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
    int lowerIndex_;                        ///< Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
    int upperIndex_;                        ///< Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
    std::shared_ptr<double[]> eigenVals_;   ///< Storage for eigenvalues.
    std::shared_ptr<double[]> eigenVecs_;   ///< Storage for eigenvectors.

};


/**
 * @brief Base class for solvers of the symmetric tridiagonal eigenvalue problem.
 * 
 * Derived classes provide interfaces to different algorithms to solve the eigenvalue problem
 * @f\[ 
 *      \mathbf{T} \mathbf{v}_j = \lambda_j \mathbf{v}_j,
 * @f\]
 * where @f$ \mathbf{T} @f$ is a real symmetric tridiagonal matrix, and @f$ \lambda_j @f$ and @f$ \mathbf{v}_j @f$
 * are the @f$ j @f$-th eigenvalue and eigenvector, respectively.
 * 
 */
class TridiagonalSymmetricEigenSolver : public RealEigenSolver
{

public:
    /**
     * @brief Construct a new `TridiagonalSymmetricEigenSolver` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    TridiagonalSymmetricEigenSolver(
                        int maxSize, 
                        EigenProblemType problemType, 
                        EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                        double absTol=EigenDefaultParams::absTol,
                        double lowerBound=EigenDefaultParams::lowerBound,
                        double upperBound=EigenDefaultParams::upperBound, 
                        int lowerIndex=EigenDefaultParams::lowerIndex, 
                        int upperIndex=EigenDefaultParams::upperIndex
                        );

    /**
     * @brief Destroy the `TridiagonalSymmetricEigenSolver` object.
     * 
     */
    virtual ~TridiagonalSymmetricEigenSolver();


    /**
     * @brief Compute eigenvalues and optionally eigenvectors of real symmetric tridiagonal matrix.
     * 
     * @param mainDiag Main diagonal.
     * @param firstDiag First upper/lower diagonal.
     * @return int Error flag.
     */
    virtual int compute(double *mainDiag, double *firstDiag) = 0;


protected:

    double *mainDiag_;                      ///< Storage for main diagonal (necessary since it is overwritten in LAPACK functions).
    double *firstDiag_;                     ///< Storage for first upper/lower diagonal (necessary since it is overwritten in LAPACK functions).

};


/**
 * @brief Class serving as interface to the Intel MKL / LAPACK routine `DSTEVR`, computing the eigenvalues
 * and optionally eigenvectors of a real tridiagonal symmetric matrix using the Relatively-Robust-Representations
 * algorithm @cite Dhillon1997.
 * 
 */
class EigenTridiagonalSymmLapackRRR : public TridiagonalSymmetricEigenSolver 
{
public:
    /**
     * @brief Construct a new `EigenTridiagonalSymmLapackRRR` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    EigenTridiagonalSymmLapackRRR (
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
            );

    /**
     * @brief Destroy the `EigenTridiagonalSymmLapackRRR` object.
     * 
     */
    virtual ~EigenTridiagonalSymmLapackRRR();

    virtual int compute(double*, double*);


private:

    lapack_int *lapackIntWorkspace_;        ///< Lapack integer workspace corresponding to th LAPACK variable `isuppz` in function `dstevr`.

    static bool inline registered_ = 
        REGISTER_TYPE(RealEigenSolverFactory, EigenTridiagonalSymmLapackRRR);

};


/**
 * @brief Class serving as interface to the Intel MKL / LAPACK routine `DSTEV`, computing the eigenvalues
 * and optionally eigenvectors of a real tridiagonal symmetric matrix using the QR-algorithm.
 * 
 */
class EigenTridiagonalSymmLapackQR : public TridiagonalSymmetricEigenSolver
{
public:
    /**
     * @brief Construct a new `EigenTridiagonalSymmLapackQR` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    EigenTridiagonalSymmLapackQR (
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
            );

    /**
     * @brief Destroy the `EigenTridiagonalSymmLapackQR` object.
     * 
     */
    virtual ~EigenTridiagonalSymmLapackQR();

    virtual int compute(double*, double*);


private:

    static inline bool registered_ = REGISTER_TYPE(RealEigenSolverFactory, EigenTridiagonalSymmLapackQR);

};


/**
 * @brief Class serving as interface to the Intel MKL / LAPACK routine `DSTEVD`, computing the eigenvalues
 * and optionally eigenvectors of a real tridiagonal symmetric matrix using the divide-and-conquer-algorithm.
 * 
 */
class EigenTridiagonalSymmLapackDivideAndConquer : public TridiagonalSymmetricEigenSolver
{
public:
    /**
     * @brief Construct a new `EigenTridiagonalSymmLapackDivideAndConquer` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    EigenTridiagonalSymmLapackDivideAndConquer (
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
            );

    /**
     * @brief Destroy the `EigenTridiagonalSymmLapackDivideAndConquer` object.
     * 
     */
    virtual ~EigenTridiagonalSymmLapackDivideAndConquer();

    virtual int compute(double*, double*);


private:

    static inline bool registered_ = 
        REGISTER_TYPE(RealEigenSolverFactory, EigenTridiagonalSymmLapackDivideAndConquer);

};


/**
 * @brief Class serving as interface to the Intel MKL / LAPACK routine `DSTEVDX`, computing the eigenvalues
 * and optionally eigenvectors of a real tridiagonal symmetric matrix using the bisection-algorithm 
 * with inverse iteration (if eigenvectors are required).
 * 
 */
class EigenTridiagonalSymmLapackBisection : public TridiagonalSymmetricEigenSolver
{
public:
    /**
     * @brief Construct a new `EigenTridiagonalSymmLapackBisection` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    EigenTridiagonalSymmLapackBisection (
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
            );

    /**
     * @brief Destroy the EigenTridiagonalSymmLapackBisection object
     * 
     */
    virtual ~EigenTridiagonalSymmLapackBisection();

    virtual int compute(double*, double*);


private:

    lapack_int *lapackIntWorkspace_;        ///< Lapack integer workspace corresponding to th LAPACK variable `ifail` in function `dstevx`.

    static inline bool registered_ = 
        REGISTER_TYPE(RealEigenSolverFactory, EigenTridiagonalSymmLapackBisection);

};


/**
 * @brief Class for the solution of eigenvalue problems using the Eigen3 library.
 * 
 */
class EigenTridiagonalSymmEigenlibQR : public TridiagonalSymmetricEigenSolver
{

public:
    /**
     * @brief Construct a new `EigenTridiagonalSymmEigenlibQR` object.
     * 
     * @param maxSize Maximum matrix size.
     * @param problemType Problem type in terms of target quantities, see enum EigenProblemType.
     * @param eigenValSet Specification of which eigenvalues shall be computed, see enum EigenValSetType.
     * @param absTol The absolute error tolerance to which each eigenvalue is required (only needed by certain algorithms).
     * @param lowerBound Lower bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param upperBound Upper bound of the interval to be searched for eigenvalues (needed if eigenValSet==EigenValSetType::Interval).
     * @param lowerIndex Index of the smallest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * @param upperIndex Index of the largest eigenvalue to be returned (needed if eigenValSet==EigenValSetType::Indices).
     * 
     */
    EigenTridiagonalSymmEigenlibQR (
                int maxSize, 
                EigenProblemType problemType, 
                EigenValSetType eigenValSet=EigenDefaultParams::eigenValSet, 
                double absTol=EigenDefaultParams::absTol,
                double lowerBound=EigenDefaultParams::lowerBound,
                double upperBound=EigenDefaultParams::upperBound, 
                int lowerIndex=EigenDefaultParams::lowerIndex, 
                int upperIndex=EigenDefaultParams::upperIndex
            );

    /**
     * @brief Destroy the `EigenTridiagonalSymmEigenlibQR` object.
     * 
     */
    virtual ~EigenTridiagonalSymmEigenlibQR();

    virtual int compute(double*, double*);

    virtual void setProblemType(EigenProblemType problemType);

private:

    Eigen::SelfAdjointEigenSolver<Matrix> selfAdjointEigenSolver_;      ///< Symmetric eigenvalue solver in Eigen library.
    int eigenProblemType_;          ///< Eigen's decomposition type corresponding to EigenProblemType
    VectorMap mainDiagMap_;         ///< Eigen::Map mapping Vector to raw pointer containing main diagonal data and vice versa.
    VectorMap firstDiagMap_;        ///< Eigen::Map mapping Vector to raw pointer containing sub-/superdiagonal data and vice versa.
    VectorMap eigenValsMap_;        ///< Eigen::Map mapping eigenvalues from Vector to raw pointer and vice versa
    MatrixMap eigenVecsMap_;        ///< Eigen::Map mapping eigenvectors from Matrix to raw pointer and vice versa

    static inline bool registered_ = REGISTER_TYPE(RealEigenSolverFactory, EigenTridiagonalSymmEigenlibQR);

};

#endif  // EIGENSOLVE_HPP