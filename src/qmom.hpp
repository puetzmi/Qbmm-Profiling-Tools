/**
 * @file qmom.hpp
 * @author Michele Puetz
 * @brief Classes implementing the quadrature method of moments and derived methods.
 * @date 2022-05-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef QMOM_HPP
#define QMOM_HPP

#include <mkl.h>
#include <functional>

#include "global_defs.hpp"
#include "eigensolve.hpp"
#include "linsolve.hpp"
#include "core_inversion.hpp"
#include "factory.hpp"


class Qmom;     // forward declaration

/// Factory for `Qmom` objects
using QmomFactory = 
    Factory<Qmom,
        int,    // nMoments
        std::shared_ptr<CoreInversionAlgorithm>,    //coreInversion
        std::shared_ptr<RealEigenSolver>,   // eigenSolver
        std::shared_ptr<LinearSolver>,  // linearSolver,
        std::function<
        int(double * const, double * const, int, int, double*)> //momentsRateOfChangeFunction
    >;

/**
 * @brief Class implementing methods derived from the quadrature method of moments @cite mcgraw1997.
 * 
 * The essential procedure performed within this class is to 
 * 
 * -# take a set of @f$ 2N @f$ integer moments defined as
 * @f\[
 *  m_k = \int_{\Omega} n(\xi) \xi^k \mathrm{d} \xi, \quad k=0,1,\dots,2N-1,
 * @f\]
 * where @f$ n(\xi) @f$ is a number density function with support @f$ \Omega @f$,
 * 
 * -# compute from the given moments the Jacobi matrix
 * @f\[
 * \mathbf{J} = 
 *  \begin{pmatrix}
 *  \alpha_0 & \gamma_1 &  \\
 *  \gamma_1 & \alpha_1 & \ddots  \\
 *   & \ddots & \ddots & \gamma_{N-1} \\
 *   & & \gamma_{N-1} & \alpha_{N-1}
 *  \end{pmatrix},
 * @f\]
 * where @f$ \gamma_k =  \sqrt{\beta_k} @f$ and @f$ \alpha_k @f$ and @f$ \beta_k @f$ are the recurrence coefficients
 * of the orthogonal polynomials with respect to @f$ n(\xi) @f$ (see e.g. @cite Gautschi2004),
 * 
 * -# compute the quadrature nodes @f$ \xi_j @f$ and weights @f$ w_j, \; j=1,\dots,M @f$ from the Jacobi matrix, where
 * @f$ M @f$ depends on the specific method,
 * 
 * -# use the obtained quadrature to approximate an integral expression of the form
 * @f\[
 *  \int_{\Omega} n(\xi) g(\xi) \mathrm{d} \xi \approx \sum_{j} w_j g(\xi_j),
 * @f\]
 * where @f$ g @f$ is an arbitrary function.
 */
class Qmom {

protected:

    int nMoments_;                              ///< Number of moments.
    int N_;                                     ///< N as defined above.
    int nNodes_;                                ///< Number of quadrature nodes (depends on specific QMOM type).

    const std::shared_ptr<CoreInversionAlgorithm> 
        coreInversion_;                         ///< Core inversion algorithm to compute the Jacobi matrix from moments.

    const std::shared_ptr<RealEigenSolver> 
        eigenSolver_;                           ///< Eigenvalue algorithm.

    const std::shared_ptr<LinearSolver> 
        linearSolver_;                          ///< Linear solver (may be null-pointer if weights are to be computed from eigenvectors of the Jacobi matrix).

    const bool computeWeightsFromEigenVectors_; ///< If true eigenvectors are computed to obtain weights, linear Vandermonde system is solved otherwise.

    std::shared_ptr<double[]> quadratureNodes_; ///< Quadrature nodes. 
    std::shared_ptr<double[]> 
        quadratureWeights_;                     ///< Quadrature weights.

    double *vandermondeMatrix_;                 ///< Vectorized Vandermonde matrix (C-style)
    double *alpha_;                             ///< Diagonal elements of Jacobi matrix.
    double *gamma_;                             ///< Sub-/superdiagonal elements of Jacobi matrix.
    std::shared_ptr<const double[]> eigenVecs_; ///< Pointer to eigenvectors in Eigen solver.
    std::shared_ptr<double[]> 
        momentsRateOfChange_;                   ///< Storage for the rate of change of 2N moments; updated in `computeMomentsRateOfChange`.
    
    /**
     * @brief Function or callable object computing the moments' rate of change given
     * the abscissas, weights and size of the quadrature, number of moments and an array
     * with allocated memory to store the result (in this order).
     * 
     */
    std::function<int(double * const ,double * const, int, int, 
        double*)> momentsRateOfChangeFunction_;


    /**
     * @brief Compute number of quadrature nodes from number of moments.
     * 
     * @param nMoments Number of moments.
     * @return int Number of quadrature nodes.
     */
    virtual int numberOfNodes(int nMoments) const = 0;

    
    /**
     * @brief Update selected rows of the Vandermonde matrix with current quadrature nodes.
     * 
     * The Vandermonde matrix is here defined in the transposed form, i.e.
     * @f\[
     *  \mathbf{V}^{\mathrm{T}} = 
     * \begin{pmatrix}
     * 1 & 1 & \dots & 1 \\
     * \xi_1 & \xi_2 & \dots & \xi_{N} \\
     * \xi_1^2 & \xi_2^2 & \dots & \xi_N^2 \\ 
     * \vdots & \vdots & \dots & \vdots \\ 
     * \xi_1^{2N-1} & \xi_2^{2N-1} & \dots & \xi_N^{2N-1}
     * \end{pmatrix}.
     * @f\]
     * 
     * @param end The row after the last to be updated, i.e. rows up to end-1 are updated.
     * @param begin First row to be updated, must satisfy `1 <= begin`, by default `begin=1` since the zeroth row contains only ones.
     * @param nodes Nodes to update Vandermonde matrix, set to the currently stored `quadratureNodes` is not given.
     */
    void updateVandermondeMatrix(int end, int begin=1, double *nodes=nullptr);


public:

    /**
     * @brief Compute reconstructed moments @f$ m_k @f$ where the range of @f$ k @f$ is provided by parameters.
     * 
     * The @f$ k @f$-th reconstructed moment is 
     * @f\[
     *  m_k = \int_{\Omega} \xi^k n(\xi) \mathrm{d} \xi,
     * @f\]
     * No moment inversion is performed here, the currently stored quadrature nodes and weights are used for the calculation.
     * @param [in] kmin Lower bound of @f$ k @f$.
     * @param [in] kmax Exclusive upper bound of @f$ k @f$ that must satisfy `begin <= end <= 2*nNodes_`.
     * @param [out] moments Reconstructed moments. Allocated memory must be at least of size `kmax - kmin`.
     * @return int Dummy error flag.
     */
    virtual int computeMoments(int kmin, int kmax, double *moments);

    /**
     * @brief Construct a new `Qmom` object.
     * 
     * @param nMoments  Number of moments.
     * @param int Number of quadrature nodes (depends on specific QMOM-type).
     * @param coreInversion Shared pointer to `CoreInversionAlgorithm` object to compute the Jacobi matrix from moments.
     * @param eigenSolver Shared pointer to `EigenSolver` object to compute the eigenvalues (and possibly -vectors) of the tridiagonal Jacobi matrix.
     * @param linearSolver Shared pointer to `LinearSolver` to compute quadrature weights from Vandermonde system, or null pointer if the weights are to be computed from the eigenvectors of the Jacobi matrix.
     * @param momentsRateOfChangeFunction Function to compute moments' rate of change, see corresponing attribute `Qmom::momantsRateOfChangeFunction_`.
     */
    Qmom(int nMoments, int nNodes, std::shared_ptr<CoreInversionAlgorithm> coreInversion, 
        std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver,
        std::function<int(double * const, double * const, int, int, double*)>momentsRateOfChangeFunction);

    
    /**
     * @brief Destroy the Qmom object.
     * 
     */
    virtual ~Qmom();


    /**
     * @brief Create new `Qmom` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `QmomFactory`.
     * @param nMoments  Number of moments.
     * @param coreInversion Shared pointer to `CoreInversionAlgorithm` object to compute the Jacobi matrix from moments.
     * @param eigenSolver Shared pointer to `EigenSolver` object to compute the eigenvalues (and possibly -vectors) of the tridiagonal Jacobi matrix.
     * @param linearSolver Shared pointer to `LinearSolver` to compute quadrature weights from Vandermonde system, or null pointer if the weights are to be computed from the eigenvectors of the Jacobi matrix.
     * @param momentsRateOfChangeFunction Function to compute moments' rate of change, see corresponing attribute `Qmom::momantsRateOfChangeFunction_`.
     * @return std::shared_ptr<Qmom> Shared pointer to new `Qmom` object.
     */
    static std::shared_ptr<Qmom> makeShared(const std::string &typeName, int nMoments, 
        std::shared_ptr<CoreInversionAlgorithm> coreInversion, std::shared_ptr<RealEigenSolver> eigenSolver,
        std::shared_ptr<LinearSolver> linearSolver,
        std::function<int(double * const, double * const, int, int, double*)>momentsRateOfChangeFunction);


    /**
     * @brief Create new `Qmom` object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `QmomFactory`.
     * @param nMoments  Number of moments.
     * @param coreInversion Shared pointer to `CoreInversionAlgorithm` object to compute the Jacobi matrix from moments.
     * @param eigenSolver Shared pointer to `EigenSolver` object to compute the eigenvalues (and possibly -vectors) of the tridiagonal Jacobi matrix.
     * @param linearSolver Shared pointer to `LinearSolver` to compute quadrature weights from Vandermonde system, or null pointer if the weights are to be computed from the eigenvectors of the Jacobi matrix.
     * @param momentsRateOfChangeFunction Function to compute moments' rate of change, see corresponing attribute `Qmom::momantsRateOfChangeFunction_`.
     * @return std::unique_ptr<Qmom> Unique pointer to new `Qmom` object.
     */
    static std::unique_ptr<Qmom> makeUnique(const std::string &typeName, int nMoments, 
        std::shared_ptr<CoreInversionAlgorithm> coreInversion, std::shared_ptr<RealEigenSolver> eigenSolver,
        std::shared_ptr<LinearSolver> linearSolver,
        std::function<int(double * const, double * const, int, int, double*)>momentsRateOfChangeFunction);


    /**
     * @brief Essential function performing all steps stated in the class description.
     * 
     * @param m Set of moments.
     * @return int Error flag.
     */
    virtual int compute(double *moments);

    /**
     * @brief Compute the tridiagonal Jacobi matrix given a moment sequence.
     * 
     * @param moments A realizable moment sequence.
     * @return int Error flag.
     */
    virtual int computeJacobiMatrix(double *moments);

    /**
     * @brief Compute quadrature nodes and weights presuming an already updated Jacobi matrix.
     * 
     * @param moments A realizable moment sequence.
     * @return int Error flag
     */
    virtual int computeQuadrature(double *moments) = 0;


    /**
     * @brief Compute @f$ \dot{m}_k = \mathrm{d}m_k/\mathrm{d}t @f$ where the range of @f$ k @f$ is provided by parameters.
     * 
     * The @f$ k @f$-th moment source term is 
     * @f\[
     *  \dot{m}_k = \int_{\Omega} \xi^k g(\xi) \mathrm{d} \xi,
     * @f\]
     * where @f$ g @f$ is a problem-dependent function considered in `momentsRateOfChangeFunction`.
     * No moment inversion is performed here, the current quadrature nodes Qmom::xi_ and Qmom::w_ are used for the calculation.
     * The result is stored in Qmom::mDot_. The moment orders @f$ k @f$ range from `begin` to `end - 1`.
     * @param nMoments Number of moments to be computed.
     * @return int Dummy error flag.
     */
    virtual int computeMomentsRateOfChange(int nMoments);

    /**
     * @brief Get shared pointer to quadrature abscissas.
     * 
     * @return std::shared_ptr<const double[]> Pointer to quadrature abscissas (read-only).
     */
    inline std::shared_ptr<const double[]> quadratureNodes() const {
        return quadratureNodes_;
    }

    /**
     * @brief Get shared pointer to quadrature weights.
     * 
     * @return std::shared_ptr<const double[]> Pointer to quadrature weights (read-only).
     */
    inline std::shared_ptr<const double[]> quadratureWeights() const {
        return quadratureWeights_;
    }

    /**
     * @brief Get shared pointer to moments' rate of change.
     * 
     * @return std::shared_ptr<const double[]> Pointer to moments' rate of change (read-only).
     */
    inline std::shared_ptr<const double[]> momentsRateOfChange() const {
        return momentsRateOfChange_;
    }

    /**
     * @brief Get number of moments.
     * 
     * @return int Number of moments
     */
    inline int numberOfMoments() const {
        return nMoments_;
    }


    /**
     * @brief Set the Jacobi matrix "manually" (e.g. to use exact recurrence coefficients).
     * 
     * @param mainDiagonal Main diagonal elements.
     * @param firstDiagonal First sub-/superdiagonal elements.
     * 
     */
    virtual void setJacobiMatrix(double *mainDiagonal, double *firstDiagonal);

};


/**
 * @brief Class implementing the standard quadrature method of moments (QMOM) @cite mcgraw1997.
 * 
 */
class QmomStd : public Qmom {

public:

    /**
     * @brief Construct a new `QmomStd` object.
     * 
     * @param nMoments  Number of moments.
     * @param coreInversion Shared pointer to `CoreInversionAlgorithm` object to compute the Jacobi matrix from moments.
     * @param eigenSolver Shared pointer to `EigenSolver` object to compute the eigenvalues (and possibly -vectors) of the tridiagonal Jacobi matrix.
     * @param linearSolver Shared pointer to `LinearSolver` to compute quadrature weights from Vandermonde system, or null pointer if the weights are to be computed from the eigenvectors of the Jacobi matrix.
     * @param momentsRateOfChangeFunction Function to compute moments' rate of change, see corresponing attribute `Qmom::momantsRateOfChangeFunction_`.
     */
    QmomStd(int nMoments, std::shared_ptr<CoreInversionAlgorithm> coreInversion, 
        std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver,
        std::function<int(double * const, double * const, int, int, double*)> momentsRateOfChangeFunction);
    
    /**
     * @brief Destroy the `QmomStd` object.
     * 
     */
    virtual ~QmomStd();

    virtual int numberOfNodes(int) const;
    virtual int computeQuadrature(double *moments);

private:

    static bool inline registered_ = REGISTER_TYPE(QmomFactory, QmomStd);

};


/**
 * @brief Class implementing the Gauss/anti-Gauss quadrature method of moments @cite Puetz2022.
 * 
 */
class QmomGaG : public Qmom {

private:

    double *work_;              //< Additional workspace of size N_ to store weights.

public:

    /**
     * @brief Construct a new `QmomGaG` object.
     * 
     * @param nMoments  Number of moments.
     * @param coreInversion Shared pointer to `CoreInversionAlgorithm` object to compute the Jacobi matrix from moments.
     * @param eigenSolver Shared pointer to `EigenSolver` object to compute the eigenvalues (and possibly -vectors) of the tridiagonal Jacobi matrix.
     * @param linearSolver Shared pointer to `LinearSolver` to compute quadrature weights from Vandermonde system, or null pointer if the weights are to be computed from the eigenvectors of the Jacobi matrix.
     * @param momentsRateOfChangeFunction Function to compute moments' rate of change, see corresponing attribute `Qmom::momantsRateOfChangeFunction_`.
     */
    QmomGaG(int nMoments, std::shared_ptr<CoreInversionAlgorithm> coreInversion, 
        std::shared_ptr<RealEigenSolver> eigenSolver, std::shared_ptr<LinearSolver> linearSolver,
        std::function<int(double * const, double * const, int, int, double*)> momentsRateOfChangeFunction);

    /**
     * @brief Destroy the `QmomGaG` object.
     * 
     */
    virtual ~QmomGaG();

    virtual int numberOfNodes(int) const;
    virtual int computeQuadrature(double *moments);


private:

    static bool inline registered_ = REGISTER_TYPE(QmomFactory, QmomGaG);

};

#endif // QMOM_HPP