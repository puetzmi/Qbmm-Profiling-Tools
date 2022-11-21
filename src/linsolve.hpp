/**
 * @file linsolve.hpp
 * @author Michele Puetz
 * @brief Classes to solve linear systems, in particular Vandermonde systems.
 * @date 2022-05-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef LINSOLVE_HPP
#define LINSOLVE_HPP

#include "global_defs.hpp"
#include <mkl.h>
#include <string>
#include "factory.hpp"
#include "eigen_defs.hpp"
#include <eigen3/Eigen/LU>


class LinearSolver;         // forward declaration

/// Factory for `LinearSolver` objects
using LinearSolverFactory = Factory<LinearSolver, int>;


/**
 * @brief Solve Vandermonde system using the efficient algorithm due to @cite Bjorck1970.
 * 
 * @param [in] n Number of variables.
 * @param [in] V Vandermonde matrix. Only the first row (with power 1) is used here.
 * @param [in] b Right-hand side vector.
 * @param [out] x Solution of the Vandermonde system.
 * @return int Error flag.
 */
int solveVandermondeSystemBjorckPereyra(int n, double *V, double *b, double *x);


/**
 * @brief Base class of linear solvers used to compute weights by solution of Vandermonde system.
 * 
 */
class LinearSolver 
{
public:
    /**
     * @brief Construct a new LinearSolver object.
     * 
     * @param size Size of the linear system.
     */
    LinearSolver(unsigned int size);

    /**
     * @brief Destroy the LinearSolver object.
     * 
     */
    virtual ~LinearSolver();

     
    /**
     * @brief Create new `LinearSolver` object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `LinearSolverFactory`.
     * @param size Size of the linear system.
     * @return std::unique_ptr<LinearSolver> Unique pointer to new `LinearSolver` object.
     */
    static std::unique_ptr<LinearSolver> makeUnique(const std::string &typeName, unsigned int size);


    /**
     * @brief Create new `LinearSolver` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the map in `LinearSolverFactory`.
     * @param size Size of the linear system.
     * @return std::shared_ptr<LinearSolver> shared pointer to new `LinearSolver` object.
     */
    static std::shared_ptr<LinearSolver> makeShared(const std::string &typeName, unsigned int size);


    /**
     * @brief Solve linear system @f$ \mathbf{A} \mathbf{x} = \mathbf{b} @f$.
     * 
     * @param [in] A Vectorized coefficient matrix (C-style).
     * @param [in] b Right-hand-side vector.
     * @param [out] x Solution of the linear system.
     * @return int Error flag.
     */
    int operator()(double *A, double *b, double *x);


    /**
     * @brief Solve linear system @f$ \mathbf{A} \mathbf{x} = \mathbf{b} @f$.
     * 
     * @param [in] A Vectorized coefficient matrix (C-style).
     * @param [in] b Right-hand-side vector.
     * @param [out] x Solution of the linear system.
     * @return int Error flag.
     */
    virtual int solve(double *A, double *b, double *x) = 0;


    /**
     * @brief Returns the begin of the required elements of the given matrix, which is an array of length `size*size`.
     * 
     * Example: If the solver needed only the second row of the given matrix A to solve the linear system
     * (due to a presumed a particular form of the matrix), the return value of this method would be 2*`size_`.
     * 
     * @return int Index of the first required element in the given array of size `size_*size_`.
     */
    virtual int beginOfElementsNeeded() const;


    /**
     * @brief Returns the end of the required elements of the given matrix (array of size n*n), 
     * i.e. the index first element that is not needed.
     * 
     * Example: If the solver needed only the first three rows of the given matrix A to solve the linear system
     * (due to a presumed a particular form of the matrix), the return value of this method would be `4*size_`.
     * 
     * @return int Index of the first not required element in the given array of size `size_*size_`.
     */
    virtual int endOfElementsNeeded() const;


protected:

    const unsigned int size_;       ///< Size of the linear system.

};


/**
 * @brief Class implementing algorithm due to @cite Bjorck1970 to solve Vandermonde system.
 * 
 */
class LinearVandermondeSolver : public LinearSolver {

public:

    /**
     * @brief Construct a new LinearVandermondeSolver object.
     * 
     * @param size Size of the linear system.
     */
    LinearVandermondeSolver(unsigned int size);

    /**
     * @brief Destroy the LinearVandermondeSolver object.
     * 
     */
    virtual ~LinearVandermondeSolver();


    /**
     * @brief Solve linear Vandermonde system @f$ \mathbf{V} \mathbf{x} = \mathbf{b} @f$.
     * 
     * @param [in] V Vectorized Vandermonde matrix (C-style, only the first row with power 1 is used).
     * @param [in] b Right-hand-side vector.
     * @param [out] x Solution of the linear system.
     * @return int Error flag.
     */
    virtual int solve(double *V, double *b, double *x);

    virtual int beginOfElementsNeeded() const;

    virtual int endOfElementsNeeded() const;


private:

    static bool inline registered_ = REGISTER_TYPE(LinearSolverFactory, LinearVandermondeSolver);

};


/**
 * @brief Class serving as interface to LAPACK GESV linear solver.
 * 
 */
class LinearLapackGesvSolver : public LinearSolver {

public:

    /**
     * @brief Construct a new LinearLapackGesvSolver object.
     * 
     * @param size Size of the linear system.
     */
    LinearLapackGesvSolver(unsigned int size);

    /**
     * @brief Destroy the LinearLapackGesvSolver object.
     * 
     */
    virtual ~LinearLapackGesvSolver();


    /**
     * @brief Solve linear system @f$ \mathbf{A} \mathbf{x} = \mathbf{b} @f$.
     * 
     * @param [in] A Vectorized coefficient matrix (C-style).
     * @param [in] b Right-hand-side vector.
     * @param [out] x Solution of the linear system.
     * @return int Error flag.
     */
    virtual int solve(double *A, double *b, double *x);


private:

    double *work_;                 ///< Workspace for the matrix, which gets overwritten in LAPACK routines
    lapack_int *pivotIndices_;     ///< Workspace for pivot element indices needed by LAPACK routines.

    static bool inline registered_ = REGISTER_TYPE(LinearSolverFactory, LinearLapackGesvSolver);

};


/**
 * @brief Class serving as interface to the linear solver in the Eige library using a partial pivoting LU decomposition.
 * 
 */
class LinearEigenlibPartialPivLuSolver : public LinearSolver {

public:

    /**
     * @brief Construct a new `LinearEigenlibPartialPivLuSolver` object.
     * 
     * @param size Size of the linear system.
     */
    LinearEigenlibPartialPivLuSolver(unsigned int size);

    /**
     * @brief Destroy `LinearEigenlibPartialPivLuSolver` object.
     * 
     */
    virtual ~LinearEigenlibPartialPivLuSolver();


    /**
     * @brief Solve linear system @f$ \mathbf{A} \mathbf{x} = \mathbf{b} @f$.
     * 
     * @param [in] A Vectorized coefficient matrix (C-style).
     * @param [in] b Right-hand-side vector.
     * @param [out] x Solution of the linear system.
     * @return int Error flag.
     */
    virtual int solve(double *A, double *b, double *x);


private:

    MatrixMap matrixMap_;                   ///< Eigen::Map mapping coefficient matrix from Matrix type to raw pointer and vice versa
    VectorMap xVectorMap_;                  ///< Eigen::Map mapping solution vector from Vector type to raw pointer and vice versa
    VectorMap bVectorMap_;                  ///< Eigen::Map mapping r.h.s vector from Vector type to raw pointer and vice versa
    std::unique_ptr<Eigen::PartialPivLU<Matrix> > partialPivLuSolver_;   ///< Pointer to Eigen::PartialPivLU object in Eigen library

    static bool inline registered_ = REGISTER_TYPE(LinearSolverFactory, LinearEigenlibPartialPivLuSolver);

};

#endif // LINSOLVE_HPP