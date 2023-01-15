/**
 * @file test_qmom.cpp
 * @author M. Puetz
 * @brief Test of the funcionality of all possible QMOM setups.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <cmath>
#include <mkl.h>
#include <vector>
#include <random>

#include "global_defs.hpp"
#include "moment_utils.hpp"
#include "qmom.hpp"
#include "linsolve.hpp"
#include "eigensolve.hpp"
#include "core_inversion.hpp"

#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>


// Dummy rate-of-change function that simply reconstructs the moments
int computeMomentsRateOfChange(double * const nodes, double * const weights, 
    int nNodes, int nMoments, double *momentsRateOfChange)
{
    computeMomentsFromQuadrature(nodes, weights, nNodes, nMoments, momentsRateOfChange);
    return 0;
}

int main ()
{

    int minNumberOfMoments = 4;
    int maxNumberOfMoments = 20;
    const int seed = DEFAULT_RANDOM_SEED;
    double relTol = 1e-10;

    std::mt19937 randomGen;
    randomGen.seed(seed);
    std::uniform_real_distribution<double> uniform(0., 1.);
    std::normal_distribution<double> normal(0., 1.);

    std::vector<std::string> coreInversionTypes = CoreInversionAlgorithmFactory::keys();
    std::vector<std::string> eigenSolverTypes = RealEigenSolverFactory::keys();
    //std::vector<std::string> linearSolverTypes; // = LinearSolverFactory::keys();
    std::vector<std::string> linearSolverTypes = LinearSolverFactory::keys();
    linearSolverTypes.push_back(LinearSolverFactory::noneKey());

    // Test QMOM for defined range of moment sequence lengths
    for (int nMom=minNumberOfMoments; nMom<=maxNumberOfMoments; nMom += 2) {

        // Generate moments from random weights and nodes
        int nNodes = nMom/2;
        double *x = static_cast<double*>(mkl_malloc(nMom/2*sizeof(double), MALLOC_ALIGN));
        double *w = static_cast<double*>(mkl_malloc(nMom/2*sizeof(double), MALLOC_ALIGN));
        double sum = 0.;
        for (int j=0; j<nNodes; j++) {
            x[j] = normal(randomGen);
            w[j] = uniform(randomGen);
            sum += w[j];
        }
        for (int j=0; j<nNodes; j++) {
            w[j] /= sum;
        }
        double *moments = static_cast<double*>(mkl_malloc(nMom*sizeof(double), MALLOC_ALIGN));
        computeMomentsFromQuadrature(x, w, nNodes, nMom, moments);
        mkl_free(x);
        mkl_free(w);

        double *momentsReconst = static_cast<double*>(mkl_malloc(nMom*sizeof(double), MALLOC_ALIGN));
        // Test for all types of core inversion algrithms, eigenvalue algorithms and linear solvers
        for (auto& coreInversionType : coreInversionTypes) {
            for (auto& eigenSolverType : eigenSolverTypes) {
                for (auto& linearSolverType : linearSolverTypes) {
                    
                    std::vector<std::string> qmomTypes = {"QmomStd", "QmomGaG"};
                    for (auto& qmomType : qmomTypes) {

                        auto coreInversionPtr = CoreInversionAlgorithm::makeShared(coreInversionType, nMom);
                        EigenProblemType problemType = EigenProblemType::EigenValsOnly;

                        // If linear solver type is none, use eigenvectors to compute weights
                        if (linearSolverType == LinearSolverFactory::noneKey()) {
                            problemType = EigenProblemType::EigenPairs;
                        }
                        auto eigenSolverPtr = RealEigenSolver::makeShared(eigenSolverType, nMom/2, problemType);
                        auto linearSolverPtr = LinearSolver::makeShared(linearSolverType, nMom/2);

                        auto qmomObj = Qmom::makeShared(qmomType, nMom, coreInversionPtr, 
                            eigenSolverPtr, linearSolverPtr, computeMomentsRateOfChange);
                        qmomObj->compute(moments);
                        qmomObj->computeMoments(0, nMom, momentsReconst);
                        qmomObj->computeMomentsRateOfChange(nMom);

                        if (qmomType == "QmomGaG") {
                            assert(qmomObj->numberOfNodes() == nMom - 1);
                        }
                        else if (qmomType == "QmomGaG") {
                            assert(qmomObj->numberOfNodes() == nMom/2);
                        }
                        // Check if the resulting moments reconstructed from the computed quadrature
                        // are approximately equal to the given moments
                        for (int k=0; k<nMom; k++) {

                            double relError = std::fabs((moments[k] - momentsReconst[k])/moments[k]);
                            assert(relError < relTol);

                            // With the rate-of-change function as defined above, the `momentsRateOfChange`
                            // variable should also equal the original moments
                            relError = std::fabs((moments[k] - qmomObj->momentsRateOfChange()[k])/moments[k]);
                            assert(relError < relTol);
                        }
                    }
                }
            }
        }

        mkl_free(moments);
        mkl_free(momentsReconst);
    }

    return 0;
}
