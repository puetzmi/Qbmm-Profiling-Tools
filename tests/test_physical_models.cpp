/**
 * @file test_physical_models.cpp
 * @author M. Puetz
 * @brief Test functionality of classes derived from `PhysicalModel` class.
 * @date 2022-11-23
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <cstdio>
#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>
#include <memory>
#include <functional>
#include <mkl.h>

#include "global_defs.hpp"
#include "core_inversion.hpp"
#include "eigensolve.hpp"
#include "linsolve.hpp"
#include "qmom.hpp"
#include "physical_models.hpp"
#include "moment_utils.hpp"


/**
 * @brief Get specified number (up to 20) of reference moments used for case 3 in @cite Puetz2022 with
 *  the parameters used by default, see `defaultModelConstants` in 'physical_models.hpp'.
 * 
 * @param nMoments Number of moments.
 * @return std::shared_ptr<double[]> Shared pointer to moments.
 */
std::shared_ptr<double[]> getFpeReferenceMoments(int nMoments)
{

    double allMoments[20] =
    {
        1.000000000000000e+00,
        0.000000000000000e+00,
        1.939393939393939e-01,
        0.000000000000000e+00,
        1.880624426078972e-01,
        0.000000000000000e+00,
        3.282544452792387e-01,
        0.000000000000000e+00,
        8.275990862797774e-01,
        0.000000000000000e+00,
        2.728569108704237e+00,
        0.000000000000000e+00,
        1.111271782454089e+01,
        0.000000000000000e+00,
        5.387984399777400e+01,
        0.000000000000000e+00,
        3.030333044238441e+02,
        0.000000000000000e+00
    };

    std::shared_ptr<double[]> moments
    (
        static_cast<double*>(mkl_malloc(nMoments*sizeof(double), MALLOC_ALIGN)),
        &mkl_free
    );

    for (int k=0; k<nMoments; k++) {
        moments[k] = allMoments[k];
    }

    return moments;
}

void testFokkerPlanckEquationConstantDrag()
{

    int nMoments = 6;
    std::string qmomType = "QmomGaG";
    double absTol = 1e-10;


    auto coreInversion = CoreInversionAlgorithm::makeShared(
        CoreInversionAlgorithmFactory::keys()[0],
        nMoments
    );

    auto eigenSolver = RealEigenSolver::makeShared(
        RealEigenSolverFactory::keys()[0],
        nMoments/2,
        EigenProblemType::EigenPairs
    );
    auto fokkerPlanck = PhysicalModelFactory::makeShared("FokkerPlanckConstantCd");

    std::function<int(double * const, double * const, int, int, double*)> 
        momentsRateOfChangeFunction
        =
        [&](double * const nodes, double * const weights, int nNodes, 
            int nMoments, double* momentsRateOfChange) {
                return fokkerPlanck->computeMomentsRateOfChange(
                    nodes, weights, nNodes, nMoments, momentsRateOfChange
                );
        };

    auto QmomPtr = Qmom::makeShared("QmomGaG", nMoments, coreInversion,
        eigenSolver, nullptr, momentsRateOfChangeFunction);

    auto moments = getFpeReferenceMoments(nMoments);

    QmomPtr->compute(moments.get());

    auto momentsRateOfChange = QmomPtr->momentsRateOfChange();

    for (int k=0; k<nMoments; k++) {
        // For this specific setup the source terms are zero (see Puetz2022)
        assert(std::abs(momentsRateOfChange[k]) < absTol);
    }
}


void testFokkerPlanckEquationVelocityDependentDrag()
{

    int nMomentsMin = 4;
    int nMomentsMax = 12;
    std::string qmomType = "QmomGaG";
    double absTol = 1e-7;

    // Use normalized Chebyshev moments as input for simplicity
    std::shared_ptr<double[]> moments
    (
        static_cast<double*>(mkl_malloc(nMomentsMax*sizeof(double), MALLOC_ALIGN)),
        &mkl_free
    );
    for (int k=0; k<nMomentsMax; k++) {
        moments[k] = getChebyshevMoment(k);
    }
    for (int k=0; k<nMomentsMax; k++) {
        moments[k] /= moments[0];
    }

    for (int nMoments=nMomentsMin; nMoments<=nMomentsMax; nMoments += 2)
    {
        auto coreInversion = CoreInversionAlgorithm::makeShared(
            CoreInversionAlgorithmFactory::keys()[0],
            nMoments
        );

        auto eigenSolver = RealEigenSolver::makeShared(
            RealEigenSolverFactory::keys()[0],
            nMoments/2,
            EigenProblemType::EigenPairs
        );
        auto fokkerPlanck = PhysicalModelFactory::makeShared("FokkerPlanckVelocityDependentCd");

        std::function<int(double * const, double * const, int, int, double*)> 
            momentsRateOfChangeFunction
            =
            [&](double * const nodes, double * const weights, int nNodes, 
                int nMoments, double* momentsRateOfChange) {
                    return fokkerPlanck->computeMomentsRateOfChange(
                        nodes, weights, nNodes, nMoments, momentsRateOfChange
                    );
            };

        auto QmomPtr = Qmom::makeShared("QmomGaG", nMoments, coreInversion,
            eigenSolver, nullptr, momentsRateOfChangeFunction);

        QmomPtr->compute(moments.get());

        auto momentsRateOfChange = QmomPtr->momentsRateOfChange();

        // Zeroth moment must remain constant
        assert(std::abs(momentsRateOfChange[0]) < absTol);

        // All odd moments remain constant (zero)
        for (int k=1; k<nMoments; k+=2) {
            assert(std::abs(momentsRateOfChange[k]) < absTol);
        }

        // With the given setup, the second moment should decrease
        assert(momentsRateOfChange[2] < 0);

    }
}


int main()
{

    testFokkerPlanckEquationConstantDrag();

    testFokkerPlanckEquationVelocityDependentDrag();

    return 0;
}