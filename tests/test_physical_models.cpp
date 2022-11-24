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


/**
 * @brief Get specified number (up to 20) of reference moments used for case 3 in @cite Puetz2022 with
 *  the parameters used by default, see `defaultModelConstants` in 'physical_models.hpp'.
 * 
 * @param nMoments Number of moments.
 * @return std::shared_ptr<double[]> Shared pointer to moments.
 */
std::shared_ptr<double[]> getReferenceMoments(int nMoments)
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

void test_fokker_planck_equation()
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
    auto fokkerPlanck = PhysicalModelFactory::makeShared("FokkerPlanckEquation");

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

    auto moments = getReferenceMoments(nMoments);

    QmomPtr->compute(moments.get());

    auto momentsRateOfChange = QmomPtr->momentsRateOfChange();

    for (int k=0; k<nMoments; k++) {
        // For this specific setup the source terms are zero (see Puetz2022)
        assert(std::abs(momentsRateOfChange[k]) < absTol);
    }
}

int main()
{

    test_fokker_planck_equation();

    return 0;
}