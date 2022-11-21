/**
 * @file eigensolver_benchmark_moments.cpp
 * @author M. Puetz
 * @brief This application compares different algorithms to compute the eigenvalues of the Jacobi matrix resulting from random moments (given in input files) in terms of performance and accuracy.
 * @date 2022-11-08
 *
 * @copyright Copyright (c) 2022
 *
 */

#include <cmath>
#include <cstdio>
#include <mkl.h>
#include <chrono>
#include <eigen3/Eigen/Eigen>
#include <climits>
#include <stdexcept>

#include "global_defs.hpp"
#include "eigen_defs.hpp"
#include "build_id.hpp"
#include "io_utils.hpp"
#include "moment_utils.hpp"
#include "eigensolve.hpp"
#include "linsolve.hpp"
#include "qmom.hpp"


int main(int argc, char *argv[])
{

    const std::string appName = "eigensolver_benchmark_moments";

    [[maybe_unused]] double doubleEpsilon = 
        std::numeric_limits<double>::epsilon();

    // Read command line parameters
    int nMoments = parseArgument<int>(argc, argv, "n_moms");
    std::printf("Running eigensolver_benchmark_moments with %d moments. (%s)\n", nMoments, getTimestamp().c_str());
    std::fflush(stdout);

    int nExecutions = parseArgument<int>(argc, argv, "n_exec");
    int nMomentSets = parseArgument<int>(argc, argv, "n_momsets");
    std::string setupFilename = parseArgument<std::string>(argc, argv, "setup_file");
    std::string momentsFilename = parseArgument<std::string>(argc, argv, "moments_file");
    std::string jacobiMatrixDiagonalFilename = parseArgument<std::string>(argc, argv, "jacobi_diagonal_file");
    std::string jacobiMatrixSuperDiagonalFilename = parseArgument<std::string>(argc, argv, "jacobi_superdiagonal_file");
    std::string quadratureNodesFilename = parseArgument<std::string>(argc, argv, "quadrature_nodes_file");
    std::string quadratureWeightsFilename = parseArgument<std::string>(argc, argv, "quadrature_weights_file");
    std::string outputFilePrefix = parseArgument<std::string>(argc, argv, "outfile_prefix");

    // Read algorithms to be tested from setup file
    auto const setup = parseSetupFile(setupFilename);
    const std::vector<std::string> &eigenSolverTypes = setup.at("eigenSolverTypes");
    const std::vector<std::string> &linearSolverTypes = setup.at("linearSolverTypes");

    // Read moments from given file
    double *moments = static_cast<double *>(
        mkl_malloc(nMomentSets * nMoments * sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(momentsFilename, nMomentSets, nMoments, moments);

    // Create vector of all setups
    std::string eigenSolverKey = "eigenSolver";
    std::string linearSolverKey = "linearSolver";
    std::vector<std::map<std::string, std::string>> configurations;
    for (const auto &eigenSolverType : eigenSolverTypes)
    {
        for (const auto &linearSolverType : linearSolverTypes) {
            configurations.push_back(
                {
                    {eigenSolverKey, eigenSolverType},
                    {linearSolverKey, linearSolverType}
                });
        }
    }

    // Get timestamp for output files
    std::string timestamp = getTimestamp();

    // Write output file with summary of cases
    std::string summaryOutputFilename = outputFilePrefix;
    if (outputFilePrefix.size() > 0)
    {
        summaryOutputFilename += "_";
    }
    summaryOutputFilename += "summary_nmom" + std::to_string(nMoments) + ".out";
    std::FILE *summaryOutputFile = 
        openOutputFile(summaryOutputFilename, appName, CMAKE_GIT_BUILD_ID);
    std::string header = "ConfigNo";
    for (const auto &type : configurations[0])
    {
        header += " " + type.first;
    }
    std::fprintf(summaryOutputFile, "%s\n", header.data());
    int configNo = 0;
    for (const auto &configuration : configurations)
    {
        std::fprintf(summaryOutputFile, "%d", ++configNo);
        for (const auto &type : configuration)
        {
            std::fprintf(summaryOutputFile, " %s", type.second.data());
        }
        std::fprintf(summaryOutputFile, "\n");
    }
    closeOutputFile(summaryOutputFile);

    // Prepare output file for data (CPU times / errors)
    std::string dataOutputFilename = outputFilePrefix;
    if (outputFilePrefix.size() > 0)
    {
        dataOutputFilename += "_";
    }
    dataOutputFilename += "data_nmom" + std::to_string(nMoments) + ".out";
    std::FILE *dataOutputFile = 
        openOutputFile(dataOutputFilename, appName, CMAKE_GIT_BUILD_ID);
    std::fprintf(dataOutputFile, "# Input files : '%s', '%s', '%s', '%s'\n",
                 setupFilename.c_str(), momentsFilename.c_str(),
                 jacobiMatrixDiagonalFilename.c_str(), jacobiMatrixSuperDiagonalFilename.c_str());
    std::fprintf(dataOutputFile, "# Case summary file : './%s'\n", summaryOutputFilename.c_str());
    std::fprintf(dataOutputFile, "# Number of moments : %d\n", nMoments);
    std::fprintf(dataOutputFile, "# Number of moment sequences : %d\n", nMomentSets);
    std::fprintf(dataOutputFile, "# Number of executions per moment sequence : %d\n#\n", nExecutions);
    header = "CaseNo ComputingTime";

    // Error types to be computed
    std::vector<std::string> outputQuantityKeys = {
                                                        "RelativeSeparationCoefficient", 
                                                        "MomentsRelError2Norm", 
                                                        "MomentsRelErrorInfNorm", 
                                                        "QuadratureNodesRelError2Norm", 
                                                        "QuadratureWeightsRelError2Norm"
                                                    };
    std::map<std::string, double> outputQuantities;
    for (const auto &key : outputQuantityKeys) {
        outputQuantities[key] = 0.;
        header += " " + key;
    }
    std::fprintf(dataOutputFile, "%s\n", header.data());

    // Allocate space for error in moments
    double *momentsError = static_cast<double *>(mkl_malloc(nMoments * sizeof(double), MALLOC_ALIGN));
    VectorMap momentsErrorMap(momentsError, nMoments);

    // Read Jacobi diagonal (alpha) and sub-/superdiagonal (gamma)
    int jacobiMatrixSize = nMoments/2;
    double *alpha = static_cast<double *>(mkl_malloc(nMomentSets * jacobiMatrixSize * sizeof(double), MALLOC_ALIGN));
    double *gamma = static_cast<double *>(mkl_malloc(nMomentSets * jacobiMatrixSize * sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(jacobiMatrixSuperDiagonalFilename, nMomentSets, jacobiMatrixSize, gamma);
    readArrayFromFile(jacobiMatrixDiagonalFilename, nMomentSets, jacobiMatrixSize, alpha);

    // Allocate space for quadrature nodes and weights
    double* quadratureNodes = static_cast<double*>(mkl_malloc(jacobiMatrixSize*sizeof(double), MALLOC_ALIGN));
    double* quadratureWeights = static_cast<double*>(mkl_malloc(jacobiMatrixSize*sizeof(double), MALLOC_ALIGN));
    // Read reference data
    double* quadratureNodesRef = static_cast<double*>(mkl_malloc(nMomentSets*jacobiMatrixSize*sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(quadratureNodesFilename, nMomentSets, jacobiMatrixSize, quadratureNodesRef);
    double* quadratureWeightsRef = static_cast<double*>(mkl_malloc(nMomentSets*jacobiMatrixSize*sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(quadratureWeightsFilename, nMomentSets, jacobiMatrixSize, quadratureWeightsRef);


    // Loop over configurations
    configNo = 0;
    for (const auto &configuration : configurations)
    {

        configNo++;
        std::printf("Configuration %d of %lu...\n", configNo, configurations.size());

        // Make QMOM object
        const std::string &eigenSolverType = configuration.at(eigenSolverKey);
        std::shared_ptr<RealEigenSolver> eigenSolver = 
            RealEigenSolver::makeShared(eigenSolverType, jacobiMatrixSize, EigenProblemType::EigenPairs, EigenValSetType::All);
        const std::string &linearSolverType = configuration.at(linearSolverKey);
        std::shared_ptr<LinearSolver> linearSolver =
            LinearSolver::makeShared(linearSolverType, jacobiMatrixSize);
        QmomStd qmom(nMoments, nullptr, eigenSolver, linearSolver);

        // Loop over all moment sequences
        for (int i = 0; i < nMomentSets; i++)
        {
            double *mom = &(moments[i * nMoments]);
            double *alphaCurr = alpha + i * jacobiMatrixSize;
            double *gammaCurr = gamma + i * jacobiMatrixSize + 1;
            double *quadratureNodesRefCurr = quadratureNodesRef + i*jacobiMatrixSize;
            double *quadratureWeightsRefCurr = quadratureWeightsRef + i*jacobiMatrixSize;

            qmom.setJacobiMatrix(alphaCurr, gammaCurr);
            // Execute number of times specified in command line parameters and compute average CPU time
            qmom.computeQuadrature(moments);    // some involved classes initialize data during the very first call
            auto begin = std::chrono::high_resolution_clock::now();
            for (int j=0; j<nExecutions; j++)
            {
                qmom.computeQuadrature(mom);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double cpuTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            cpuTime *= 1e-9 / nExecutions;

            auto quadratureNodes = qmom.quadratureNodes();
            auto quadratureWeights = qmom.quadratureWeights();

            // Compute separation of eigenvalues
            double minDistance = quadratureNodesRefCurr[1] - quadratureNodesRefCurr[0];
            double maxDistance = minDistance;
            // This works for sorted eigenvalues
            for (int j=2; j<jacobiMatrixSize; j++) {
                double distance = quadratureNodesRefCurr[j] - quadratureNodesRefCurr[j-1];
                minDistance = std::min(minDistance, distance);
                maxDistance = std::max(maxDistance, distance);
            }
            outputQuantities.at("RelativeSeparationCoefficient") = minDistance/maxDistance;

            // Compute moment errors
            qmom.computeMomentsRateOfChange(nMoments, [](double){return 1;});
            auto reconstructedMoments = qmom.momentsRateOfChange();
            for (int j=0; j<nMoments; j++) {
                momentsError[j] = std::abs((reconstructedMoments[j] - mom[j])/mom[j]);
            }
            outputQuantities.at("MomentsRelError2Norm") = momentsErrorMap.norm();
            outputQuantities.at("MomentsRelErrorInfNorm") = 
                momentsErrorMap.lpNorm<Eigen::Infinity>();

            // Compute quadrature errors 
            double &quadratureNodesRelError2Norm = outputQuantities.at("QuadratureNodesRelError2Norm");
            quadratureNodesRelError2Norm = 0;
            double &quadratureWeightsRelError2Norm = outputQuantities.at("QuadratureWeightsRelError2Norm");
            quadratureWeightsRelError2Norm = 0;
            for (int j=0; j<jacobiMatrixSize; j++) {
                double error = quadratureNodes[j] - quadratureNodesRefCurr[j];
                quadratureNodesRelError2Norm += error*error;
                error = quadratureWeights[j] - quadratureWeightsRefCurr[j];
                quadratureWeightsRelError2Norm += error*error;

            }
            outputQuantities.at("QuadratureNodesRelError2Norm") = 
                std::sqrt(quadratureNodesRelError2Norm) / cblas_dnrm2(jacobiMatrixSize, quadratureNodesRefCurr, 1);
            outputQuantities.at("QuadratureWeightsRelError2Norm") = 
                std::sqrt(quadratureWeightsRelError2Norm) / cblas_dnrm2(jacobiMatrixSize, quadratureWeightsRefCurr, 1);

            // Write data to output file
            std::fprintf(dataOutputFile, "%d.%d ", configNo, i);
            std::fprintf(dataOutputFile, "%10.9e", cpuTime);
            for (const auto &key : outputQuantityKeys)
            {
                std::fprintf(dataOutputFile, " %10.9e", outputQuantities[key]);
            }
            std::fprintf(dataOutputFile, "\n");
        }
    }

    closeOutputFile(dataOutputFile);

    mkl_free(moments);
    mkl_free(alpha);
    mkl_free(gamma);
    mkl_free(momentsError);
    mkl_free(quadratureNodes);
    mkl_free(quadratureWeights);
    mkl_free(quadratureNodesRef);
    mkl_free(quadratureWeightsRef);

    std::printf("Done. (%s)\n", getTimestamp().c_str());

    return 0;
}