/**
 * @file qmom_linsolver.cpp
 * @author M. Puetz
 * @brief This application compares different linear solvers to compute the weights corresponding to a set of moments from the Vandermonde matrix of quadrature nodes.
 * @date 2022-11-05
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

#include "build_id.hpp"
#include "io_utils.hpp"
#include "eigen_defs.hpp"
#include "moment_utils.hpp"
#include "linsolve.hpp"


int main(int argc, char *argv[])
{

    const std::string appName = "qmom_linsolver";

    [[maybe_unused]] double doubleEpsilon = 
        std::numeric_limits<double>::epsilon();

    // Read command line parameters
    int nMoments = parseArgument<int>(argc, argv, "n_moms");
    std::printf("Running `%s` with %d moments. (%s)\n", appName.data(), nMoments, getTimestamp().c_str());
    std::fflush(stdout);

    int nExecutions = parseArgument<int>(argc, argv, "n_exec");
    int nMomentSets = parseArgument<int>(argc, argv, "n_momsets");
    std::string setupFilename = parseArgument<std::string>(argc, argv, "setup_file");
    std::string momentsFilename = parseArgument<std::string>(argc, argv, "moments_file");
    std::string quadratureNodesFilename = parseArgument<std::string>(argc, argv, "quadrature_nodes_file");
    std::string quadratureWeightsFilename = parseArgument<std::string>(argc, argv, "quadrature_weights_file");
    std::string outputFilePrefix = parseArgument<std::string>(argc, argv, "outfile_prefix");


    // Read algorithms to be tested from setup file
    auto const setup = parseSetupFile(setupFilename);
    const std::vector<std::string> &linearSolverTypes = setup.at("linearSolverTypes");

    // Read moments from given file
    double *moments = static_cast<double *>(
        mkl_malloc(nMomentSets * nMoments * sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(momentsFilename, nMomentSets, nMoments, moments);

    // Create vector of all setups
    std::string linearSolverKey = "linearSolver";
    std::vector<std::map<std::string, std::string>> configurations;
    for (const auto &linearSolverType : linearSolverTypes)
    {
        configurations.push_back(
            {
                {linearSolverKey, linearSolverType},
            });
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
    std::FILE *dataOutputFile = openOutputFile(dataOutputFilename, appName, CMAKE_GIT_BUILD_ID);
    std::fprintf(dataOutputFile, "# Input files : '%s', '%s', '%s', '%s'\n",
                 setupFilename.c_str(), momentsFilename.c_str(),
                 quadratureNodesFilename.c_str(), quadratureWeightsFilename.c_str());
    std::fprintf(dataOutputFile, "# Case summary file : './%s'\n", summaryOutputFilename.c_str());
    std::fprintf(dataOutputFile, "# Number of moments : %d\n", nMoments);
    std::fprintf(dataOutputFile, "# Number of moment sequences : %d\n", nMomentSets);
    std::fprintf(dataOutputFile, "# Number of executions per moment sequence : %d\n#\n", nExecutions);
    header = "CaseNo ComputingTime";

    // Output quantities to be computed
    std::vector<std::string> outputQuantityKeys = {"VandermondeMatrixConditionNumber1Norm",
                                                   "VandermondeMatrixConditionNumber2Norm",
                                                   "VandermondeMatrixConditionNumberInfNorm",
                                                   "MomentsRelError2Norm",
                                                   "MomentsRelErrorInfNorm",
                                                   "WeightsRelError2Norm"
    };
    std::map<std::string, double> outputQuantities;
    for (const auto &outputQuantityKey : outputQuantityKeys)
    {
        outputQuantities[outputQuantityKey] = 0.;
        header += " " + outputQuantityKey;
    }
    std::fprintf(dataOutputFile, "%s\n", header.data());

    // Allocate space for quadrature weights (nodes are not computed but provided as input)
    int nQuadratureNodes = nMoments / 2;
    double *quadratureWeights = static_cast<double*>(
        mkl_malloc(nQuadratureNodes*sizeof(double), MALLOC_ALIGN));
    VectorMap quadratureWeightsMap(quadratureWeights, nQuadratureNodes);

    // Allocate space for moments and errors
    double *momentsError = static_cast<double *>(mkl_malloc(nMoments*sizeof(double), MALLOC_ALIGN));
    VectorMap momentsErrorMap(momentsError, nMoments);

    // Read reference values
    double *quadratureNodesRef = static_cast<double*>(
        mkl_malloc(nMomentSets*nQuadratureNodes*sizeof(double), MALLOC_ALIGN));
    double *quadratureWeightsRef = static_cast<double*>(
        mkl_malloc(nMomentSets*nQuadratureNodes*sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(quadratureNodesFilename, nMomentSets, nQuadratureNodes, quadratureNodesRef);
    readArrayFromFile(quadratureWeightsFilename, nMomentSets, nQuadratureNodes, quadratureWeightsRef);

    // Allocate space for Vandermonde matrix and create corresponding Eigen::Map
    int VandermondeMatrixSize = nQuadratureNodes*nQuadratureNodes;
    double *VandermondeMatrix = static_cast<double*>(
        mkl_malloc(VandermondeMatrixSize*sizeof(double), MALLOC_ALIGN));
    MatrixMap VandermondeMatrixMap(VandermondeMatrix, nQuadratureNodes, nQuadratureNodes);
    VandermondeMatrixMap.setOnes();     // First row constantly contains only ones
    
    // Loop over configurations
    configNo = 0;
    for (const auto &configuration : configurations)
    {

        configNo++;
        std::printf("Configuration %d of %lu...\n", configNo, configurations.size());

        // Make core inversion type object
        const std::string &linearSolverType = configuration.at(linearSolverKey);
        unsigned size = nQuadratureNodes;
        std::shared_ptr<LinearSolver> linearSolver =
            LinearSolver::makeShared(linearSolverType, size);

        // Loop over all moment sequences
        for (int i = 0; i < nMomentSets; i++)
        {
            double *mom = moments + i*nMoments;
            double *quadratureNodesRefCurr = quadratureNodesRef + i*nQuadratureNodes;
            VectorMap quadratureWeightsRefCurrMap(quadratureWeightsRef + i*nQuadratureNodes, 
                nQuadratureNodes);

            // Assemble Vandermonde matrix from reference quadrature nodes
            for (int j=1; j<nQuadratureNodes; j++) {
                for (int k=0; k<nQuadratureNodes; k++) {
                    VandermondeMatrixMap(j,k) = 
                        VandermondeMatrixMap(j-1,k)*quadratureNodesRefCurr[k];
                }
            }

            // Compute condition numbers
            // Spectral condition number from SVD
            Eigen::JacobiSVD<Matrix> svd(VandermondeMatrixMap);
            double conditionNumber = svd.singularValues()(0)
                / svd.singularValues()(nQuadratureNodes - 1);
            outputQuantities.at("VandermondeMatrixConditionNumber2Norm") = conditionNumber;

            // Estimate for 1-norm condition number
            int info = LAPACKE_dgecon(LAPACK_ROW_MAJOR, '1', 
                nQuadratureNodes, VandermondeMatrix, nQuadratureNodes, 
                VandermondeMatrixMap.cwiseAbs().colwise().sum().maxCoeff(), &conditionNumber);
            outputQuantities.at("VandermondeMatrixConditionNumber1Norm") = 1/conditionNumber;
            if (info != 0) {
                std::fprintf(stderr, "Warning: A problem occurred while computing the condition number (1-norm),"
                    " `LAPACKE_dgecon` returned `info=%d`.\n", info);
                std::fflush(stdout);
                std::fflush(stderr);
            }

            // Estimate for infinity-norm condition number
            info = LAPACKE_dgecon(LAPACK_ROW_MAJOR, 'I', nQuadratureNodes, VandermondeMatrix, 
                nQuadratureNodes, VandermondeMatrixMap.cwiseAbs().rowwise().sum().maxCoeff(), 
                &conditionNumber);
            outputQuantities.at("VandermondeMatrixConditionNumberInfNorm") = 1/conditionNumber;
            if (info != 0) {
                std::fprintf(stderr, "Warning: A problem occurred while computing the condition number (inf-norm),"
                    " `LAPACKE_dgecon` returned `info=%d`.\n", info);
                std::fflush(stdout);
                std::fflush(stderr);
            }

            // Execute number of times specified in command line parameters and compute average CPU time
            linearSolver->solve(VandermondeMatrix, mom, quadratureWeights); // do not measure time for initializations during first call
            auto begin = std::chrono::high_resolution_clock::now();
            for (int j=0; j<nExecutions; j++)
            {
                linearSolver->solve(VandermondeMatrix, mom, quadratureWeights);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double cpuTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            cpuTime *= 1e-9 / nExecutions;

            // Compute moment errors
            computeMomentsFromQuadrature(quadratureNodesRefCurr, quadratureWeights,
                nQuadratureNodes, nMoments, momentsError);
            for (int j=0; j<nMoments; j++) {
                momentsError[j] = std::abs((momentsError[j] - mom[j])/mom[j]);
            }
            outputQuantities.at("MomentsRelError2Norm") = momentsErrorMap.norm();
            outputQuantities.at("MomentsRelErrorInfNorm") = 
                momentsErrorMap.lpNorm<Eigen::Infinity>();

            // Compute weight error norm
            outputQuantities.at("WeightsRelError2Norm") =
                (quadratureWeightsMap - quadratureWeightsRefCurrMap).norm();

            // Write data to output file
            std::fprintf(dataOutputFile, "%d.%d ", configNo, i);
            std::fprintf(dataOutputFile, "%10.9e", cpuTime);
            for (const auto &key : outputQuantityKeys)
            {
                std::fprintf(dataOutputFile, " %10.9e", outputQuantities.at(key));
            }
            std::fprintf(dataOutputFile, "\n");
        }
    }

    closeOutputFile(dataOutputFile);

    mkl_free(moments);
    mkl_free(quadratureNodesRef);
    mkl_free(quadratureWeights);
    mkl_free(quadratureWeightsRef);
    mkl_free(momentsError);
    mkl_free(VandermondeMatrix);

    std::printf("Done. (%s)\n", getTimestamp().c_str());

    return 0;
}
