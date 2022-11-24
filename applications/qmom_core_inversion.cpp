/**
 * @file qmom_core_inversion.cpp
 * @author M. Puetz
 * @brief This application compares different algorithms for the computation of the Jacobi matrix from moments in terms of performance and accuracy.
 * @date 2022-10-26
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
#include "io_utils.hpp"
#include "moment_utils.hpp"
#include "build_id.hpp"
#include "core_inversion.hpp"

int main(int argc, char *argv[])
{

    [[maybe_unused]] double doubleEpsilon = 
        std::numeric_limits<double>::epsilon();

    const std::string appName = "qmom_core_inversion";

    // Read command line parameters
    int nMoments = parseArgument<int>(argc, argv, "n_moms");
    std::printf("Running `%s` with %d moments. (%s)\n", appName.data(), nMoments, getTimestamp().c_str());
    std::fflush(stdout);

    int nExecutions = parseArgument<int>(argc, argv, "n_exec");
    int nMomentSets = parseArgument<int>(argc, argv, "n_momsets");
    std::string setupFilename = parseArgument<std::string>(argc, argv, "setup_file");
    std::string momentsFilename = parseArgument<std::string>(argc, argv, "moments_file");
    std::string jacobiMatrixDiagonalFilename = parseArgument<std::string>(argc, argv, "jacobi_diagonal_file");
    std::string jacobiMatrixSuperDiagonalFilename = parseArgument<std::string>(argc, argv, "jacobi_superdiagonal_file");
    std::string outputFilePrefix = parseArgument<std::string>(argc, argv, "outfile_prefix");

    bool computeErrors = true;  // compute errors and output errors by default 
    try {
        computeErrors = parseArgument<bool>(argc, argv, "compute_errors");
    }
    catch (std::runtime_error &e) {
        std::string msg = "Warning: Parameter `compute_errors` not found. Using default `compute_errors=true'.";
        std::fprintf(stderr, "%s\n", msg.c_str());
        std::fflush(stderr);
    }

    // Read algorithms to be tested from setup file
    auto const setup = parseSetupFile(setupFilename);
    const std::vector<std::string> &coreInversionTypes = setup.at("coreInversionTypes");

    // Read moments from given file
    double *moments = static_cast<double *>(
        mkl_malloc(nMomentSets * nMoments * sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(momentsFilename, nMomentSets, nMoments, moments);

    // Create vector of all setups
    std::string coreInversionKey = "coreInversion";
    std::vector<std::map<std::string, std::string>> configurations;
    for (const auto &coreInversionType : coreInversionTypes)
    {
        configurations.push_back(
            {
                {coreInversionKey, coreInversionType},
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
    std::vector<std::string> errorQuantities = {"JacobiMatrix"};
    std::vector<std::string> errorNormTypes = {"Frobenius"};
    std::vector<std::string> errorKeys;
    std::map<std::string, double> errors;
    if (computeErrors)
    {
        for (const auto &errorQuantity : errorQuantities)
        {
            for (const auto &errorNormType : errorNormTypes)
            {
                std::string key = errorQuantity + "RelError" + errorNormType + "Norm";
                errorKeys.push_back(key);
                errors[key] = 0.;
                header += " " + key;
            }
        }
    }
    std::fprintf(dataOutputFile, "%s\n", header.data());

    // Allocate space for Jacobi matrix diagonals
    int jacobiMatrixSize = nMoments / 2;
    double *alpha = static_cast<double *>(mkl_malloc(jacobiMatrixSize * sizeof(double), MALLOC_ALIGN));
    double *gamma = static_cast<double *>(mkl_malloc((jacobiMatrixSize - 1) * sizeof(double), MALLOC_ALIGN));

    // Allocate space for error in moments
    double *momentsError = static_cast<double *>(mkl_malloc(nMoments * sizeof(double), MALLOC_ALIGN));

    // Read reference values
    double *alphaRef = static_cast<double *>(mkl_malloc(nMomentSets * jacobiMatrixSize * sizeof(double), MALLOC_ALIGN));
    double *gammaRef = static_cast<double *>(mkl_malloc(nMomentSets * (jacobiMatrixSize - 1) * sizeof(double), MALLOC_ALIGN));
    // The zeroth gamma coefficient is stored in the reference files but not needed here. So the file contents are
    // temporarily stored in alphaRef (with has sufficient size) and then assigned to gammaRef
    readArrayFromFile(jacobiMatrixSuperDiagonalFilename, nMomentSets, jacobiMatrixSize, alphaRef);
    for (int i = 0; i < nMomentSets; i++)
    {
        double *alphaRefRowBegin = alphaRef + i * jacobiMatrixSize;
        double *gammaRefRowBegin = gammaRef + i * (jacobiMatrixSize - 1);
        for (int j = 0; j < jacobiMatrixSize - 1; j++)
        {
            gammaRefRowBegin[j] = alphaRefRowBegin[j + 1];
        }
    }
    readArrayFromFile(jacobiMatrixDiagonalFilename, nMomentSets, jacobiMatrixSize, alphaRef);

    // Loop over configurations
    configNo = 0;
    for (const auto &configuration : configurations)
    {

        configNo++;
        std::printf("Configuration %d of %lu...\n", configNo, configurations.size());

        // Make core inversion type object
        const std::string &coreInversionType = configuration.at(coreInversionKey);
        std::shared_ptr<CoreInversionAlgorithm> coreInversion =
            CoreInversionAlgorithm::makeShared(coreInversionType, nMoments);

        // Loop over all moment sequences
        for (int i = 0; i < nMomentSets; i++)
        {
            double *mom = &(moments[i * nMoments]);
            double *alphaRefCurr = alphaRef + i * jacobiMatrixSize;
            double *gammaRefCurr = gammaRef + i * (jacobiMatrixSize - 1);

            // Execute number of times specified in command line parameters and compute average CPU time
            auto begin = std::chrono::high_resolution_clock::now();
            for (int j=0; j<nExecutions; j++)
            {
                coreInversion->compute(mom, alpha, gamma);
            }
            auto end = std::chrono::high_resolution_clock::now();
            double cpuTime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
            cpuTime *= 1e-9 / nExecutions;

            // Compute errors (add computation here if any errors are added to `errors`)
            if (computeErrors) {
                // Frobenius norm of the error in the tridiagonal Jacobi matrix
                double &jacobiMatrixRelErrorFrobeniusNorm = errors.at("JacobiMatrixRelErrorFrobeniusNorm");
                jacobiMatrixRelErrorFrobeniusNorm = alpha[0] - alphaRefCurr[0];
                jacobiMatrixRelErrorFrobeniusNorm *= jacobiMatrixRelErrorFrobeniusNorm;
                double jacobiMatrixRefFrobeniusNorm = alphaRefCurr[0]*alphaRefCurr[0];
                for (int k=1; k<jacobiMatrixSize; k++) {
                    // diagonal elements
                    double squaredElement = alpha[k] - alphaRefCurr[k];
                    squaredElement *= squaredElement;
                    jacobiMatrixRelErrorFrobeniusNorm += squaredElement;
                    squaredElement = alphaRefCurr[k]*alphaRefCurr[k];
                    jacobiMatrixRefFrobeniusNorm += squaredElement;

                    // sub-superdiagonal elements
                    squaredElement = gamma[k-1] - gammaRefCurr[k-1];
                    squaredElement *= squaredElement;
                    jacobiMatrixRelErrorFrobeniusNorm += 2*squaredElement;
                    squaredElement = gammaRefCurr[k-1]*gammaRefCurr[k-1];
                    jacobiMatrixRefFrobeniusNorm += 2*squaredElement;
                }
                jacobiMatrixRelErrorFrobeniusNorm /= jacobiMatrixRefFrobeniusNorm;
                jacobiMatrixRelErrorFrobeniusNorm = std::sqrt(jacobiMatrixRelErrorFrobeniusNorm);
            }

            // Write data to output file
            std::fprintf(dataOutputFile, "%d.%d ", configNo, i);
            std::fprintf(dataOutputFile, "%10.9e", cpuTime);
            for (const auto &errorKey : errorKeys)
            {
                std::fprintf(dataOutputFile, " %10.9e", errors.at(errorKey));
            }
            std::fprintf(dataOutputFile, "\n");
        }
    }

    closeOutputFile(dataOutputFile);

    mkl_free(moments);
    mkl_free(alpha);
    mkl_free(alphaRef);
    mkl_free(gamma);
    mkl_free(gammaRef);
    mkl_free(momentsError);

    std::printf("Done. (%s)\n", getTimestamp().c_str());

    return 0;
}
