/**
 * @file qmom_subroutines_profiler.cpp
 * @author M. Puetz
 * @brief This application measures the computation times for each subroutine in the QMOM/GaG-QMOM algorithms.
 * @date 2022-12-23
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
#include "physical_models.hpp"
#include "qmom.hpp"
#include "qmom_profiler.hpp"


int main(int argc, char *argv[])
{

    const std::string appName = "qmom_subroutines_profiler";

    [[maybe_unused]] double doubleEpsilon = 
        std::numeric_limits<double>::epsilon();

    // Tolerances for verification of results
    double absTol = 1e-10;
    double relTol = 1e-6;

    // Read command line parameters
    int nMoments = parseArgument<int>(argc, argv, "n_moms");
    std::printf("Running `%s` with %d moments. (%s)\n", appName.data(), nMoments, getTimestamp().data());
    std::fflush(stdout);

    int nExecutions = parseArgument<int>(argc, argv, "n_exec");
    int nMomentSets = parseArgument<int>(argc, argv, "n_momsets");
    std::string setupFilename = parseArgument<std::string>(argc, argv, "setup_file");
    std::string momentsFilename = parseArgument<std::string>(argc, argv, "moments_file");
    std::string outputFilePrefix = parseArgument<std::string>(argc, argv, "outfile_prefix");

    // When solving linear system the error can get relatively large with a lot of moments...
    relTol *= nMoments*nMoments;

    // Read algorithms / physical models to be tested from setup file
    auto const setup = parseSetupFile(setupFilename);
    const std::vector<std::string> &qmomTypes = setup.at("qmomTypes");
    const std::vector<std::string> &coreInversionTypes = setup.at("coreInversionTypes");
    const std::vector<std::string> &eigenSolverTypes = setup.at("eigenSolverTypes");
    const std::vector<std::string> &linearSolverTypes = setup.at("linearSolverTypes");
    const std::vector<std::string> &physicalModelTypes = setup.at("physicalModelTypes");

    // Read moments from given file
    double *moments = static_cast<double *>(
        mkl_malloc(nMomentSets * nMoments * sizeof(double), MALLOC_ALIGN));
    readArrayFromFile(momentsFilename, nMomentSets, nMoments, moments);

    // Create vector of all configurations
    std::string qmomTypeKey = "qmomType";
    std::string coreInversionKey = "coreInversion";
    std::string eigenSolverKey = "eigenSolver";
    std::string linearSolverKey = "linearSolver";
    std::string physicalModelKey = "physicalModel";
    std::vector<std::map<std::string, std::string>> configurations;
    for (const auto &qmomType : qmomTypes) {
        for (const auto &coreInversionType : coreInversionTypes) {
            for (const auto &eigenSolverType : eigenSolverTypes) {
                for (const auto &linearSolverType : linearSolverTypes) {
                    for (const auto &physicalModelType : physicalModelTypes) {

                        configurations.push_back(
                            {
                                {qmomTypeKey, qmomType},
                                {coreInversionKey, coreInversionType},
                                {eigenSolverKey, eigenSolverType},
                                {linearSolverKey, linearSolverType},
                                {physicalModelKey, physicalModelType}
                            }
                        );
                    }
                }
            }
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
    std::fprintf(dataOutputFile, "# Input files : '%s', '%s'\n",
                 setupFilename.data(), momentsFilename.data());
    std::fprintf(dataOutputFile, "# Case summary file : './%s'\n", summaryOutputFilename.data());
    std::fprintf(dataOutputFile, "# Number of moments : %d\n", nMoments);
    std::fprintf(dataOutputFile, "# Number of moment sequences : %d\n", nMomentSets);
    std::fprintf(dataOutputFile, "# Number of executions per moment sequence : %d\n#\n", nExecutions);
    header = "CaseNo";
    std::map<std::string, double> outputQuantities;
    for (const auto &subroutineName : QmomProfiler::subroutineNames()) {
        header += " " + subroutineName;
    }
    std::fprintf(dataOutputFile, "%s\n", header.data());

    int jacobiMatrixSize = nMoments/2;

    double *reconstructedMoments = static_cast<double*>(
            mkl_malloc(nMoments*sizeof(double), MALLOC_ALIGN)
    );

    // Loop over configurations
    configNo = 0;
    for (const auto &configuration : configurations)
    {

        configNo++;
        std::printf("Configuration %d of %lu...\n", configNo, configurations.size());

        // Make QMOM profiler object
        const std::string &coreInversionType = configuration.at(coreInversionKey);
        std::shared_ptr<CoreInversionAlgorithm> coreInversion =
        CoreInversionAlgorithm::makeShared(coreInversionType, nMoments);

        const std::string &eigenSolverType = configuration.at(eigenSolverKey);
        std::shared_ptr<RealEigenSolver> eigenSolver = 
            RealEigenSolver::makeShared(eigenSolverType, jacobiMatrixSize, EigenProblemType::EigenPairs, EigenValSetType::All);

        const std::string &linearSolverType = configuration.at(linearSolverKey);
        std::shared_ptr<LinearSolver> linearSolver =
            LinearSolver::makeShared(linearSolverType, jacobiMatrixSize);

        const std::string &physicalModelType = configuration.at(physicalModelKey);
        std::shared_ptr<PhysicalModel> physicalModel =
            PhysicalModel::makeShared(physicalModelType);

        std::function<int(double * const, double * const, int, int, double*)> 
            computeMomentsRateOfChange
            =
            [&](double * const nodes, double * const weights, int nNodes, 
                int nMoments, double* momentsRateOfChange) {
                    return physicalModel->computeMomentsRateOfChange(
                        nodes, weights, nNodes, nMoments, momentsRateOfChange
                    );
            };
        
        const std::string &qmomType = configuration.at(qmomTypeKey);
        std::shared_ptr<Qmom> qmomPtr = Qmom::makeShared(qmomType, nMoments, 
            coreInversion, eigenSolver, linearSolver, computeMomentsRateOfChange);
        
        QmomProfiler qmomProfiler(nExecutions, qmomPtr);

        // Loop over all moment sequences
        for (int i = 0; i < nMomentSets; i++)
        {
            double *mom = &(moments[i * nMoments]);

            qmomProfiler.compute(mom);

            int nNodes = qmomPtr->numberOfNodes();
            computeMomentsFromQuadrature(
            qmomPtr->quadratureNodes().get(),
                    qmomPtr->quadratureWeights().get(), nNodes, nMoments, reconstructedMoments);
                    
            // Ensure correctness of the results
            for (int k=0; k<nMoments; k++) {
                bool ok = std::abs(mom[k] - reconstructedMoments[k]) < absTol + relTol*std::abs(mom[k]);
                // this should never happen
                if (!ok) {
                    throw std::runtime_error("Moment comparison failed.");
                }
            }

            // Write data to output file
            std::fprintf(dataOutputFile, "%d.%d ", configNo, i);
            for (const auto &subroutineName : qmomProfiler.subroutineNames())
            {
                std::fprintf(dataOutputFile, " %10.9e", qmomProfiler.cpuTimes().at(subroutineName));
            }
            std::fprintf(dataOutputFile, "\n");
        }
    }

    closeOutputFile(dataOutputFile);

    mkl_free(moments);
    mkl_free(reconstructedMoments);

    std::printf("Done. (%s)\n", getTimestamp().data());

    return 0;
}
