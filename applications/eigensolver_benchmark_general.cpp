/**
 * @file eigensolver_benchmark_general.cpp
 * @author M. Puetz
 * @brief This application compares the performance and accuracy of different eigen solvers for different tridiagonal Jacobi matrices.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <cstdio>
#include <vector>
#include <map>
#include <random>
#include <chrono>
#include "build_id.hpp"
#include "io_utils.hpp"
#include "eigensolve.hpp"


int main (int argc, char *argv[]) {

    const std::string appName = "eigensolver_benchmark_general";

    // Default parameters
    int defaultRandomSeed = DEFAULT_RANDOM_SEED;
    std::string defaultOutputFilename = "eigensolver_benchmark.out";


    // Read parameters from command line arguments
    int randomSeed;
    std::string outputFilename;

    try {
        randomSeed = parseArgument<int>(argc, argv, "random_seed");
    }
    catch (std::runtime_error const &e) {
        randomSeed = defaultRandomSeed;
    } 

    try {
        outputFilename = parseArgument<std::string>(argc, argv, "outfile");
    }
    catch (std::runtime_error const &e) {
        outputFilename = defaultOutputFilename;
    } 

    int nMatrices = parseArgument<int>(argc, argv, "n_matrices");
    int nMin = parseArgument<int>(argc, argv, "n_min");
    int nMax = parseArgument<int>(argc, argv, "n_max");
    int nSizes = parseArgument<int>(argc, argv, "n_sizes");
    int nExecutions = parseArgument<int>(argc, argv, "n_exec");

    bool computeEigenvectors = parseArgument<bool>(argc, argv, "compute_eigenvectors");
    EigenProblemType problemType = EigenProblemType::EigenValsOnly;
    if (computeEigenvectors) {
        problemType = EigenProblemType::EigenPairs;
    }

    // Get all available EigenSolver types
    const std::vector<std::string> &eigenSolverTypes = RealEigenSolverFactory::keys();

    // Prepare file for output
    std::FILE *outputFile = openOutputFile(outputFilename, appName, CMAKE_GIT_BUILD_ID);
    std::fprintf(outputFile, "n");
    for (const auto &eigenSolverType : eigenSolverTypes) {
        std::fprintf(outputFile, " %s", eigenSolverType.data());
    }
    std::fprintf(outputFile, "\n");


    // Create random symmetric tridiagonal matrix
    double mu = 0.;     // mean
    double sigma = 1.;  // standard deviation
    std::mt19937 randomNumberGenerator(randomSeed);
    std::normal_distribution<double> normalDistribution(mu, sigma);
    std::gamma_distribution<double> gammaDistribution(0.5);

    double *mainDiagonals = static_cast<double*>(mkl_malloc(nMatrices*nMax*sizeof(double), MALLOC_ALIGN));
    double *superDiagonals = static_cast<double*>(mkl_malloc(nMatrices*(nMax-1)*sizeof(double), MALLOC_ALIGN));

    // Generate all values of main and sub-/superdiagonal elements and write to output files
    std::printf("Generating random symmetric tridiagonal matrices...\n");
    std::FILE *mainDiagonalOutputFile = 
        openOutputFile("diag0.out", appName, CMAKE_GIT_BUILD_ID);
    std::FILE *superDiagonalOutputFile =
        openOutputFile("diag1.out", appName, CMAKE_GIT_BUILD_ID);
    for (int i=0; i<nMatrices; i++) {
        double *mainDiagonal = mainDiagonals + i*nMax;
        double *superDiagonal = superDiagonals + i*(nMax - 1);
        for (int j=0; j<nMax-1; j++) {
            mainDiagonal[j] = normalDistribution(randomNumberGenerator);
            superDiagonal[j] = gammaDistribution(randomNumberGenerator);
            std::fprintf(mainDiagonalOutputFile, "%16.15e ", mainDiagonal[j]);
            std::fprintf(superDiagonalOutputFile, "%16.15e ", superDiagonal[j]);
        }
        mainDiagonal[nMax-1] = normalDistribution(randomNumberGenerator);
        std::fprintf(mainDiagonalOutputFile, "%16.15e\n", mainDiagonal[nMax-1]);
        std::fprintf(superDiagonalOutputFile, "\n");
    }

    closeOutputFile(mainDiagonalOutputFile);
    closeOutputFile(superDiagonalOutputFile);

    // Linearly uniform spacing
    std::vector<int> nAll(nSizes);
    double delta = (nMax - nMin)/(nSizes - 1);
    for (int i=0; i<nSizes; i++) {
        nAll[i] = static_cast<int>(nMin + i*delta + 0.5);   // round to int
    }

    // Run eigenvalue algorithms and measure CPU times
    std::printf("Measuring computing times of eigenvalue algorithms...\n");
    for (const int n : nAll) {
        std::printf("n = %d\n", n);

        for (int i=0; i<nMatrices; i++) {

            std::fprintf(outputFile, "%d", n);
            
            double *mainDiagonal = mainDiagonals + i*nMax;
            double *superDiagonal = superDiagonals + i*(nMax - 1);

            for (auto const &eigenSolverType : eigenSolverTypes) {
                auto eigenSolver = RealEigenSolver::makeUnique(eigenSolverType, n, problemType);

                // This is done because some classes initialize variables or allocate memory during the first call
                eigenSolver->compute(mainDiagonal, superDiagonal);

                auto begin = std::chrono::high_resolution_clock::now();
                for (int i=0; i<nExecutions; i++) {
                    eigenSolver->compute(mainDiagonal, superDiagonal);
                }
                auto end = std::chrono::high_resolution_clock::now();
                double runtime = 1e-9*std::chrono::duration_cast<std::chrono::nanoseconds>(end - begin).count();
                std::fprintf(outputFile, " %10.9e", runtime/nExecutions);
            }
            std::fprintf(outputFile, "\n");
        }
    }
    std::printf("... done!");

    mkl_free(mainDiagonals);
    mkl_free(superDiagonals);

    closeOutputFile(outputFile);

    return 0;
}
