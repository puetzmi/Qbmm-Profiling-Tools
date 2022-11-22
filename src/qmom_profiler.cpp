/**
 * @file qmom_profiler.cpp
 * @author M. Puetz
 * @brief Implementation of QMOM profiler routines.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include <map>
#include <memory>
#include <chrono>
#include "qmom_profiler.hpp"



QmomProfiler::QmomProfiler(int numberOfExecutions, std::shared_ptr<Qmom> qmomPtr)
    :
    nExecutions_(numberOfExecutions),
    cpuTimes_(),
    clock_(),
    qmomPtr_(qmomPtr)
{
    for (auto& key : subroutineNames_) {
        cpuTimes_[key] = 0.;
    }
}

QmomProfiler::~QmomProfiler()
{
}


template <class... ArgTypes>
double QmomProfiler::qmomSubroutineRuntime(std::function<int(ArgTypes...)> callableSubroutine, 
    ArgTypes... args)
{
    [[maybe_unused]] int info;
    double totalRuntime, averageRuntime;

    auto start = clock_.now();
    for (int i=0; i<nExecutions_; i++) {
        info = callableSubroutine(args...);
    }
    auto end = clock_.now();

    // runtime in seconds
    totalRuntime = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
    averageRuntime = totalRuntime/nExecutions_;

    return averageRuntime;
}


int QmomProfiler::computeJacobiMatrix(double *moments, bool resetCpuTime)
{

    const std::string cpuTimesKey = "CoreInversion";

    std::function<int(double*)> callableSubroutine = 
        std::bind(&Qmom::computeJacobiMatrix, qmomPtr_, std::placeholders::_1);

    if (resetCpuTime)
        cpuTimes_.at(cpuTimesKey) = 0;

    cpuTimes_.at(cpuTimesKey) += qmomSubroutineRuntime<double*>(callableSubroutine, moments);

    return 0;
}


int QmomProfiler::computeQuadrature(double *moments, bool resetCpuTime)
{

    const std::string cpuTimesKey = "Quadrature";

    std::function<int(double*)> callableSubroutine = 
        std::bind(&Qmom::computeQuadrature, qmomPtr_, std::placeholders::_1);

    if (resetCpuTime)
        cpuTimes_.at(cpuTimesKey) = 0;

    cpuTimes_.at(cpuTimesKey) += qmomSubroutineRuntime<double*>(callableSubroutine, moments);

    return 0;
}


int QmomProfiler::computeMomentsRateOfChange(int nMoments, bool resetCpuTime)
{
    const std::string cpuTimesKey = "Integrate";

    std::function<int(int)> callableSubroutine = std::bind(
        &Qmom::computeMomentsRateOfChange, qmomPtr_, std::placeholders::_1);

    if (resetCpuTime)
        cpuTimes_.at(cpuTimesKey) = 0;

    cpuTimes_.at(cpuTimesKey) += 
        qmomSubroutineRuntime<int>(callableSubroutine, nMoments);

    return 0;
}


int QmomProfiler::compute(double *moments)
{

    int info = 0;
    int nMoments = qmomPtr_->numberOfMoments();

    // First subroutine: Computation of the Jacobi matrix from moments
    info += computeJacobiMatrix(moments, true);

    // Second subroutine: Computation of the quadrature from the Jacobi matrix
    info += computeQuadrature(moments, true);

    // Third subroutine: Computation of the moments' rate of change from quadrature and given function
    info += computeMomentsRateOfChange(nMoments, true);

    return info;

}