/**
 * @file qmom_profiler.hpp
 * @author M. Puetz
 * @brief Class for performance profiling of different implementations of the quadrautre method of moments and derived methods.
 * @date 2022-10-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef QMOM_PROFILER_HPP
#define QMOM_PROFILER_HPP

#include "global_defs.hpp"
#include "qmom.hpp"
#include <vector>
#include <string>
#include <map>
#include <chrono>
#include <memory>
#include <functional>


/**
 * @brief Class serving as a wrapper for QMOM object with additional attributes and methods for performance profiling.
 * 
 */
class QmomProfiler
{

private:

    static inline const std::vector<std::string> subroutineNames_ = 
    {
            "CoreInversion",
            "Quadrature",
            "Integrate"
    };                                          ///< Names of subroutines.

    int nExecutions_;                           ///< The number of times subroutines are executed to compute an average CPU time.
    std::map<std::string, double> cpuTimes_;    ///< Map holding the measured CPU times per subroutine.
    std::chrono::high_resolution_clock clock_;  ///< Clock to measure CPU times.
    std::shared_ptr<Qmom> qmomPtr_;             ///< Pointer to QMOM object.


    /**
     * @brief Get constant pointer to `Qmom` object.
     * 
     * @return std::shared_ptr<const Qmom> Read-only shared pointer to `Qmom` object.
     */
    std::shared_ptr<const Qmom> qmomPtr() const
    {
        return qmomPtr_;
    }


    /**
     * @brief Template to measure runtimes of any of the QMOM subroutines
     * 
     * @tparam ArgTypes Argument types of the QMOM subroutine function called.
     * @param callableSubroutine QMOM subroutine with the argument types `ArgTypes`.
     * @param args Arguments passed to callable subroutine.
     * @return double Average runtime.
     */
    template <class... ArgTypes>
    inline double qmomSubroutineRuntime(std::function<int(ArgTypes...)> callableSubroutine,
        ArgTypes... args);


public:

    /**
     * @brief Construct a new `QmomProfiler` object.
     * 
     * @param numberOfExecutions Number of executions of subroutines to compute average CPU time.
     * @param qmomPtr Shared Pointer to QMOM object.
     */
    QmomProfiler(int numberOfExecutions, std::shared_ptr<Qmom> qmomPtr);

    /**
     * @brief Destroy the `QmomProfiler` object.
     * 
     */
    virtual ~QmomProfiler();

    /**
     * @brief Get read-only reference to `cpuTimes_` map;
     * 
     * @return const std::map<std::string, double>& Map holding the CPU times of each subroutine.
     */
    const std::map<std::string, double>& cpuTimes() {
        return cpuTimes_;
    }

    /**
     * @brief Get subroutine names.
     * 
     * @return const std::vector<std::string>& Reference (read-only) to vector of subroutine names.
     */
    static const std::vector<std::string>& subroutineNames() {
        return subroutineNames_;
    }


    /**
     * @brief Execute all QMOM subroutines or derived algorithms and measure runtimes.
     * 
     * @param moments Realizable moment sequence.
     * @return int Error flag.
     */
    virtual int compute(double *moments);


    /**
     * @brief Interface for `Qmom::computeJacobiMatrix` with runtime measurement.
     * 
     * @param moments A realizable moment sequence.
     * @param resetCpuTime Indicates if stored CPU times are supposed to be reset (optional, default `true`).
     * @return int Error flag.
     */
    virtual int computeJacobiMatrix(double *moments, bool resetCpuTime=true);


    /**
     * @brief Interface method for `Qmom::computeQuadrature` with runtime measurement.
     * 
     * @param moments A realizable moment sequence.
     * @param resetCpuTime Indicates if stored CPU times are supposed to be reset (optional, default `true`).
     * @return int Error flag
     */
    virtual int computeQuadrature(double *moments, bool resetCpuTime);


    /**
     * @brief Interface method for `Qmom::computeMomentsRateOfChange` with runtime measurement.
     * @param nMoments Number of moments to be computed.
     * @param resetCpuTime Indicates if stored CPU times are supposed to be reset (optional, default `true`).
     * @return int Dummy error flag.
     */
    virtual int computeMomentsRateOfChange(int nMoments, bool resetCpuTime=true);

};

#endif // QMOM_PROFILER_HPP