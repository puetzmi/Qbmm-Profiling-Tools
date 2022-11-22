/**
 * @file physical_models.hpp
 * @author M. Puetz
 * @brief Some examples for the problem-dependent function in the moment transport equations,
 * more information is given in the `PhysicalModel` class description.
 *
 * @date 2022-11-21
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef PHYSICAL_MODELS_H
#define PHYSICAL_MODELS_H

#include <cmath>

#include "global_defs.hpp"
#include "factory.hpp"


/// @brief Default model constants in one place for convenience
namespace defaultModelConstants {

    // Constants in advection and diffusion terms
    double phi = 16;
    double gamma = 330;

}

class PhysicalModel;         // forward declaration

/// Factory for `PhysicalModel` objects
using PhysicalModelFactory = Factory<PhysicalModel>;


/**
 * @brief Abstract base class for physical models.
 * 
 * This class provides the interface for the computation of physical source terms. More precisely,
 * considering the moment equation
 * @f\[
 *  \frac{\mathrm{d}m_k}{\mathrm{d} t} = \int_{\Omega} n(\xi) \xi^k g(\xi) \mathrm{d} \xi
 *      \quad k=0,1,\dots,M-1
 * @f\]
 * and the corresponding quadrature approximation
 * @f\[
 *  \frac{\mathrm{d}m_k}{\mathrm{d} t} \approx \sum\limits_{j=0} w_j \xi_j^k g(\xi_j), 
 *      \quad k=0,1,\dots,M-1
 * @f\]
 * derived classes compute @f$ \xi_j^k g(\xi_j), \forall k \in \{0,1,\dots,M-1 \} @f$ given a single
 * value @f$ x_j @f$ and an arbitrary number of moments @f$ M @f$.
 * 
 */
class PhysicalModel
{

public:

    /**
     * @brief Construct a new `PhysicalModel` object.
     * 
     */
    PhysicalModel();

    /**
     * @brief Destroy the `PhysicalModel` object.
     * 
     */
    virtual ~PhysicalModel();


    /**
     * @brief Create new `PhysicalModel` object of specified type (using factory) and return unique pointer to it.
     * 
     * @param typeName Type name corresponding to key in the static map in `PhysicalModelFactory`.
     * @return std::unique_ptr<PhysicalModel> Unique pointer to new `PhysicalModel` object.
     */
    static std::unique_ptr<PhysicalModel> makeUnique(const std::string &typeName, unsigned int size);
    

    /**
     * @brief Create new `PhysicalModel` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the static map in `PhysicalModelFactory`.
     * @return std::shared_ptr<PhysicalModel> Shared pointer to new `PhysicalModel` object.
     */
    static std::shared_ptr<PhysicalModel> makeShared(const std::string &typeName, unsigned int size);


    /**
     * @brief Compute the moments' rate of change given a single node value and the
     *  number of moments and adds the result to given array.
     * 
     * @param [in] nodeValue Internal coordinate value at single node.
     * @param [in] nMoments Number of moments.
     * @param [out] momentsRateOfChange Moments' rate of change with the computed results added.
     * @return int Dummy error flag.
     */
    virtual int computeMomentsRateOfChange(double nodeValue, int nMoments,
        double *momentsRateOfChange) const = 0;


    /**
     * @brief Alternative way to call `computeMomentsRateOfChange`.
     * 
     */
    int operator()(double v, int nMoments, double *momentsRateOfChange) const
    {
        return this->computeMomentsRateOfChange(v, nMoments, momentsRateOfChange);
    }

};


/**
 * @brief Class to compute moment source terms based on the one dimensional Fokker-Planck equation
 * as described in @cite Puetz2022. 
 * 
 * The moments' rate of change is computed from the quadrature rule
 * @f\[
 *  \frac{\mathrm{d}m_k}{\mathrm{d} t} \approx \sum\limits_{j=0} w_j \xi_j^k g(\xi_j), 
 *      \quad k=0,1,\dots,M-1
 * @f\]
 * where
 * @f\[
 *  \xi^k g(\xi) = k \text{sgn}(\xi) \left[ -\gamma \xi^{k+1} + \frac{\phi^2}{4} \xi^{k-1}
 *      + (k-1) \frac{\phi^2}{2} \xi^{k-1} \right],
 * @f\]
 * see @cite Puetz2022, Eq. (34).
 * 
 * 
 */
class FokkerPlanckEquation : public PhysicalModel 
{

public:

    /**
     * @brief Construct a new `FokkerPlanckEquation` object.
     * 
     */
    FokkerPlanckEquation();

    /**
     * @brief Construct a new `FokkerPlanckEquation` object.
     * 
     * @param gamma Constant in the advection terms (see `FokkerPlanckEquation` class description).
     * @param phi Constant in the diffusion terms (see `FokkerPlanckEquation` class description).
     */
    FokkerPlanckEquation(double gamma, double phi);

    /**
     * @brief Destroy the `FokkerPlanckEquation` object.
     * 
     */
    virtual ~FokkerPlanckEquation();


    virtual int computeMomentsRateOfChange(double v, int nMoments, 
        double *momentsRateOfChange) const;


private:

    double gamma_;              ///< Constant in advection term
    double phi_;                ///< Constant in diffusion term
    double phiSquared_;         ///< Squared constant in diffusion term

    static bool inline registered_ = REGISTER_TYPE(PhysicalModelFactory, FokkerPlanckEquation);

};

#endif // PHYSICAL_MODELS_H