/**
 * @file physical_models.hpp
 * @author M. Puetz
 * @brief Some examples for the problem-dependent function in the moment transport equations.
 * 
 * @f$ \mathrm{d}m_k / \mathrm{d} t @f$ = \int_{\Omega} \xi^k g(\xi) n(\xi) \mathrm{d} \xi @f$.
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


class PhysicalModel
{

public:

    PhysicalModel();
    virtual ~PhysicalModel();

    static std::unique_ptr<PhysicalModel> makeUnique(const std::string &typeName, unsigned int size);
    static std::shared_ptr<PhysicalModel> makeShared(const std::string &typeName, unsigned int size);

    virtual int computeMomentsRateOfChange(double v, int nMoments,
        double *momentsRateOfChange) const = 0;

    int operator()(double v, int nMoments, double *momentsRateOfChange) const
    {
        return this->computeMomentsRateOfChange(v, nMoments, momentsRateOfChange);
    }

};


class FokkerPlanckEquation : public PhysicalModel 
{

public:

    FokkerPlanckEquation();
    FokkerPlanckEquation(double gamma, double phi);
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