/**
 * @file physical_models.cpp
 * @author M. Puetz
 * @brief Implementation of the `PhysicalModel` class and derived classes.
 * @date 2022-11-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include <stdexcept>
#include "physical_models.hpp"


// Sign function for convenience
template<class T=double>
inline int signFunction(T x, double tol=1e-15) {
    return static_cast<int>((x - tol > 0) - (x + tol < 0));
}

PhysicalModel::PhysicalModel()
{
}

PhysicalModel::~PhysicalModel()
{
}


FokkerPlanckConstantCd::FokkerPlanckConstantCd(double gamma, double phi)
    :
    PhysicalModel(),
    gamma_(gamma),
    phi_(phi)
{

    auto errorMessage = [](const std::string&& s){
        return std::string("The parameter `" + s + "` must be non-negative.");
    };

    if (gamma_ < 0) {
        throw std::runtime_error(errorMessage("gamma"));
    }

    if (phi_ < 0) {
        throw std::runtime_error(errorMessage("phi"));
    }
}


FokkerPlanckConstantCd::FokkerPlanckConstantCd()
    :
    FokkerPlanckConstantCd(
        defaultModelConstants::gamma, defaultModelConstants::phi
    )
{
}


FokkerPlanckConstantCd::~FokkerPlanckConstantCd()
{
}
    

int FokkerPlanckConstantCd::computeMomentsRateOfChange(double * const abscissas,
    double * const weights, int nNodes, int nMoments,
    double *momentsRateOfChange) const
{
    for (int j=0; j<nNodes; j++) {

        double abscissa = abscissas[j];

        // Sign of internal coordinate value at node
        int sign = signFunction(abscissa);
        if (sign == 0) {
            continue;       // no addition to moment source term for this node
        }

        // advection terms
        double dmdt = weights[j] * (-gamma_*sign*abscissa*abscissa + 0.25*phi_*phi_*sign);
        for (int k=1; k<nMoments; k++) {
            momentsRateOfChange[k] += k*dmdt;
            dmdt *= abscissa;
        }

        // diffusion term
        dmdt = weights[j] * (0.5*phi_*phi_*sign*abscissa);
        for (int k=2; k<nMoments; k++) {
            momentsRateOfChange[k] += k*(k-1)*dmdt;
            dmdt *= abscissa;
        }
    }

    return 0;
}


FokkerPlanckVelocityDependentCd::FokkerPlanckVelocityDependentCd
    (
        double dragPreReynoldsFactor,
        double dragReynoldsExponent,
        double dragConstant,
        double fluidDensity,
        double particleDensity,
        double turbulentKineticEnergy,
        double dynamicViscosity,
        double particleDiameter
    )
    :
    PhysicalModel(),
    dragPreReynoldsFactor_(dragPreReynoldsFactor),
    dragReynoldsExponent_(dragReynoldsExponent),
    dragConstant_(dragConstant),
    fluidDensity_(fluidDensity),
    particleDensity_(particleDensity),
    turbulentKineticEnergy_(turbulentKineticEnergy),
    dynamicViscosity_(dynamicViscosity),
    particleDiameter_(particleDiameter)
{

    auto errorMessage = [](const std::string&& s){
        return std::string("The parameter `" + s + "` must be non-negative.");
    };

    if (dragPreReynoldsFactor_ < 0) {
        throw std::runtime_error(errorMessage("dragPreReynoldsFactor"));
    }

    if (dragReynoldsExponent_ < 0) {
        throw std::runtime_error(errorMessage("dragReynoldsExponent"));
    }

    if (dragConstant_ < 0) {
        throw std::runtime_error(errorMessage("dragConstant"));
    }

    if (fluidDensity_ < 0) {
        throw std::runtime_error(errorMessage("fluidDensity"));
    }

    if (particleDensity_ < 0) {
        throw std::runtime_error(errorMessage("particleDensity"));
    }

    if (turbulentKineticEnergy_ < 0) {
        throw std::runtime_error(errorMessage("turbulentKineticEnergy"));
    }

    if (dynamicViscosity_ < 0) {
        throw std::runtime_error(errorMessage("dynamicViscosity"));
    }

    if (particleDiameter_ < 0) {
        throw std::runtime_error(errorMessage("particleDiameter"));
    }
}



FokkerPlanckVelocityDependentCd::FokkerPlanckVelocityDependentCd()
    :
    FokkerPlanckVelocityDependentCd(
        defaultModelConstants::dragPreReynoldsFactor,
        defaultModelConstants::dragReynoldsExponent,
        defaultModelConstants::dragConstant,
        defaultModelConstants::fluidDensity,
        defaultModelConstants::particleDensity,
        defaultModelConstants::turbulentKineticEnergy,
        defaultModelConstants::dynamicViscosity,
        defaultModelConstants::particleDiameter
    )
{
}


FokkerPlanckVelocityDependentCd::~FokkerPlanckVelocityDependentCd()
{
}
    

int FokkerPlanckVelocityDependentCd::computeMomentsRateOfChange(double * const abscissas,
    double * const weights, int nNodes, int nMoments,
    double *momentsRateOfChange) const
{
    for (int j=0; j<nNodes; j++) {

        const double abscissa = abscissas[j];

        double dragCoefficient = computeDragCoefficient(abscissa);
        double tauStarI = 0.75*fluidDensity_*dragCoefficient
            / (particleDensity_*particleDiameter_);

        // advection terms
        //  TODO:
        double driftCoefficient = tauStarI*abscissa;
        double dmdt = weights[j]*driftCoefficient;
        for (int k=1; k<nMoments; k++) {
            momentsRateOfChange[k] += k*dmdt;
            dmdt *= abscissa;
        }

        // diffusion term
        //  TODO:
        double diffusionCoefficient = turbulentKineticEnergy_*tauStarI*abscissa;
        dmdt = weights[j]*diffusionCoefficient;
        for (int k=2; k<nMoments; k++) {
            momentsRateOfChange[k] += k*(k-1)*dmdt;
            dmdt *= abscissa;
        }
    }

    return 0;
}