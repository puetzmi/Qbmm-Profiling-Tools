/**
 * @file physical_models.cpp
 * @author M. Puetz
 * @brief Implementation of the `PhysicalModel` class and derived classes.
 * @date 2022-11-22
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#include <functional>
#include <vector>
#include <stdexcept>
#include <mkl.h>
#include "collision_integral_functions.hpp"
#include "physical_models.hpp"
#include "constants.hpp"


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

std::shared_ptr<PhysicalModel> PhysicalModel::makeShared(const std::string &typeName)
{
    return PhysicalModelFactory::makeShared(typeName);
}

std::unique_ptr<PhysicalModel> PhysicalModel::makeUnique(const std::string &typeName)
{
    return PhysicalModelFactory::makeUnique(typeName);
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

        const double velocity = abscissas[j];

        int signVelocity = signFunction(velocity);
        if (signVelocity == 0)
            continue;

        double absVelocity = signVelocity*velocity;

        double reynoldsNumber = absVelocity*particleDiameter_*fluidDensity_/dynamicViscosity_;
        double reynoldsFactor = dragPreReynoldsFactor_* std::pow(reynoldsNumber, dragReynoldsExponent_);

        // pre-factor of both advection and diffusion terms
        double commonPreFactor = 18*dynamicViscosity_
                                / (particleDiameter_*particleDiameter_*particleDensity_*velocity);

        // contribution from advection terms (drift and spurious drift)
        double advectionCoefficient = commonPreFactor
                                    * (-velocity*velocity*(1 + reynoldsFactor)
                                        + 0.25*reynoldsFactor*dragReynoldsExponent_*turbulentKineticEnergy_);
        double dmdt = weights[j]*advectionCoefficient;
        for (int k=1; k<nMoments; k++) {
            momentsRateOfChange[k] += k*dmdt;
            dmdt *= velocity;
        }

        // contribution from diffusion term
        double diffusionCoefficient = 0.5*commonPreFactor*velocity*turbulentKineticEnergy_*(1 + reynoldsFactor);
        dmdt = weights[j]*diffusionCoefficient;
        for (int k=2; k<nMoments; k++) {
            momentsRateOfChange[k] += k*(k-1)*dmdt;
            dmdt *= velocity;
        }
    }

    return 0;
}


HardSphereCollision1D::HardSphereCollision1D()
    :
    HardSphereCollision1D
    (
        defaultModelConstants::coefficientOfRestitution,
        defaultModelConstants::particleDiameter,
        nMomentsMax_
    )
{
}


HardSphereCollision1D::HardSphereCollision1D(double coefficientOfRestitution,
    double particleDiameter, int nMomentsMax)
    :
    coefficientOfRestitution_(coefficientOfRestitution),
    particleDiameter_(particleDiameter),
    I0Functions_(initializeI0Functions(nMomentsMax)),
    I1Functions_(initializeI1Functions(nMomentsMax)),
    omega_(0.5*(coefficientOfRestitution_ + 1)), // Fox2010, Eq. (28)
    omegaPower_(static_cast<double*>(mkl_malloc(nMomentsMax*sizeof(double), MALLOC_ALIGN))),
    g1Power_(static_cast<double*>(mkl_malloc(nMomentsMax*sizeof(double), MALLOC_ALIGN))),
    gPower_(static_cast<double*>(mkl_malloc(nMomentsMax*sizeof(double), MALLOC_ALIGN))),
    v1Power_(static_cast<double*>(mkl_malloc(nMomentsMax*sizeof(double), MALLOC_ALIGN)))
{

    if (particleDiameter_ <= 0)
        throw std::runtime_error(
            "The paramter `particleDiameter` must be positive."
        );
    
    if (coefficientOfRestitution < 0 or coefficientOfRestitution_ > 1)
        throw std::runtime_error(
            "The paramter `coefficientOfRestitution` must be within [0,1]."
        );

    omegaPower_[0] = 1;
    for (int k=1; k<nMomentsMax; k++) {
        omegaPower_[k] = omegaPower_[k-1]*omega_;
    }

    g1Power_[0] = 1;
    gPower_[0] = 1;
    v1Power_[0] = 1;
}


HardSphereCollision1D::~HardSphereCollision1D()
{
    mkl_free(omegaPower_);
    mkl_free(g1Power_);
    mkl_free(gPower_);
    mkl_free(v1Power_);
}


int HardSphereCollision1D::computeMomentsRateOfChange
(
    double * const abscissas, double * const weights, 
    int nNodes, int nMoments, double *momentsRateOfChange
) const
{

    double m0 = 0;
    for (int i=0; i< nNodes; i++) {
        m0 += weights[i];
    }

    double volumeFraction = m0*constants::pi/6;
    for (int i=0; i<3; i++) {
        volumeFraction *= particleDiameter_;
    }

    // see Marchisio2013, Chapter 6.1
    double c = std::min<double>(volumeFraction/0.63, 1);
    double g0 = (2 - c) / (2*(1 - c)*(1 - c)*(1 - c));


    // Fox2010, Eq. (59) without advection and collisional flux term
    for (int i=0; i<nNodes; i++) {

        // update powers of velocity
        v1Power_[1] = abscissas[i];
        for (int k=2; k<nMoments; k++) {
            v1Power_[k] = v1Power_[k-1]*v1Power_[1];
        }

        for (int j=0; j<nNodes; j++) {

            // update powers of relative velocity and its magnitude
            g1Power_[1] = v1Power_[1] - abscissas[j];
            gPower_[1] = std::abs(g1Power_[1]);
            for (int k=2; k<nMoments; k++) {
                g1Power_[k] = g1Power_[k-1]*g1Power_[1];
                gPower_[k] = gPower_[k-1]*gPower_[1];
            }

            // compute moment source terms
            for (int k=1; k<nMoments; k++) {
                momentsRateOfChange[k] += g0*(
                    6/particleDiameter_*weights[i]*weights[j]*gPower_[1]*I0Functions_[k]()
                );
            }
        }
    }

    return 0;
}