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
    constexpr double phi = 16;
    constexpr double gamma = 330;

    // Constants in the computation of Reynolds-dependent
    // drag coefficient: 
    // Cd = 24/Re*(1 + preReynoldsFactor*Re^reynoldsExponent) 
    //      + dragConstant
    constexpr double dragPreReynoldsFactor = 1/6;
    constexpr double dragReynoldsExponent = 2/3;
    constexpr double dragConstant = 0;

    // Modified constants for velocity-dependent-drag FPE
    // using the constants that `gamma` and `phi` are
    // based on [Puetz2022]
    constexpr double fluidDensity = 1.;
    constexpr double turbulentKineticEnergy = 1.;
    constexpr double particleDiameter = 1e-5;

    // Further constants needed for the calculation of drag coefficient
    constexpr double particleDensity = 100*fluidDensity;
    constexpr double dynamicViscosity = 2e-5;


    // Constants in collision terms
    constexpr double coefficientOfRestitution = 0.8;
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
    static std::unique_ptr<PhysicalModel> makeUnique(const std::string &typeName);
    

    /**
     * @brief Create new `PhysicalModel` object of specified type (using factory) and return shared pointer to it.
     * 
     * @param typeName Type name corresponding to key in the static map in `PhysicalModelFactory`.
     * @return std::shared_ptr<PhysicalModel> Shared pointer to new `PhysicalModel` object.
     */
    static std::shared_ptr<PhysicalModel> makeShared(const std::string &typeName);


    /**
     * @brief Compute the rate of change of a specified number of moments given a set of abscissas
     *  and weights
     * 
     * @param abscissas Array of abscissas.
     * @param weights Array of weights corresponding to abscissas.
     * @param nNodes Total number of nodes.
     * @param nMoments Total number of moments.
     * @param momentsRateOfChange Allocated memory for `nMoments` doubles to store the results, which are added to the existing values.
     * @return int Error flag.
     */
    virtual int computeMomentsRateOfChange(double * const abscissas, double * const weights, int nNodes,
        int nMoments, double *momentsRateOfChange) const = 0;


    /**
     * @brief Alternative way to call `computeMomentsRateOfChange`.
     * 
     */
    int operator()(double * const abscissas, double * const weights, int nNodes,
        int nMoments, double *momentsRateOfChange) const
    {
        return this->computeMomentsRateOfChange(abscissas, weights, nNodes, nMoments, momentsRateOfChange);
    }

};


/**
 * @brief Class to compute moment source terms based on the one dimensional Fokker-Planck equation
 * using a constant drag coefficient as described in @cite Puetz2022. 
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
class FokkerPlanckConstantCd : public PhysicalModel 
{

public:

    /**
     * @brief Construct a new `FokkerPlanckConstantCd` object.
     * 
     */
    FokkerPlanckConstantCd();

    /**
     * @brief Construct a new `FokkerPlanckConstantCd` object.
     * 
     * @param gamma Constant in the advection terms (see `FokkerPlanckConstantCd` class description).
     * @param phi Constant in the diffusion terms (see `FokkerPlanckConstantCd` class description).
     */
    FokkerPlanckConstantCd(double gamma, double phi);

    /**
     * @brief Destroy the `FokkerPlanckConstantCd` object.
     * 
     */
    virtual ~FokkerPlanckConstantCd();


    virtual int computeMomentsRateOfChange(double * const abscissas, double * const weights, int nNodes,
        int nMoments, double *momentsRateOfChange) const;


private:

    double gamma_;              ///< Constant in advection term
    double phi_;                ///< Constant in diffusion term

    static bool inline registered_ = REGISTER_TYPE(PhysicalModelFactory, FokkerPlanckConstantCd);

};


/**
 * @brief Class to compute moment source terms based on the one dimensional
 * Fokker-Planck equation as described in @cite Puetz2022, extended by a
 * Reynolds-number/velocity-dependent drag coefficient.
 *
 */
class FokkerPlanckVelocityDependentCd : public PhysicalModel 
{

public:

    /**
     * @brief Construct a new `FokkerPlanckVelocityDependentCd` object.
     * 
     */
    FokkerPlanckVelocityDependentCd();

    /**
     * @brief Construct a new `FokkerPlanckVelocityDependentCd` object.
     * 
     * @param dragPreReynoldsFactor Pre-factor in Reynolds-number-dependent drag coefficient formula.
     * @param dragReynoldsExponent Exponent in Reynolds-number-dependent drag coefficient formula.
     * @param dragConstant Constant in Reynolds-number-dependent drag coefficient formula.
     * @param fluidDensity Density of the carrier fluid.
     * @param particleDensity Uniform particle density.
     * @param turbulentKineticEnergy Constant turbulent kinetic energy.
     * @param dynamicViscosity Dynamic viscosity of the carrier fluid.
     * @param particleDiameter Uniform particle diameter.
     */
    FokkerPlanckVelocityDependentCd(
            double dragPreReynoldsFactor,
            double dragReynoldsExponent,
            double dragConstant,
            double fluidDensity,
            double particleDensity,
            double turbulentKineticEnergy,
            double dynamicViscosity,
            double particleDiameter
    );

    /**
     * @brief Destroy the `FokkerPlanckVelocityDependentCd` object.
     * 
     */
    virtual ~FokkerPlanckVelocityDependentCd();


    virtual int computeMomentsRateOfChange(double * const abscissas, double * const weights, int nNodes,
        int nMoments, double *momentsRateOfChange) const;


private:

    // Constants in the computation of Reynolds-dependent
    // drag coefficient: 
    // Cd = 24/Re*(1 + preReynoldsFactor*Re^reynoldsExponent) 
    //      + dragConstant
    const double dragPreReynoldsFactor_;        ///< Pre-factor in Reynolds-number-dependent drag coefficient formula
    const double dragReynoldsExponent_;         ///< Exponent in Reynolds-number-dependent drag coefficient formula
    const double dragConstant_;                 ///< Constant in Reynolds-number-dependent drag coefficient formula

    // Additional material constants used for the computation
    // of the drag coefficient
    const double fluidDensity_;                 ///< Density of the carrier fluid
    const double particleDensity_;              ///< Uniform particle density
    const double turbulentKineticEnergy_;       ///< Constant turbulent kinetic energy
    const double dynamicViscosity_;             ///< Dynamic viscosity of the carrier fluid
    const double particleDiameter_;             ///< Uniform particle diameter


    static bool inline registered_ = REGISTER_TYPE(PhysicalModelFactory, FokkerPlanckVelocityDependentCd);

};


/**
 * @brief Hard-sphere collision as described in @cite Fox2010 and @cite Marchisio2013, reduced to one dimension.
 * 
 */
class HardSphereCollision1D : public PhysicalModel
{

public:

    /**
     * @brief Construct a new `HardSphereCollision1D` object.
     * 
     */
    HardSphereCollision1D();

    /**
     * @brief Construct a new `HardSphereCollision1D` object.
     * 
     * @param coefficientOfRestitution Coefficient of restitiution, see @cite FOx2010.
     * @param particleDiameter Uniform particle diameter.
     * @param nMomentsMax Maximum number of moments.
     */
    HardSphereCollision1D
    (
        double coefficientOfRestitution,
        double particleDiameter,
        int nMomentsMax
    );

    /**
     * @brief Destroy the `HardSphereCollision1D` object.
     * 
     */
    ~HardSphereCollision1D();


    virtual int computeMomentsRateOfChange(double * const abscissas, double * const weights, int nNodes,
        int nMoments, double *momentsRateOfChange) const;


private:

    /**
     * @brief Compute the expression $I^{(0)}_{k,0,0}$ in $k$th order collision
     * integral, see @cite Fox2010.
     *
     * @tparam momentOrder Moment order (one-dimensional).
     * @return double Value of $I^{(0)}_{k,0,0}$, where $k$ is the given moment
     * order.
     */
    template <int momentOrder>
    double computeI0();
    
    /**
     * @brief Compute the expression $I^{(1)}_{k,0,0}$ in $k$th order collision
     * integral, see @cite Fox2010.
     *
     * @tparam momentOrder Moment order (one-dimensional).
     * @return double Value of $I^{(1)}_{k,0,0}$, where $k$ is the given moment
     * order.
     */
    template <int momentOrder>
    double computeI1();

    /**
     * @brief Initialize vector of functions to compute the expression
     * $I_{k,0,0}^(0)$ in the collision integral, see @cite Fox2010.
     *
     * @param nMoments Number of moments.
     * @return std::vector<std::function<double()> > Vector of functions.
     */
    std::vector<std::function<double()> >
        initializeI0Functions(int nMoments);

    /**
     * @brief Initialize vector of functions to compute the expression
     * $I_{k,0,0}^(1)$ in the collision integral, see @cite Fox2010.
     *
     * @param nMoments Number of moments.
     * @return std::vector<std::function<double()> > Vector of functions.
     */
    std::vector<std::function<double()> >
        initializeI1Functions(int nMoments);

    double coefficientOfRestitution_;               ///< Coefficient of restitution (measure of elasticity), must be in [0, 1].
    double particleDiameter_;                       ///< Uniform particle diameter.
    std::vector<std::function<double()> > 
        I0Functions_;                               ///< Functions to compute I0 from particle velocities (see @cite Fox2010).
    std::vector<std::function<double()> > 
        I1Functions_;                               ///< Functions to compute I1 from particle velocities (see @cite Fox2010).
    double omega_;                                  ///< Measure of elasticity in I0 functions, see @cite Fox2010.
    double *omegaPower_;                            ///< Array of powers of `omega_`.
    double *g1Power_;                               ///< Array of powers of relative velocity (1D -> here only first component).
    double *gPower_;                                ///< Array of powers of relative velocity magnitude (1D -> equals magnitude of first component).
    double *v1Power_;                               ///< Array of powers of particle velocity (1D -> here only first component).

    
    static constexpr int nMomentsMax_ = 20;         ///< Maximum allowed number of moments.

    static bool inline registered_ = REGISTER_TYPE(PhysicalModelFactory, HardSphereCollision1D);

};

#endif // PHYSICAL_MODELS_H