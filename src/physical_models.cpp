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


FokkerPlanckEquation::FokkerPlanckEquation(double gamma, double phi)
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


FokkerPlanckEquation::FokkerPlanckEquation()
    :
    FokkerPlanckEquation(
        defaultModelConstants::gamma, defaultModelConstants::phi
    )
{
}


FokkerPlanckEquation::~FokkerPlanckEquation()
{
}
    

int FokkerPlanckEquation::computeMomentsRateOfChange(double * const abscissas,
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