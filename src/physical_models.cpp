#include "physical_models.hpp"
#include <stdexcept>


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
    phi_(phi),
    phiSquared_(phi*phi)
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
    

int FokkerPlanckEquation::computeMomentsRateOfChange(double nodeValue, int nMoments, 
    double *momentsRateOfChange) const
{

        // Sign of internal coordinate value at node
        int sign = (nodeValue > 0) - (nodeValue < 0);

        // No change if value is 'numerically zero'
        if (sign*nodeValue < std::numeric_limits<double>::epsilon()) {
            return 0;
        }

        // advection terms
        double dmdt = -gamma_*sign*nodeValue + 0.25*phiSquared_*sign;
        for (int k=1; k<nMoments; k++) {
            momentsRateOfChange[k] += k*dmdt;
            dmdt *= nodeValue;
        }

        // diffusion term
        dmdt = 0.5*phiSquared_*sign*nodeValue;
        for (int k=2; k<nMoments; k++) {
            momentsRateOfChange[k] += k*(k-1)*dmdt;
            dmdt *= nodeValue;
        }

        return 0;
}