/**
 * @file collision_integral_functions.hpp
 * @author M. Puetz
 * @brief This file contains the analytical expressions for I(0) in the
 * collision integrals (see @cite Fox2010 and @cite Marchisio2013) in one
 * spatial dimension up to 19th order moments.
 * @date 2022-12-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#ifndef COLLISION_INTEGRAL_FUNCTIONS_HPP
#define COLLISION_INTEGRAL_FUNCTIONS_HPP

#include "physical_models.hpp"


template<>
double HardSphereCollision1D::computeI0<0>()
{
    return 0;
}


template<>
double HardSphereCollision1D::computeI0<1>()
{
    return -g1Power_[1]*omegaPower_[1]/2;
}


template<>
double HardSphereCollision1D::computeI0<2>()
{
    return omegaPower_[2]*(3*g1Power_[2] + gPower_[2])/12 -
        g1Power_[1]*omegaPower_[1]*v1Power_[1];
}


template<>
double HardSphereCollision1D::computeI0<3>()
{
    return -omegaPower_[3]*(g1Power_[1]*gPower_[2] + g1Power_[3])/8 +
        omegaPower_[2]*(3*g1Power_[2] + gPower_[2])/4*v1Power_[1] -
        3*g1Power_[1]*omegaPower_[1]/2*v1Power_[2];
}


template<>
double HardSphereCollision1D::computeI0<4>()
{
    return omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/80 - omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/2*v1Power_[1] + omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/2*v1Power_[2] - 2*g1Power_[1]*omegaPower_[1]*v1Power_[3];
}


template<>
double HardSphereCollision1D::computeI0<5>()
{
    return -omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/96 +
        omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[1] - 5*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/4*v1Power_[2] + 5*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/6*v1Power_[3] - 5*g1Power_[1]*omegaPower_[1]/2*v1Power_[4];
}


template<>
double HardSphereCollision1D::computeI0<6>()
{
    return omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/448 -
        omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2] +
        3*g1Power_[5])/16*v1Power_[1] +
        3*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[2] - 5*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/2*v1Power_[3] + 5*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[4] - 3*g1Power_[1]*omegaPower_[1]*v1Power_[5];
}


template<>
double HardSphereCollision1D::computeI0<7>()
{
    return -omegaPower_[7]*(g1Power_[1]*gPower_[6] +
        7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] + g1Power_[7])/128 +
        omegaPower_[6]*(21*g1Power_[2]*gPower_[4] + 35*g1Power_[4]*gPower_[2] +
        7*g1Power_[6] + gPower_[6])/64*v1Power_[1] -
        7*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2] +
        3*g1Power_[5])/32*v1Power_[2] +
        7*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[3] - 35*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/8*v1Power_[4] + 7*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[5] - 7*g1Power_[1]*omegaPower_[1]/2*v1Power_[6];
}


template<>
double HardSphereCollision1D::computeI0<8>()
{
    return omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/2304 - omegaPower_[7]*(g1Power_[1]*gPower_[6] +
        7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/16*v1Power_[1] + omegaPower_[6]*(21*g1Power_[2]*gPower_[4]
        + 35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/16*v1Power_[2]
        - 7*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/12*v1Power_[3] +
        7*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/8*v1Power_[4] - 7*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])*v1Power_[5] + 7*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/3*v1Power_[6] - 4*g1Power_[1]*omegaPower_[1]*v1Power_[7];
}


template<>
double HardSphereCollision1D::computeI0<9>()
{
    return -omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/2560 +
        omegaPower_[8]*(36*g1Power_[2]*gPower_[6] + 126*g1Power_[4]*gPower_[4] +
        84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] + gPower_[8])/256*v1Power_[1]
        - 9*omegaPower_[7]*(g1Power_[1]*gPower_[6] + 7*g1Power_[3]*gPower_[4] +
        7*g1Power_[5]*gPower_[2] + g1Power_[7])/32*v1Power_[2] +
        3*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] + 35*g1Power_[4]*gPower_[2]
        + 7*g1Power_[6] + gPower_[6])/16*v1Power_[3] -
        21*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/16*v1Power_[4] +
        63*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/40*v1Power_[5] - 21*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/2*v1Power_[6] + 3*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])*v1Power_[7] - 9*g1Power_[1]*omegaPower_[1]/2*v1Power_[8];
}


template<>
double HardSphereCollision1D::computeI0<10>()
{
    return omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8]
        + 330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/11264 -
        omegaPower_[9]*(5*g1Power_[1]*gPower_[8] + 60*g1Power_[3]*gPower_[6] +
        126*g1Power_[5]*gPower_[4] + 60*g1Power_[7]*gPower_[2] +
        5*g1Power_[9])/256*v1Power_[1] +
        5*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] + 126*g1Power_[4]*gPower_[4]
        + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/256*v1Power_[2] - 15*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/16*v1Power_[3] +
        15*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] + 35*g1Power_[4]*gPower_[2]
        + 7*g1Power_[6] + gPower_[6])/32*v1Power_[4] -
        21*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/8*v1Power_[5] +
        21*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/8*v1Power_[6] - 15*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])*v1Power_[7] + 15*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[8] - 5*g1Power_[1]*omegaPower_[1]*v1Power_[9];
}


template<>
double HardSphereCollision1D::computeI0<11>()
{
    return -omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11]
        + 55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/6144 +
        omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/1024*v1Power_[1] -
        11*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] + 60*g1Power_[3]*gPower_[6]
        + 126*g1Power_[5]*gPower_[4] + 60*g1Power_[7]*gPower_[2] +
        5*g1Power_[9])/512*v1Power_[2] +
        55*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/768*v1Power_[3] - 165*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/64*v1Power_[4] +
        33*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] + 35*g1Power_[4]*gPower_[2]
        + 7*g1Power_[6] + gPower_[6])/32*v1Power_[5] -
        77*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/16*v1Power_[6] +
        33*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/8*v1Power_[7] - 165*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/8*v1Power_[8] + 55*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/12*v1Power_[9] -
        11*g1Power_[1]*omegaPower_[1]/2*v1Power_[10];
}


template<>
double HardSphereCollision1D::computeI0<12>()
{
    return omegaPower_[12]*(286*g1Power_[10]*gPower_[2] +
        13*g1Power_[12] + 78*g1Power_[2]*gPower_[10] +
        715*g1Power_[4]*gPower_[8] + 1716*g1Power_[6]*gPower_[6] +
        1287*g1Power_[8]*gPower_[4] + gPower_[12])/53248 -
        omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/512*v1Power_[1]
        + 3*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/512*v1Power_[2] -
        11*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] + 60*g1Power_[3]*gPower_[6]
        + 126*g1Power_[5]*gPower_[4] + 60*g1Power_[7]*gPower_[2] +
        5*g1Power_[9])/128*v1Power_[3] +
        55*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/256*v1Power_[4] - 99*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/16*v1Power_[5] +
        33*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] + 35*g1Power_[4]*gPower_[2]
        + 7*g1Power_[6] + gPower_[6])/16*v1Power_[6] -
        33*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/4*v1Power_[7] +
        99*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[8] - 55*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])/2*v1Power_[9] + 11*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/2*v1Power_[10] - 6*g1Power_[1]*omegaPower_[1]*v1Power_[11];
}


template<>
double HardSphereCollision1D::computeI0<13>()
{
    return -omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] + 1001*g1Power_[9]*gPower_[4])/57344 +
        omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/4096*v1Power_[1] -
        13*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/1024*v1Power_[2]
        + 13*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/512*v1Power_[3] -
        143*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] + 60*g1Power_[3]*gPower_[6]
        + 126*g1Power_[5]*gPower_[4] + 60*g1Power_[7]*gPower_[2] +
        5*g1Power_[9])/512*v1Power_[4] +
        143*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/256*v1Power_[5] - 429*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/32*v1Power_[6] +
        429*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/112*v1Power_[7]
        - 429*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/32*v1Power_[8] +
        143*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[9] - 143*omegaPower_[3]*(g1Power_[1]*gPower_[2]
        + g1Power_[3])/4*v1Power_[10] + 13*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/2*v1Power_[11] -
        13*g1Power_[1]*omegaPower_[1]/2*v1Power_[12];
}


template<>
double HardSphereCollision1D::computeI0<14>()
{
    return omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/245760 - omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/4096*v1Power_[1] +
        7*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/4096*v1Power_[2] -
        91*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/1536*v1Power_[3]
        + 91*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/1024*v1Power_[4] -
        1001*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/1280*v1Power_[5] +
        1001*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/768*v1Power_[6] - 429*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/16*v1Power_[7] +
        429*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/64*v1Power_[8] -
        1001*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/48*v1Power_[9] +
        1001*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/80*v1Power_[10] - 91*omegaPower_[3]*(g1Power_[1]*gPower_[2]
        + g1Power_[3])/2*v1Power_[11] + 91*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/12*v1Power_[12] - 7*g1Power_[1]*omegaPower_[1]*v1Power_[13];
}


template<>
double HardSphereCollision1D::computeI0<15>()
{
    return -omegaPower_[15]*(g1Power_[1]*gPower_[14] +
        273*g1Power_[11]*gPower_[4] + 35*g1Power_[13]*gPower_[2] + g1Power_[15]
        + 35*g1Power_[3]*gPower_[12] + 273*g1Power_[5]*gPower_[10] +
        715*g1Power_[7]*gPower_[8] + 715*g1Power_[9]*gPower_[6])/32768 +
        omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/16384*v1Power_[1] -
        15*omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/8192*v1Power_[2] +
        35*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/4096*v1Power_[3] -
        455*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/2048*v1Power_[4]
        + 273*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/1024*v1Power_[5] -
        1001*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/512*v1Power_[6] +
        715*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/256*v1Power_[7] -
        6435*omegaPower_[7]*(g1Power_[1]*gPower_[6] + 7*g1Power_[3]*gPower_[4] +
        7*g1Power_[5]*gPower_[2] + g1Power_[7])/128*v1Power_[8] +
        715*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/64*v1Power_[9] -
        1001*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/32*v1Power_[10] +
        273*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/16*v1Power_[11] - 455*omegaPower_[3]*(g1Power_[1]*gPower_[2]
        + g1Power_[3])/8*v1Power_[12] + 35*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[13] -
        15*g1Power_[1]*omegaPower_[1]/2*v1Power_[14];
}


template<>
double HardSphereCollision1D::computeI0<16>()
{
    return omegaPower_[16]*(19448*g1Power_[10]*gPower_[6] +
        6188*g1Power_[12]*gPower_[4] + 680*g1Power_[14]*gPower_[2] +
        17*g1Power_[16] + 136*g1Power_[2]*gPower_[14] +
        2380*g1Power_[4]*gPower_[12] + 12376*g1Power_[6]*gPower_[10] +
        24310*g1Power_[8]*gPower_[8] + gPower_[16])/1114112 -
        omegaPower_[15]*(g1Power_[1]*gPower_[14] + 273*g1Power_[11]*gPower_[4] +
        35*g1Power_[13]*gPower_[2] + g1Power_[15] + 35*g1Power_[3]*gPower_[12] +
        273*g1Power_[5]*gPower_[10] + 715*g1Power_[7]*gPower_[8] +
        715*g1Power_[9]*gPower_[6])/2048*v1Power_[1] +
        omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/2048*v1Power_[2] -
        5*omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/512*v1Power_[3] +
        35*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/1024*v1Power_[4] -
        91*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/128*v1Power_[5]
        + 91*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/128*v1Power_[6] -
        143*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] + 60*g1Power_[3]*gPower_[6]
        + 126*g1Power_[5]*gPower_[4] + 60*g1Power_[7]*gPower_[2] +
        5*g1Power_[9])/32*v1Power_[7] +
        715*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/128*v1Power_[8] - 715*omegaPower_[7]*(g1Power_[1]*gPower_[6]
        + 7*g1Power_[3]*gPower_[4] + 7*g1Power_[5]*gPower_[2] +
        g1Power_[7])/8*v1Power_[9] +
        143*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/8*v1Power_[10] -
        91*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] + 10*g1Power_[3]*gPower_[2]
        + 3*g1Power_[5])/2*v1Power_[11] +
        91*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/4*v1Power_[12] - 70*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])*v1Power_[13] + 10*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])*v1Power_[14] - 8*g1Power_[1]*omegaPower_[1]*v1Power_[15];
}


template<>
double HardSphereCollision1D::computeI0<17>()
{
    return -omegaPower_[17]*(9*g1Power_[1]*gPower_[16] +
        15912*g1Power_[11]*gPower_[6] + 4284*g1Power_[13]*gPower_[4] +
        408*g1Power_[15]*gPower_[2] + 9*g1Power_[17] +
        408*g1Power_[3]*gPower_[14] + 4284*g1Power_[5]*gPower_[12] +
        15912*g1Power_[7]*gPower_[10] + 24310*g1Power_[9]*gPower_[8])/1179648 +
        omegaPower_[16]*(19448*g1Power_[10]*gPower_[6] +
        6188*g1Power_[12]*gPower_[4] + 680*g1Power_[14]*gPower_[2] +
        17*g1Power_[16] + 136*g1Power_[2]*gPower_[14] +
        2380*g1Power_[4]*gPower_[12] + 12376*g1Power_[6]*gPower_[10] +
        24310*g1Power_[8]*gPower_[8] + gPower_[16])/65536*v1Power_[1] -
        17*omegaPower_[15]*(g1Power_[1]*gPower_[14] +
        273*g1Power_[11]*gPower_[4] + 35*g1Power_[13]*gPower_[2] + g1Power_[15]
        + 35*g1Power_[3]*gPower_[12] + 273*g1Power_[5]*gPower_[10] +
        715*g1Power_[7]*gPower_[8] +
        715*g1Power_[9]*gPower_[6])/4096*v1Power_[2] +
        17*omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/6144*v1Power_[3] -
        85*omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/2048*v1Power_[4] +
        119*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/1024*v1Power_[5] -
        1547*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/768*v1Power_[6]
        + 221*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/128*v1Power_[7] -
        2431*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/256*v1Power_[8] +
        12155*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/1152*v1Power_[9] -
        2431*omegaPower_[7]*(g1Power_[1]*gPower_[6] + 7*g1Power_[3]*gPower_[4] +
        7*g1Power_[5]*gPower_[2] + g1Power_[7])/16*v1Power_[10] +
        221*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/8*v1Power_[11] -
        1547*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/24*v1Power_[12] +
        119*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/4*v1Power_[13] - 85*omegaPower_[3]*(g1Power_[1]*gPower_[2] +
        g1Power_[3])*v1Power_[14] + 34*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/3*v1Power_[15] -
        17*g1Power_[1]*omegaPower_[1]/2*v1Power_[16];
}


template<>
double HardSphereCollision1D::computeI0<18>()
{
    return omegaPower_[18]*(92378*g1Power_[10]*gPower_[8] +
        50388*g1Power_[12]*gPower_[6] + 11628*g1Power_[14]*gPower_[4] +
        969*g1Power_[16]*gPower_[2] + 19*g1Power_[18] +
        171*g1Power_[2]*gPower_[16] + 3876*g1Power_[4]*gPower_[14] +
        27132*g1Power_[6]*gPower_[12] + 75582*g1Power_[8]*gPower_[10] +
        gPower_[18])/4980736 - omegaPower_[17]*(9*g1Power_[1]*gPower_[16] +
        15912*g1Power_[11]*gPower_[6] + 4284*g1Power_[13]*gPower_[4] +
        408*g1Power_[15]*gPower_[2] + 9*g1Power_[17] +
        408*g1Power_[3]*gPower_[14] + 4284*g1Power_[5]*gPower_[12] +
        15912*g1Power_[7]*gPower_[10] +
        24310*g1Power_[9]*gPower_[8])/65536*v1Power_[1] +
        9*omegaPower_[16]*(19448*g1Power_[10]*gPower_[6] +
        6188*g1Power_[12]*gPower_[4] + 680*g1Power_[14]*gPower_[2] +
        17*g1Power_[16] + 136*g1Power_[2]*gPower_[14] +
        2380*g1Power_[4]*gPower_[12] + 12376*g1Power_[6]*gPower_[10] +
        24310*g1Power_[8]*gPower_[8] + gPower_[16])/65536*v1Power_[2] -
        51*omegaPower_[15]*(g1Power_[1]*gPower_[14] +
        273*g1Power_[11]*gPower_[4] + 35*g1Power_[13]*gPower_[2] + g1Power_[15]
        + 35*g1Power_[3]*gPower_[12] + 273*g1Power_[5]*gPower_[10] +
        715*g1Power_[7]*gPower_[8] +
        715*g1Power_[9]*gPower_[6])/2048*v1Power_[3] +
        51*omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/4096*v1Power_[4] -
        153*omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/1024*v1Power_[5] +
        357*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/1024*v1Power_[6] -
        663*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/128*v1Power_[7]
        + 1989*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/512*v1Power_[8] -
        2431*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/128*v1Power_[9] +
        2431*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/128*v1Power_[10] -
        1989*omegaPower_[7]*(g1Power_[1]*gPower_[6] + 7*g1Power_[3]*gPower_[4] +
        7*g1Power_[5]*gPower_[2] + g1Power_[7])/8*v1Power_[11] +
        663*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/16*v1Power_[12]
        - 357*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/4*v1Power_[13] +
        153*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/4*v1Power_[14] - 102*omegaPower_[3]*(g1Power_[1]*gPower_[2]
        + g1Power_[3])*v1Power_[15] + 51*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[16] - 9*g1Power_[1]*omegaPower_[1]*v1Power_[17];
}


template<>
double HardSphereCollision1D::computeI0<19>()
{
    return -omegaPower_[19]*(5*g1Power_[1]*gPower_[18] +
        41990*g1Power_[11]*gPower_[8] + 19380*g1Power_[13]*gPower_[6] +
        3876*g1Power_[15]*gPower_[4] + 285*g1Power_[17]*gPower_[2] +
        5*g1Power_[19] + 285*g1Power_[3]*gPower_[16] +
        3876*g1Power_[5]*gPower_[14] + 19380*g1Power_[7]*gPower_[12] +
        41990*g1Power_[9]*gPower_[10])/2621440 +
        omegaPower_[18]*(92378*g1Power_[10]*gPower_[8] +
        50388*g1Power_[12]*gPower_[6] + 11628*g1Power_[14]*gPower_[4] +
        969*g1Power_[16]*gPower_[2] + 19*g1Power_[18] +
        171*g1Power_[2]*gPower_[16] + 3876*g1Power_[4]*gPower_[14] +
        27132*g1Power_[6]*gPower_[12] + 75582*g1Power_[8]*gPower_[10] +
        gPower_[18])/262144*v1Power_[1] -
        19*omegaPower_[17]*(9*g1Power_[1]*gPower_[16] +
        15912*g1Power_[11]*gPower_[6] + 4284*g1Power_[13]*gPower_[4] +
        408*g1Power_[15]*gPower_[2] + 9*g1Power_[17] +
        408*g1Power_[3]*gPower_[14] + 4284*g1Power_[5]*gPower_[12] +
        15912*g1Power_[7]*gPower_[10] +
        24310*g1Power_[9]*gPower_[8])/131072*v1Power_[2] +
        57*omegaPower_[16]*(19448*g1Power_[10]*gPower_[6] +
        6188*g1Power_[12]*gPower_[4] + 680*g1Power_[14]*gPower_[2] +
        17*g1Power_[16] + 136*g1Power_[2]*gPower_[14] +
        2380*g1Power_[4]*gPower_[12] + 12376*g1Power_[6]*gPower_[10] +
        24310*g1Power_[8]*gPower_[8] + gPower_[16])/65536*v1Power_[3] -
        969*omegaPower_[15]*(g1Power_[1]*gPower_[14] +
        273*g1Power_[11]*gPower_[4] + 35*g1Power_[13]*gPower_[2] + g1Power_[15]
        + 35*g1Power_[3]*gPower_[12] + 273*g1Power_[5]*gPower_[10] +
        715*g1Power_[7]*gPower_[8] +
        715*g1Power_[9]*gPower_[6])/8192*v1Power_[4] +
        969*omegaPower_[14]*(3003*g1Power_[10]*gPower_[4] +
        455*g1Power_[12]*gPower_[2] + 15*g1Power_[14] +
        105*g1Power_[2]*gPower_[12] + 1365*g1Power_[4]*gPower_[10] +
        5005*g1Power_[6]*gPower_[8] + 6435*g1Power_[8]*gPower_[6] +
        gPower_[14])/20480*v1Power_[5] -
        969*omegaPower_[13]*(7*g1Power_[1]*gPower_[12] +
        182*g1Power_[11]*gPower_[2] + 7*g1Power_[13] +
        182*g1Power_[3]*gPower_[10] + 1001*g1Power_[5]*gPower_[8] +
        1716*g1Power_[7]*gPower_[6] +
        1001*g1Power_[9]*gPower_[4])/2048*v1Power_[6] +
        969*omegaPower_[12]*(286*g1Power_[10]*gPower_[2] + 13*g1Power_[12] +
        78*g1Power_[2]*gPower_[10] + 715*g1Power_[4]*gPower_[8] +
        1716*g1Power_[6]*gPower_[6] + 1287*g1Power_[8]*gPower_[4] +
        gPower_[12])/1024*v1Power_[7] -
        12597*omegaPower_[11]*(3*g1Power_[1]*gPower_[10] + 3*g1Power_[11] +
        55*g1Power_[3]*gPower_[8] + 198*g1Power_[5]*gPower_[6] +
        198*g1Power_[7]*gPower_[4] + 55*g1Power_[9]*gPower_[2])/1024*v1Power_[8]
        + 4199*omegaPower_[10]*(11*g1Power_[10] + 55*g1Power_[2]*gPower_[8] +
        330*g1Power_[4]*gPower_[6] + 462*g1Power_[6]*gPower_[4] +
        165*g1Power_[8]*gPower_[2] + gPower_[10])/512*v1Power_[9] -
        46189*omegaPower_[9]*(5*g1Power_[1]*gPower_[8] +
        60*g1Power_[3]*gPower_[6] + 126*g1Power_[5]*gPower_[4] +
        60*g1Power_[7]*gPower_[2] + 5*g1Power_[9])/1280*v1Power_[10] +
        4199*omegaPower_[8]*(36*g1Power_[2]*gPower_[6] +
        126*g1Power_[4]*gPower_[4] + 84*g1Power_[6]*gPower_[2] + 9*g1Power_[8] +
        gPower_[8])/128*v1Power_[11] -
        12597*omegaPower_[7]*(g1Power_[1]*gPower_[6] + 7*g1Power_[3]*gPower_[4]
        + 7*g1Power_[5]*gPower_[2] + g1Power_[7])/32*v1Power_[12] +
        969*omegaPower_[6]*(21*g1Power_[2]*gPower_[4] +
        35*g1Power_[4]*gPower_[2] + 7*g1Power_[6] + gPower_[6])/16*v1Power_[13]
        - 969*omegaPower_[5]*(3*g1Power_[1]*gPower_[4] +
        10*g1Power_[3]*gPower_[2] + 3*g1Power_[5])/8*v1Power_[14] +
        969*omegaPower_[4]*(10*g1Power_[2]*gPower_[2] + 5*g1Power_[4] +
        gPower_[4])/20*v1Power_[15] - 969*omegaPower_[3]*(g1Power_[1]*gPower_[2]
        + g1Power_[3])/8*v1Power_[16] + 57*omegaPower_[2]*(3*g1Power_[2] +
        gPower_[2])/4*v1Power_[17] -
        19*g1Power_[1]*omegaPower_[1]/2*v1Power_[18];
}


inline std::vector<std::function<double()> >
    HardSphereCollision1D::initializeI0Functions(int nMoments)
{

    if (nMoments <= 0) {
        std::string errorMessage = 
            "The parameter `nMoments` must be a positive integer.";
        throw std::runtime_error(errorMessage);
    }

    std::vector<std::function<double()> > functions;

    static_assert(nMomentsMax_ > 0);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<0>, this)
    );
    if (nMoments == 1)
        return functions;

    static_assert(nMomentsMax_ > 1);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<1>, this)
    );
    if (nMoments == 2)
        return functions;

    static_assert(nMomentsMax_ > 2);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<2>, this)
    );
    if (nMoments == 3)
        return functions;

    static_assert(nMomentsMax_ > 3);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<3>, this)
    );
    if (nMoments == 4)
        return functions;

    static_assert(nMomentsMax_ > 4);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<4>, this)
    );
    if (nMoments == 5)
        return functions;

    static_assert(nMomentsMax_ > 5);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<5>, this)
    );
    if (nMoments == 6)
        return functions;

    static_assert(nMomentsMax_ > 6);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<6>, this)
    );
    if (nMoments == 7)
        return functions;

    static_assert(nMomentsMax_ > 7);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<7>, this)
    );
    if (nMoments == 8)
        return functions;

    static_assert(nMomentsMax_ > 8);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<8>, this)
    );
    if (nMoments == 9)
        return functions;

    static_assert(nMomentsMax_ > 9);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<9>, this)
    );
    if (nMoments == 10)
        return functions;

    static_assert(nMomentsMax_ > 10);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<10>, this)
    );
    if (nMoments == 11)
        return functions;

    static_assert(nMomentsMax_ > 11);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<11>, this)
    );
    if (nMoments == 12)
        return functions;

    static_assert(nMomentsMax_ > 12);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<12>, this)
    );
    if (nMoments == 13)
        return functions;

    static_assert(nMomentsMax_ > 13);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<13>, this)
    );
    if (nMoments == 14)
        return functions;

    static_assert(nMomentsMax_ > 14);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<14>, this)
    );
    if (nMoments == 15)
        return functions;

    static_assert(nMomentsMax_ > 15);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<15>, this)
    );
    if (nMoments == 16)
        return functions;

    static_assert(nMomentsMax_ > 16);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<16>, this)
    );
    if (nMoments == 17)
        return functions;

    static_assert(nMomentsMax_ > 17);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<17>, this)
    );
    if (nMoments == 18)
        return functions;

    static_assert(nMomentsMax_ > 18);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<18>, this)
    );
    if (nMoments == 19)
        return functions;

    static_assert(nMomentsMax_ > 19);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI0<19>, this)
    );
    if (nMoments == 20)
        return functions;
    

    static_assert(nMomentsMax_ == 20);
    std::string errorMessage = 
        "The maximum allowed number of moments `nMoments` is "
        + std::to_string(nMomentsMax_);
    throw std::runtime_error(errorMessage);
}


#endif // COLLISION_INTEGRAL_FUNCTIONS_HPP