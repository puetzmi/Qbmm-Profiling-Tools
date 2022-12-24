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

template<>
double HardSphereCollision1D::computeI1<0>()
{
    return 0;
}


template<>
double HardSphereCollision1D::computeI1<1>()
{
    return -2*omegaPower_[1]*(2*g1Power_[2] + gPower_[2])/15;
}


template<>
double HardSphereCollision1D::computeI1<2>()
{
    return -4*omegaPower_[1]*v1Power_[1]*(2*g1Power_[2] + gPower_[2])/15
        + 2*omegaPower_[2]*(3*g1Power_[1]*gPower_[2] + 2*g1Power_[3])/35;
}


template<>
double HardSphereCollision1D::computeI1<3>()
{
    return -2*omegaPower_[1]*v1Power_[2]*(2*g1Power_[2] + gPower_[2])/5
        + 6*omegaPower_[2]*v1Power_[1]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 2*omegaPower_[3]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/315;
}


template<>
double HardSphereCollision1D::computeI1<4>()
{
    return -8*omegaPower_[1]*v1Power_[3]*(2*g1Power_[2] + gPower_[2])/15
        + 12*omegaPower_[2]*v1Power_[2]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 8*omegaPower_[3]*v1Power_[1]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/315 +
        2*omegaPower_[4]*(8*g1Power_[5] + 5*gPower_[1]*(3*g1Power_[1]*gPower_[3]
        + 8*g1Power_[3]*gPower_[1]))/693;
}


template<>
double HardSphereCollision1D::computeI1<5>()
{
    return -2*omegaPower_[1]*v1Power_[4]*(2*g1Power_[2] + gPower_[2])/3
        + 4*omegaPower_[2]*v1Power_[3]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/7 - 4*omegaPower_[3]*v1Power_[2]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/63 +
        10*omegaPower_[4]*v1Power_[1]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/693
        - 2*omegaPower_[5]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/3003;
}


template<>
double HardSphereCollision1D::computeI1<6>()
{
    return -4*omegaPower_[1]*v1Power_[5]*(2*g1Power_[2] + gPower_[2])/5
        + 6*omegaPower_[2]*v1Power_[4]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/7 - 8*omegaPower_[3]*v1Power_[3]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/63 +
        10*omegaPower_[4]*v1Power_[2]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/231
        - 4*omegaPower_[5]*v1Power_[1]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/1001 + 2*omegaPower_[6]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/6435;
}


template<>
double HardSphereCollision1D::computeI1<7>()
{
    return -14*omegaPower_[1]*v1Power_[6]*(2*g1Power_[2] +
        gPower_[2])/15 + 6*omegaPower_[2]*v1Power_[5]*(3*g1Power_[1]*gPower_[2]
        + 2*g1Power_[3])/5 - 2*omegaPower_[3]*v1Power_[4]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/9 +
        10*omegaPower_[4]*v1Power_[3]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/99 -
        2*omegaPower_[5]*v1Power_[2]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/143 + 14*omegaPower_[6]*v1Power_[1]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/6435 - 2*omegaPower_[7]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/109395;
}


template<>
double HardSphereCollision1D::computeI1<8>()
{
    return -16*omegaPower_[1]*v1Power_[7]*(2*g1Power_[2] +
        gPower_[2])/15 + 8*omegaPower_[2]*v1Power_[6]*(3*g1Power_[1]*gPower_[2]
        + 2*g1Power_[3])/5 - 16*omegaPower_[3]*v1Power_[5]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/45 +
        20*omegaPower_[4]*v1Power_[4]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/99 -
        16*omegaPower_[5]*v1Power_[3]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/429 + 56*omegaPower_[6]*v1Power_[2]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/6435 -
        16*omegaPower_[7]*v1Power_[1]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/109395 +
        2*omegaPower_[8]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/230945;
}


template<>
double HardSphereCollision1D::computeI1<9>()
{
    return -6*omegaPower_[1]*v1Power_[8]*(2*g1Power_[2] + gPower_[2])/5
        + 72*omegaPower_[2]*v1Power_[7]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 8*omegaPower_[3]*v1Power_[6]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/15 +
        4*omegaPower_[4]*v1Power_[5]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/11 -
        12*omegaPower_[5]*v1Power_[4]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/143 + 56*omegaPower_[6]*v1Power_[3]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/2145 -
        8*omegaPower_[7]*v1Power_[2]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/12155 +
        18*omegaPower_[8]*v1Power_[1]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/230945 -
        2*omegaPower_[9]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/969969;
}


template<>
double HardSphereCollision1D::computeI1<10>()
{
    return -4*omegaPower_[1]*v1Power_[9]*(2*g1Power_[2] + gPower_[2])/3
        + 2*omegaPower_[10]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/2028117 +
        18*omegaPower_[2]*v1Power_[8]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/7 - 16*omegaPower_[3]*v1Power_[7]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/21 +
        20*omegaPower_[4]*v1Power_[6]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/33 -
        24*omegaPower_[5]*v1Power_[5]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/143 + 28*omegaPower_[6]*v1Power_[4]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/429 -
        16*omegaPower_[7]*v1Power_[3]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/7293 +
        18*omegaPower_[8]*v1Power_[2]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/46189 -
        20*omegaPower_[9]*v1Power_[1]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/969969;
}


template<>
double HardSphereCollision1D::computeI1<11>()
{
    return -22*omegaPower_[1]*v1Power_[10]*(2*g1Power_[2] +
        gPower_[2])/15 + 22*omegaPower_[10]*v1Power_[1]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/2028117 -
        2*omegaPower_[11]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/16900975 +
        22*omegaPower_[2]*v1Power_[9]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/7 - 22*omegaPower_[3]*v1Power_[8]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/21 +
        20*omegaPower_[4]*v1Power_[7]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/21 -
        4*omegaPower_[5]*v1Power_[6]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/13 + 28*omegaPower_[6]*v1Power_[5]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/195 -
        4*omegaPower_[7]*v1Power_[4]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/663 +
        6*omegaPower_[8]*v1Power_[3]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/4199 -
        10*omegaPower_[9]*v1Power_[2]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/88179;
}


template<>
double HardSphereCollision1D::computeI1<12>()
{
    return -8*omegaPower_[1]*v1Power_[11]*(2*g1Power_[2] + gPower_[2])/5
        + 44*omegaPower_[10]*v1Power_[2]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/676039 -
        24*omegaPower_[11]*v1Power_[1]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/16900975 +
        2*omegaPower_[12]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/35102025 +
        132*omegaPower_[2]*v1Power_[10]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 88*omegaPower_[3]*v1Power_[9]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/63 +
        10*omegaPower_[4]*v1Power_[8]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/7 -
        48*omegaPower_[5]*v1Power_[7]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/91 + 56*omegaPower_[6]*v1Power_[6]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/195 -
        16*omegaPower_[7]*v1Power_[5]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/1105 +
        18*omegaPower_[8]*v1Power_[4]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/4199 -
        40*omegaPower_[9]*v1Power_[3]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/88179;
}


template<>
double HardSphereCollision1D::computeI1<13>()
{
    return -26*omegaPower_[1]*v1Power_[12]*(2*g1Power_[2] +
        gPower_[2])/15 + 44*omegaPower_[10]*v1Power_[3]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/156009 -
        12*omegaPower_[11]*v1Power_[2]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/1300075 +
        26*omegaPower_[12]*v1Power_[1]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/35102025 -
        2*omegaPower_[13]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/145422675 +
        156*omegaPower_[2]*v1Power_[11]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 572*omegaPower_[3]*v1Power_[10]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/315 +
        130*omegaPower_[4]*v1Power_[9]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/63 -
        6*omegaPower_[5]*v1Power_[8]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/7 + 8*omegaPower_[6]*v1Power_[7]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/15 -
        8*omegaPower_[7]*v1Power_[6]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/255 +
        18*omegaPower_[8]*v1Power_[5]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/1615 -
        10*omegaPower_[9]*v1Power_[4]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/6783;
}


template<>
double HardSphereCollision1D::computeI1<14>()
{
    return -28*omegaPower_[1]*v1Power_[13]*(2*g1Power_[2] +
        gPower_[2])/15 + 22*omegaPower_[10]*v1Power_[4]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/22287 -
        8*omegaPower_[11]*v1Power_[3]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/185725 +
        26*omegaPower_[12]*v1Power_[2]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/5014575 -
        28*omegaPower_[13]*v1Power_[1]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/145422675 +
        6.65468391008397e-9*omegaPower_[14]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) +
        26*omegaPower_[2]*v1Power_[12]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/5 - 104*omegaPower_[3]*v1Power_[11]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/45 +
        26*omegaPower_[4]*v1Power_[10]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/9 -
        4*omegaPower_[5]*v1Power_[9]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/3 + 14*omegaPower_[6]*v1Power_[8]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/15 -
        16*omegaPower_[7]*v1Power_[7]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/255 +
        42*omegaPower_[8]*v1Power_[6]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/1615 -
        4*omegaPower_[9]*v1Power_[5]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/969;
}


template<>
double HardSphereCollision1D::computeI1<15>()
{
    return 18*g1Power_[1]*gPower_[2]*omegaPower_[2] +
        12*g1Power_[3]*omegaPower_[2]*v1Power_[13] -
        2*omegaPower_[1]*v1Power_[14]*(2*g1Power_[2] + gPower_[2]) +
        22*omegaPower_[10]*v1Power_[5]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/7429 -
        6*omegaPower_[11]*v1Power_[4]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/37145 +
        26*omegaPower_[12]*v1Power_[3]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/1002915 -
        14*omegaPower_[13]*v1Power_[2]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/9694845 +
        9.98202586512596e-8*omegaPower_[14]*v1Power_[1]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) -
        3.15823057106504e-19*omegaPower_[15]*(20922789887999.0*g1Power_[16] +
        638512875*gPower_[1]*(82001920*g1Power_[10]*gPower_[5] +
        22364160*g1Power_[12]*gPower_[3] + 1966080*g1Power_[14]*gPower_[1] +
        823680*g1Power_[2]*gPower_[13] + 13453440*g1Power_[4]*gPower_[11] +
        64576512*g1Power_[6]*gPower_[9] + 115315200*g1Power_[8]*gPower_[7] +
        6435*gPower_[15])) - 26*omegaPower_[3]*v1Power_[12]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/9 +
        130*omegaPower_[4]*v1Power_[11]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/33 -
        2*omegaPower_[5]*v1Power_[10]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5])) + 14*omegaPower_[6]*v1Power_[9]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/9 -
        2*omegaPower_[7]*v1Power_[8]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/17 +
        18*omegaPower_[8]*v1Power_[7]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/323 -
        10*omegaPower_[9]*v1Power_[6]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/969;
}


template<>
double HardSphereCollision1D::computeI1<16>()
{
    return -32*omegaPower_[1]*v1Power_[15]*(2*g1Power_[2] +
        gPower_[2])/15 + 176*omegaPower_[10]*v1Power_[6]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/22287 -
        96*omegaPower_[11]*v1Power_[5]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/185725 +
        104*omegaPower_[12]*v1Power_[4]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/1002915 -
        224*omegaPower_[13]*v1Power_[3]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/29084535 +
        7.98562069210077e-7*omegaPower_[14]*v1Power_[2]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) -
        5.05316891370406e-18*omegaPower_[15]*v1Power_[1]*(20922789887999.0*g1Power_[16]
        + 638512875*gPower_[1]*(82001920*g1Power_[10]*gPower_[5] +
        22364160*g1Power_[12]*gPower_[3] + 1966080*g1Power_[14]*gPower_[1] +
        823680*g1Power_[2]*gPower_[13] + 13453440*g1Power_[4]*gPower_[11] +
        64576512*g1Power_[6]*gPower_[9] + 115315200*g1Power_[8]*gPower_[7] +
        6435*gPower_[15])) +
        9.02351591732868e-21*omegaPower_[16]*(355687428095815.0*g1Power_[17] +
        17*gPower_[1]*(69850115960625.0*g1Power_[1]*gPower_[15] +
        80918889891840000.0*g1Power_[11]*gPower_[5] +
        18673589975040000.0*g1Power_[13]*gPower_[3] +
        1422749712383992.0*g1Power_[15]*gPower_[1] +
        2980271614320000.0*g1Power_[3]*gPower_[13] +
        29206661820336000.0*g1Power_[5]*gPower_[11] +
        100137126241152000.0*g1Power_[7]*gPower_[9] +
        139079342001600000.0*g1Power_[9]*gPower_[7])) +
        48*omegaPower_[2]*v1Power_[14]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/7 - 32*omegaPower_[3]*v1Power_[13]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/9 +
        520*omegaPower_[4]*v1Power_[12]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/99 -
        32*omegaPower_[5]*v1Power_[11]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/11 + 112*omegaPower_[6]*v1Power_[10]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/45 -
        32*omegaPower_[7]*v1Power_[9]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/153 +
        36*omegaPower_[8]*v1Power_[8]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/323 -
        160*omegaPower_[9]*v1Power_[7]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/6783;
}


template<>
double HardSphereCollision1D::computeI1<17>()
{
    return -34*omegaPower_[1]*v1Power_[16]*(2*g1Power_[2] +
        gPower_[2])/15 + 176*omegaPower_[10]*v1Power_[7]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/9177 -
        16*omegaPower_[11]*v1Power_[6]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/10925 +
        104*omegaPower_[12]*v1Power_[5]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/294975 -
        56*omegaPower_[13]*v1Power_[4]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/1710855 +
        4.5251850588571e-6*omegaPower_[14]*v1Power_[3]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) -
        4.29519357664845e-17*omegaPower_[15]*v1Power_[2]*(20922789887999.0*g1Power_[16]
        + 638512875*gPower_[1]*(82001920*g1Power_[10]*gPower_[5] +
        22364160*g1Power_[12]*gPower_[3] + 1966080*g1Power_[14]*gPower_[1] +
        823680*g1Power_[2]*gPower_[13] + 13453440*g1Power_[4]*gPower_[11] +
        64576512*g1Power_[6]*gPower_[9] + 115315200*g1Power_[8]*gPower_[7] +
        6435*gPower_[15])) +
        1.53399770594587e-19*omegaPower_[16]*v1Power_[1]*(355687428095815.0*g1Power_[17]
        + 17*gPower_[1]*(69850115960625.0*g1Power_[1]*gPower_[15] +
        80918889891840000.0*g1Power_[11]*gPower_[5] +
        18673589975040000.0*g1Power_[13]*gPower_[3] +
        1422749712383992.0*g1Power_[15]*gPower_[1] +
        2980271614320000.0*g1Power_[3]*gPower_[13] +
        29206661820336000.0*g1Power_[5]*gPower_[11] +
        100137126241152000.0*g1Power_[7]*gPower_[9] +
        139079342001600000.0*g1Power_[9]*gPower_[7])) -
        2.43878808576451e-22*omegaPower_[17]*(6402373705847470.0*g1Power_[18] +
        51*gPower_[1]*(1502056893617280000.0*g1Power_[10]*gPower_[7] +
        728270009026558180.0*g1Power_[12]*gPower_[5] +
        144053408378885280.0*g1Power_[14]*gPower_[3] +
        9603560558585937.0*g1Power_[16]*gPower_[1] +
        3771906261873750.0*g1Power_[2]*gPower_[15] +
        80467333586640000.0*g1Power_[4]*gPower_[13] +
        525719912766048000.0*g1Power_[6]*gPower_[11] +
        1351851204255552000.0*g1Power_[8]*gPower_[9] +
        23283371986875.0*gPower_[17])) +
        272*omegaPower_[2]*v1Power_[15]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 272*omegaPower_[3]*v1Power_[14]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/63 +
        680*omegaPower_[4]*v1Power_[13]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/99 -
        136*omegaPower_[5]*v1Power_[12]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/33 + 1904*omegaPower_[6]*v1Power_[11]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/495 -
        16*omegaPower_[7]*v1Power_[10]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/45 +
        4*omegaPower_[8]*v1Power_[9]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/19 -
        20*omegaPower_[9]*v1Power_[8]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/399;
}


template<>
double HardSphereCollision1D::computeI1<18>()
{
    return -12*omegaPower_[1]*v1Power_[17]*(2*g1Power_[2] +
        gPower_[2])/5 + 132*omegaPower_[10]*v1Power_[8]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/3059 -
        288*omegaPower_[11]*v1Power_[7]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/76475 +
        104*omegaPower_[12]*v1Power_[6]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/98325 -
        112*omegaPower_[13]*v1Power_[5]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/950475 +
        2.0363332764857e-5*omegaPower_[14]*v1Power_[4]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) -
        2.57711614598907e-16*omegaPower_[15]*v1Power_[3]*(20922789887999.0*g1Power_[16]
        + 638512875*gPower_[1]*(82001920*g1Power_[10]*gPower_[5] +
        22364160*g1Power_[12]*gPower_[3] + 1966080*g1Power_[14]*gPower_[1] +
        823680*g1Power_[2]*gPower_[13] + 13453440*g1Power_[4]*gPower_[11] +
        64576512*g1Power_[6]*gPower_[9] + 115315200*g1Power_[8]*gPower_[7] +
        6435*gPower_[15])) +
        1.38059793535129e-18*omegaPower_[16]*v1Power_[2]*(355687428095815.0*g1Power_[17]
        + 17*gPower_[1]*(69850115960625.0*g1Power_[1]*gPower_[15] +
        80918889891840000.0*g1Power_[11]*gPower_[5] +
        18673589975040000.0*g1Power_[13]*gPower_[3] +
        1422749712383992.0*g1Power_[15]*gPower_[1] +
        2980271614320000.0*g1Power_[3]*gPower_[13] +
        29206661820336000.0*g1Power_[5]*gPower_[11] +
        100137126241152000.0*g1Power_[7]*gPower_[9] +
        139079342001600000.0*g1Power_[9]*gPower_[7])) -
        4.38981855437611e-21*omegaPower_[17]*v1Power_[1]*(6402373705847470.0*g1Power_[18]
        + 51*gPower_[1]*(1502056893617280000.0*g1Power_[10]*gPower_[7] +
        728270009026558180.0*g1Power_[12]*gPower_[5] +
        144053408378885280.0*g1Power_[14]*gPower_[3] +
        9603560558585937.0*g1Power_[16]*gPower_[1] +
        3771906261873750.0*g1Power_[2]*gPower_[15] +
        80467333586640000.0*g1Power_[4]*gPower_[13] +
        525719912766048000.0*g1Power_[6]*gPower_[11] +
        1351851204255552000.0*g1Power_[8]*gPower_[9] +
        23283371986875.0*gPower_[17])) +
        1.87599083520347e-23*omegaPower_[18]*(40548366800373008.0*g1Power_[19] +
        19*gPower_[1]*(7520529151760625.0*g1Power_[1]*gPower_[17] +
        44105852421670993590.0*g1Power_[11]*gPower_[7] +
        18094708685814014660.0*g1Power_[13]*gPower_[5] +
        3101950060424664996.0*g1Power_[15]*gPower_[3] +
        182467650613731819.0*g1Power_[17]*gPower_[1] +
        406108574195073750.0*g1Power_[3]*gPower_[15] +
        5198189749696944000.0*g1Power_[5]*gPower_[13] +
        24258218831919072000.0*g1Power_[7]*gPower_[11] +
        48516437663838144000.0*g1Power_[9]*gPower_[9])) +
        306*omegaPower_[2]*v1Power_[16]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 544*omegaPower_[3]*v1Power_[15]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/105 +
        680*omegaPower_[4]*v1Power_[14]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/77 -
        816*omegaPower_[5]*v1Power_[13]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/143 + 952*omegaPower_[6]*v1Power_[12]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/165 -
        32*omegaPower_[7]*v1Power_[11]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/55 +
        36*omegaPower_[8]*v1Power_[10]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/95 -
        40*omegaPower_[9]*v1Power_[9]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/399;
}


template<>
double HardSphereCollision1D::computeI1<19>()
{
    return -38*omegaPower_[1]*v1Power_[18]*(2*g1Power_[2] +
        gPower_[2])/15 + 44*omegaPower_[10]*v1Power_[9]*(256*g1Power_[11] +
        11*gPower_[1]*(63*g1Power_[1]*gPower_[9] + 1050*g1Power_[3]*gPower_[7] +
        3360*g1Power_[5]*gPower_[5] + 2880*g1Power_[7]*gPower_[3] +
        640*g1Power_[9]*gPower_[1]))/483 -
        36*omegaPower_[11]*v1Power_[8]*(1024*g1Power_[12] +
        33*gPower_[1]*(1024*g1Power_[10]*gPower_[1] + 504*g1Power_[2]*gPower_[9]
        + 4200*g1Power_[4]*gPower_[7] + 8960*g1Power_[6]*gPower_[5] +
        5760*g1Power_[8]*gPower_[3] + 7*gPower_[11]))/4025 +
        104*omegaPower_[12]*v1Power_[7]*(1024*g1Power_[13] +
        39*gPower_[1]*(77*g1Power_[1]*gPower_[11] + 1024*g1Power_[11]*gPower_[1]
        + 1848*g1Power_[3]*gPower_[9] + 9240*g1Power_[5]*gPower_[7] +
        14080*g1Power_[7]*gPower_[5] + 7040*g1Power_[9]*gPower_[3]))/36225 -
        56*omegaPower_[13]*v1Power_[6]*(2048*g1Power_[14] +
        13*gPower_[1]*(59136*g1Power_[10]*gPower_[3] +
        7168*g1Power_[12]*gPower_[1] + 3234*g1Power_[2]*gPower_[11] +
        38808*g1Power_[4]*gPower_[9] + 129360*g1Power_[6]*gPower_[7] +
        147840*g1Power_[8]*gPower_[5] + 33*gPower_[13]))/150075 +
        7.73806645064564e-5*omegaPower_[14]*v1Power_[5]*(2048*g1Power_[15] +
        gPower_[1]*(6435*g1Power_[1]*gPower_[13] +
        1048320*g1Power_[11]*gPower_[3] + 107520*g1Power_[13]*gPower_[1] +
        210210*g1Power_[3]*gPower_[11] + 1513512*g1Power_[5]*gPower_[9] +
        3603600*g1Power_[7]*gPower_[7] + 3203200*g1Power_[9]*gPower_[5])) -
        1.22413016934481e-15*omegaPower_[15]*v1Power_[4]*(20922789887999.0*g1Power_[16]
        + 638512875*gPower_[1]*(82001920*g1Power_[10]*gPower_[5] +
        22364160*g1Power_[12]*gPower_[3] + 1966080*g1Power_[14]*gPower_[1] +
        823680*g1Power_[2]*gPower_[13] + 13453440*g1Power_[4]*gPower_[11] +
        64576512*g1Power_[6]*gPower_[9] + 115315200*g1Power_[8]*gPower_[7] +
        6435*gPower_[15])) +
        8.74378692389149e-18*omegaPower_[16]*v1Power_[3]*(355687428095815.0*g1Power_[17]
        + 17*gPower_[1]*(69850115960625.0*g1Power_[1]*gPower_[15] +
        80918889891840000.0*g1Power_[11]*gPower_[5] +
        18673589975040000.0*g1Power_[13]*gPower_[3] +
        1422749712383992.0*g1Power_[15]*gPower_[1] +
        2980271614320000.0*g1Power_[3]*gPower_[13] +
        29206661820336000.0*g1Power_[5]*gPower_[11] +
        100137126241152000.0*g1Power_[7]*gPower_[9] +
        139079342001600000.0*g1Power_[9]*gPower_[7])) -
        4.17032762665731e-20*omegaPower_[17]*v1Power_[2]*(6402373705847470.0*g1Power_[18]
        + 51*gPower_[1]*(1502056893617280000.0*g1Power_[10]*gPower_[7] +
        728270009026558180.0*g1Power_[12]*gPower_[5] +
        144053408378885280.0*g1Power_[14]*gPower_[3] +
        9603560558585937.0*g1Power_[16]*gPower_[1] +
        3771906261873750.0*g1Power_[2]*gPower_[15] +
        80467333586640000.0*g1Power_[4]*gPower_[13] +
        525719912766048000.0*g1Power_[6]*gPower_[11] +
        1351851204255552000.0*g1Power_[8]*gPower_[9] +
        23283371986875.0*gPower_[17])) +
        3.56438258688659e-22*omegaPower_[18]*v1Power_[1]*(40548366800373008.0*g1Power_[19]
        + 19*gPower_[1]*(7520529151760625.0*g1Power_[1]*gPower_[17] +
        44105852421670993590.0*g1Power_[11]*gPower_[7] +
        18094708685814014660.0*g1Power_[13]*gPower_[5] +
        3101950060424664996.0*g1Power_[15]*gPower_[3] +
        182467650613731819.0*g1Power_[17]*gPower_[1] +
        406108574195073750.0*g1Power_[3]*gPower_[15] +
        5198189749696944000.0*g1Power_[5]*gPower_[13] +
        24258218831919072000.0*g1Power_[7]*gPower_[11] +
        48516437663838144000.0*g1Power_[9]*gPower_[9])) -
        4.57558740293529e-25*omegaPower_[19]*(810967336218252428.0*g1Power_[20]
        + 19*gPower_[1]*(1940657506553522492736.0*g1Power_[10]*gPower_[9] +
        1470195080722378921820.0*g1Power_[12]*gPower_[7] +
        516991676737504537760.0*g1Power_[14]*gPower_[5] +
        77548751510680883625.0*g1Power_[16]*gPower_[3] +
        4054836680258497880.0*g1Power_[18]*gPower_[1] +
        1504105830352125000.0*g1Power_[2]*gPower_[17] +
        40610857419507375000.0*g1Power_[4]*gPower_[15] +
        346545983313129600000.0*g1Power_[6]*gPower_[13] +
        1212910941595954110510.0*g1Power_[8]*gPower_[11] +
        7520529151760625.0*gPower_[19])) +
        342*omegaPower_[2]*v1Power_[17]*(3*g1Power_[1]*gPower_[2] +
        2*g1Power_[3])/35 - 646*omegaPower_[3]*v1Power_[16]*(8*g1Power_[4] +
        3*gPower_[1]*(8*g1Power_[2]*gPower_[1] + gPower_[3]))/105 +
        2584*omegaPower_[4]*v1Power_[15]*(8*g1Power_[5] +
        5*gPower_[1]*(3*g1Power_[1]*gPower_[3] + 8*g1Power_[3]*gPower_[1]))/231
        - 7752*omegaPower_[5]*v1Power_[14]*(16*g1Power_[6] +
        5*gPower_[1]*(18*g1Power_[2]*gPower_[3] + 24*g1Power_[4]*gPower_[1] +
        gPower_[5]))/1001 + 18088*omegaPower_[6]*v1Power_[13]*(16*g1Power_[7] +
        7*gPower_[1]*(5*g1Power_[1]*gPower_[5] + 30*g1Power_[3]*gPower_[3] +
        24*g1Power_[5]*gPower_[1]))/2145 -
        152*omegaPower_[7]*v1Power_[12]*(128*g1Power_[8] +
        7*gPower_[1]*(160*g1Power_[2]*gPower_[5] + 480*g1Power_[4]*gPower_[3] +
        256*g1Power_[6]*gPower_[1] + 5*gPower_[7]))/165 +
        36*omegaPower_[8]*v1Power_[11]*(128*g1Power_[9] +
        3*gPower_[1]*(105*g1Power_[1]*gPower_[7] + 1120*g1Power_[3]*gPower_[5] +
        2016*g1Power_[5]*gPower_[3] + 768*g1Power_[7]*gPower_[1]))/55 -
        4*omegaPower_[9]*v1Power_[10]*(256*g1Power_[10] +
        3*gPower_[1]*(1050*g1Power_[2]*gPower_[7] + 5600*g1Power_[4]*gPower_[5]
        + 6720*g1Power_[6]*gPower_[3] + 1920*g1Power_[8]*gPower_[1] +
        21*gPower_[9]))/21;
}


inline std::vector<std::function<double()> >
    HardSphereCollision1D::initializeI1Functions(int nMoments)
{

    if (nMoments <= 0) {
        std::string errorMessage = 
            "The parameter `nMoments` must be a positive integer.";
        throw std::runtime_error(errorMessage);
    }

    std::vector<std::function<double()> > functions;

    static_assert(nMomentsMax_ > 0);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<0>, this)
    );
    if (nMoments == 1)
        return functions;

    static_assert(nMomentsMax_ > 1);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<1>, this)
    );
    if (nMoments == 2)
        return functions;

    static_assert(nMomentsMax_ > 2);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<2>, this)
    );
    if (nMoments == 3)
        return functions;

    static_assert(nMomentsMax_ > 3);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<3>, this)
    );
    if (nMoments == 4)
        return functions;

    static_assert(nMomentsMax_ > 4);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<4>, this)
    );
    if (nMoments == 5)
        return functions;

    static_assert(nMomentsMax_ > 5);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<5>, this)
    );
    if (nMoments == 6)
        return functions;

    static_assert(nMomentsMax_ > 6);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<6>, this)
    );
    if (nMoments == 7)
        return functions;

    static_assert(nMomentsMax_ > 7);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<7>, this)
    );
    if (nMoments == 8)
        return functions;

    static_assert(nMomentsMax_ > 8);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<8>, this)
    );
    if (nMoments == 9)
        return functions;

    static_assert(nMomentsMax_ > 9);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<9>, this)
    );
    if (nMoments == 10)
        return functions;

    static_assert(nMomentsMax_ > 10);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<10>, this)
    );
    if (nMoments == 11)
        return functions;

    static_assert(nMomentsMax_ > 11);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<11>, this)
    );
    if (nMoments == 12)
        return functions;

    static_assert(nMomentsMax_ > 12);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<12>, this)
    );
    if (nMoments == 13)
        return functions;

    static_assert(nMomentsMax_ > 13);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<13>, this)
    );
    if (nMoments == 14)
        return functions;

    static_assert(nMomentsMax_ > 14);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<14>, this)
    );
    if (nMoments == 15)
        return functions;

    static_assert(nMomentsMax_ > 15);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<15>, this)
    );
    if (nMoments == 16)
        return functions;

    static_assert(nMomentsMax_ > 16);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<16>, this)
    );
    if (nMoments == 17)
        return functions;

    static_assert(nMomentsMax_ > 17);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<17>, this)
    );
    if (nMoments == 18)
        return functions;

    static_assert(nMomentsMax_ > 18);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<18>, this)
    );
    if (nMoments == 19)
        return functions;

    static_assert(nMomentsMax_ > 19);
    functions.push_back(
        std::bind(&HardSphereCollision1D::computeI1<19>, this)
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