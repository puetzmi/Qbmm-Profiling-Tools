/**
 * @file test_cholesky.cpp
 * @author M. Puetz
 * @brief Test of build ID assignment by CMake from Git version.
 * @date 2022-11-18
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <cstdio>
#include "build_id.hpp"

#ifdef NDEBUG
    #undef NDEBUG
#endif
#include <cassert>



int main()
{

    std::printf("BUILD_ID: %s\n", CMAKE_GIT_BUILD_ID);

    assert(sizeof(CMAKE_GIT_BUILD_ID) > 1);
    assert(CMAKE_GIT_BUILD_ID[0] != '@');

    return 0;
}
