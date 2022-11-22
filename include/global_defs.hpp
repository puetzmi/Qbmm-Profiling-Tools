/**
 * @file global_defs.hpp
 * @author Michele Puetz
 * @brief Essential definitions required in several files.
 * @version 0.1
 * @date 2022-05-10
 * 
 * @copyright Copyright (c) 2022
 * 
 */


#ifndef DEFS_HPP
#define DEFS_HPP


#ifdef MALLOC_ALIGN_SIZE
    constexpr int MALLOC_ALIGN = MALLOC_ALIGN_SIZE;
#elif __AVX2__
    constexpr int MALLOC_ALIGN = 32;
#else
    constexpr int MALLOC_ALIGN = 16;
#endif

static_assert(MALLOC_ALIGN == 16 || MALLOC_ALIGN == 32 
            || MALLOC_ALIGN == 64 || MALLOC_ALIGN == 128);


constexpr int DEFAULT_RANDOM_SEED = 125125;

#endif // DEFS_HPP