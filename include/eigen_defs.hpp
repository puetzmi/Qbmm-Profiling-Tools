/**
 * @file eigen_defs.hpp
 * @author M. Puetz
 * @brief Some convenient type definitions for use with the Eigen library.
 * @date 2022-09-20
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#ifndef EIGEN_DEFS_HPP
#define EIGEN_DEFS_HPP

#include "global_defs.hpp"
#include <eigen3/Eigen/Core>

using Vector = Eigen::Matrix<double, Eigen::Dynamic, 1>;
using VectorMap = Eigen::Map<Vector, MALLOC_ALIGN>;
using Matrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MatrixMap = Eigen::Map<Matrix, MALLOC_ALIGN>;
using SelfAdjointViewRowMajor = Eigen::SelfAdjointView<Matrix, Eigen::UpLoType::Upper>;

#endif // EIGEN_DEFS_HPP