# pragma once
#include <Eigen/Dense>

using vector_2t = Eigen::Matrix<double, 2, 1>;
using vector_4t = Eigen::Matrix<double, 4, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

struct Obstacle
{
    matrix_t A_obstacle;
    vector_t b_obstacle;
};