#pragma once
#ifndef KERNEL_H
#define KERNEL_H

#include <vector>
#include "obstacle.h"
#include <Eigen/Core>

namespace Kernel
{
    void GraphQP_ObstacleMembershipHeuristic(std::vector<Obstacle> obstacles, const std::vector<Eigen::MatrixXd>& edges, int_vector_t& member);
    void MPC_GetActiveConstraints(std::vector<Obstacle> obstacles, const vector_t &sol, vector_t &A1, vector_t &A2, vector_t &b, vector_t &dist);
}

#endif