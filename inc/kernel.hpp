#pragma once
#ifndef KERNEL_H
#define KERNEL_H

#include <vector>
#include "obstacle.h"
#include <Eigen/Core>

namespace Kernel
{
    void GraphQP_ObstacleMembershipHeuristic(Obstacle obstacle, const std::vector<Eigen::MatrixXd>& edges, int_vector_t& member);
}

#endif