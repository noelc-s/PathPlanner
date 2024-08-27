#include <iostream>
#include "kernel.hpp"
#include "Types.h"
#include <Eigen/Core>

// Main
int main(int argc, char ** argv)
{
    ObstacleCollector obstacles;
    std::vector<matrix_t> edges;
    matrix_t edge_(4,4);
    std::vector<int> membership;
    edge_ << 0.5, 0.5, 0.5, 0.5,
            0, 0, 0, 0,
            0, 0, 0, 0,
            0, 0, 0, 0;
    edges.push_back(edge_);
    edge_ << 0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    edges.push_back(edge_);
    edge_ << 0.2, 0, 0.5, 0,
        0, 0, 0.3, 0,
        0, 0, 0, 0,
        0, 0, 0, 0;
    for (int i = 0; i < 40000; i++) {
        edges.push_back(edge_);
    }
    for (int i = 0; i < edges.size(); i++) {
        membership.push_back(3);
    }
            
    // std::cout << obstacles.obstacles[0].A * edges[0].block(0,0,4,1) - obstacles.obstacles[0].b << std::endl;
    // std::cout << std::endl;

    for (int i = 0; i < 10000; i ++) {
        std::cout << i << std::endl;
        Kernel::GraphQP_ObstacleMembershipHeuristic(obstacles.obstacles[0], edges, membership);
    }

    // for (auto mem : membership)
    //     std::cout << mem << std::endl;

	return 0;
}