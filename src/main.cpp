#include <iostream>
#include "kernel.hpp"
#include "Types.h"
#include "utils.h"
#include <Eigen/Core>

// Main
int main(int argc, char ** argv)
{
    ObstacleCollector obstacles;
    std::vector<matrix_t> edges;
    matrix_t edge_(4,4);
    // edge_ << 0.5, 0.5, 0.5, 0.5,
    //         0, 0, 0, 0,
    //         0, 0, 0, 0,
    //         0, 0, 0, 0;
    // edges.push_back(edge_);
    // edge_ << 0, 0, 0, 0,
    //     0, 0, 0, 0,
    //     0, 0, 0, 0,
    //     0, 0, 0, 0;
    // edges.push_back(edge_);
    // edge_ << 0.2, 0, 0.5, 0,
    //     0, 0, 0.3, 0,
    //     0, 0, 0, 0,
    //     0, 0, 0, 0;
    edge_.setZero();
    for (int i = 0; i < 8; i++) {
        edges.push_back(edge_);
    }
    int_vector_t membership(edges.size());
    for (int i = 0; i < edges.size(); i++) {
        membership[i] = 3;
    }
            
    std::cout << obstacles.obstacles[0].A * edges[0].block(0,0,4,1) - obstacles.obstacles[0].b << std::endl;
    std::cout << std::endl;



    int num_edges = edges.size();
    int num_obstacles = obstacles.obstacles.size();
    // Prepare edge data
    double obstacle_A_flat[num_obstacles * 16];
    double obstacle_b_flat[num_obstacles * 4];
    double obstacle_Adj_flat[num_obstacles * 16];
    double obstacle_v_flat[num_obstacles * 8];
    double edges_flat[num_edges * 16];
    int d_member[num_obstacles * num_edges];


    for (int i = 0; i < num_obstacles * num_edges; i++) {
        d_member[i] = 3;
    }


    for (int o = 0; o < num_obstacles; o++) {
        // copy obstacle A
        for (int col = 0; col < 4; col++)
        {
            for (int row = 0; row < 4; row++)
            {
                obstacle_A_flat[o * 16 + row + col*4] = obstacles.obstacles[o].A(row, col);
                obstacle_Adj_flat[o * 16 + row + col*4] = obstacles.obstacles[o].Adjacency(row, col);
            }
            obstacle_b_flat[o * 4 + col] = obstacles.obstacles[o].b(col);
        }
        for (int col = 0; col < 2; col++)
        {
            for (int row = 0; row < 4; row++)
            {
                obstacle_v_flat[o * 8 + row + col*4] = obstacles.obstacles[o].v(row, col);
            }
        }
    }
    // copy edges
    for (int i = 0; i < num_edges; i++)
    {
        Eigen::MatrixXd mat = edges[i];
        for (int col = 0; col < 4; col++)
        {
            for (int row = 0; row < 4; row++)
            {
                edges_flat[i * 16 + row + col*4] = mat(row, col);
            }
        }
    }




    Timer timer(true);
    for (int i = 0; i < 10; i ++) {
        std::cout << i << std::endl;
        timer.start();
        // Kernel::GraphQP_ObstacleMembershipHeuristicSharedMemory(obstacle_A_flat, obstacle_b_flat, obstacle_Adj_flat, obstacle_v_flat, edges_flat, d_member, num_edges, num_obstacles);
        // memcpy(membership.data(), d_member, num_obstacles * num_edges * sizeof(int));

        Kernel::GraphQP_ObstacleMembershipHeuristicFlat(obstacle_A_flat, obstacle_b_flat, obstacle_Adj_flat, obstacle_v_flat, edges_flat, membership, num_edges, num_obstacles);

        // Kernel::GraphQP_ObstacleMembershipHeuristic(obstacles.obstacles, edges, membership);
        timer.time("GPU Call: ");
    }

    for (int i = 0; i < membership.size(); i++)
        std::cout << membership[i] << std::endl;

    // for (int i = 0; i < 100; i++) {
    //     timer.start();
    //     int_vector_t Membership(obstacles.obstacles.size() * edges.size());
    //     timer.time("Dynamically allocate membership: ");
    // }

	return 0;
}