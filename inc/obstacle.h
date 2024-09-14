#pragma once
#include "Types.h"
#include <random>
#include "maze.hpp"

struct Obstacle
{
    matrix_t A;
    vector_t b;
    vector_t center;
    matrix_t v; // vertex representation
    matrix_t Adjacency; // points to faces (1 at (i,j) if point i touches face j)
};

class ObstacleCollector {
public:
    ObstacleCollector();

    void updateObstaclePositions(scalar_t t);
    void updateObstaclePositions(int o, scalar_t x, scalar_t y);
    
    std::vector<Obstacle> obstacles;
    std::vector<vector_t> b_obs;
    std::vector<scalar_t> freq;
};

void getSeparatingHyperplane(Obstacle obstacle, vector_t x, vector_t &A_hyp, scalar_t &b_hyp);
