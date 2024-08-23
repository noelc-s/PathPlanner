#pragma once
#include "Types.h"
#include <random>

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

    void updateObstaclePositions(double t);
    void updateObstaclePositions(int o, double x, double y);
    
    std::vector<Obstacle> obstacles;
    std::vector<vector_t> b_obs;
    std::vector<double> freq;
};

void getSeparatingHyperplane(Obstacle obstacle, vector_t x, matrix_t &A_hyp, vector_t &b_hyp);