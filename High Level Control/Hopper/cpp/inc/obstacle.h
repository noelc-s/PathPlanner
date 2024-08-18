#pragma once
#include "Types.h"
#include <random>

struct Obs
{
    matrix_t A;
    vector_t b;
    vector_t center;
    matrix_t v; // vertex representation
    matrix_t Adjacency; // points to faces (1 at (i,j) if point i touches face j)
};

class Obstacle {
public:
    Obstacle();

    void updateObstaclePositions(double t);
    
    std::vector<Obs> obstacles;
    std::vector<vector_t> b_obs;
    std::vector<double> freq;
};