#include "../inc/obstacle.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

ObstacleCollector::ObstacleCollector() {
    Obstacle obstacle;
    const int num_obstacle_faces = 4;
    obstacle.A.resize(num_obstacle_faces, 4);
    obstacle.v.resize(num_obstacle_faces, 2);
    obstacle.center.resize(2);
    obstacle.b.resize(num_obstacle_faces);
    obstacle.A << 1, 0, 0, 0,
        -1, 0, 0, 0,
        0, 1, 0, 0,
        0, -1, 0, 0;
    obstacle.Adjacency.resize(num_obstacle_faces, num_obstacle_faces);

    obstacle.Adjacency << 0, 1, 0, 1,
        0, 1, 1, 0,
        1, 0, 1, 0,
        1, 0, 0, 1;
    int rows = 20;                // Number of rows
    int cols = 30;                // Number of columns
    double cellWidth = 0.4;       // Width of each cell in meters
    double cellHeight = 0.4;      // Height of each cell in meters
    double wallThickness = 0.05;   // Thickness of walls in meters
    double density = 0.5;         // Proportion of walls that are active
    Maze maze(rows, cols, cellWidth, cellHeight, wallThickness);
    auto vertices = maze.computeWallVertices();
    for (int i = 0; i < vertices.size(); i++) {
        obstacle.v = vertices[i];
        obstacle.v.col(0) += -6*vector_t::Ones(4);
        obstacle.v.col(1) += -4*vector_t::Ones(4);
        obstacle.b << obstacle.v(3, 0), -obstacle.v(0, 0), obstacle.v(1, 1), -obstacle.v(0, 1);
        obstacle.center << 0,0;
        obstacles.push_back(obstacle);
    }

    // obstacle.Adjacency << 1, 0, 1, 0,
    //     1, 0, 0, 1,
    //     0, 1, 0, 1,
    //     0, 1, 1, 0;
    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::uniform_real_distribution<> dis_x(-2, 2);
    // std::uniform_real_distribution<> dis_y(-2, 2);
    // std::uniform_real_distribution<> dist_freq(1, 3);

    // scalar_t x_rand = dis_x(gen);
    // scalar_t y_rand = dis_y(gen);
    // obstacle.center << 0.5,0;
    // // obstacle.center << x_rand, y_rand;
    // scalar_t obstacle_size;


    // for (int i = 0; i < 19; i++){
    // x_rand = dis_x(gen);
    // y_rand = dis_y(gen);
    // obstacle.center << x_rand + sgn(x_rand)*0.3, y_rand + sgn(y_rand)*0.3;
    // // obstacle.center << -0.5, 0;
    // obstacle_size = 0.1;
    // obstacle.b << obstacle_size + obstacle.center(0), obstacle_size - obstacle.center(0), obstacle_size + obstacle.center(1), obstacle_size - obstacle.center(1);
    // obstacle.v << obstacle_size + obstacle.center(0), obstacle_size + obstacle.center(1),
    //     obstacle_size + obstacle.center(0), -obstacle_size + obstacle.center(1),
    //     -obstacle_size + obstacle.center(0), -obstacle_size + obstacle.center(1),
    //     -obstacle_size + obstacle.center(0), obstacle_size + obstacle.center(1);
    // obstacle.center << 0,0; // center should be difference in motion from IC
    // b_obs.push_back(obstacle.b);
    // obstacles.push_back(obstacle);
    // freq.push_back(dist_freq(gen));
    // }

    // obstacle.center << 0.5,0;
    // obstacle_size = 0.1;
    // obstacle.b << obstacle_size + obstacle.center(0), obstacle_size - obstacle.center(0), obstacle_size + obstacle.center(1), obstacle_size - obstacle.center(1);
    // obstacle.v << obstacle_size + obstacle.center(0), obstacle_size + obstacle.center(1),
    //     obstacle_size + obstacle.center(0), -obstacle_size + obstacle.center(1),
    //     -obstacle_size + obstacle.center(0), -obstacle_size + obstacle.center(1),
    //     -obstacle_size + obstacle.center(0), obstacle_size + obstacle.center(1);
    // obstacle.center << 0,0; // center should be difference in motion from IC
    // b_obs.push_back(obstacle.b);
    // obstacles.push_back(obstacle);
    // freq.push_back(dist_freq(gen));

}

void ObstacleCollector::updateObstaclePositions(scalar_t t) {
    for (int o = 0; o < obstacles.size(); o++)
    {
        scalar_t x_add = 0 * cos(2 * 3.14 * freq[o] * t);
        scalar_t y_add = 2 * sin(2 * 3.14 * freq[o] * t);
        obstacles[o].center << x_add, y_add;
        obstacles[o].b << b_obs[o](0) + obstacles[o].center(0), b_obs[o](1) - obstacles[o].center(0), b_obs[o](2) + obstacles[o].center(1), b_obs[o](3) - obstacles[o].center(1);
        obstacles[o].v << obstacles[o].b(0), obstacles[o].b(2),
            obstacles[o].b(0), -obstacles[o].b(3),
            -obstacles[o].b(1), -obstacles[o].b(3),
            -obstacles[o].b(1), obstacles[o].b(2);
    }
}

void ObstacleCollector::updateObstaclePositions(int o, scalar_t x, scalar_t y) {

        obstacles[o].center << x, y;
        obstacles[o].b << b_obs[o](0) + obstacles[o].center(0), b_obs[o](1) - obstacles[o].center(0), b_obs[o](2) + obstacles[o].center(1), b_obs[o](3) - obstacles[o].center(1);
        obstacles[o].v << obstacles[o].b(0), obstacles[o].b(2),
            obstacles[o].b(0), -obstacles[o].b(3),
            -obstacles[o].b(1), -obstacles[o].b(3),
            -obstacles[o].b(1), obstacles[o].b(2);
}

void getSeparatingHyperplane(Obstacle obstacle, vector_t x, vector_t &A_hyp, scalar_t &b_hyp)
{
    int closest_point = -1;
    scalar_t closest_dist = 1e3;
    scalar_t dist_to_point;
    Eigen::Array<bool, Eigen::Dynamic, 1> inds;

    // (obstacle.v.block(0,0,obstacle.v.rows(),2).rowwise() - x.transpose()).rowwise().squaredNorm().minCoeff(&closest_point);

    for (int j = 0; j < obstacle.v.rows(); j++) {
        dist_to_point = (x - obstacle.v.block(j,0,1,2).transpose()).squaredNorm();
        if (dist_to_point < closest_dist) {
            closest_point = j;
            closest_dist = dist_to_point;
        }
    }
    vector_t faces = obstacle.Adjacency.block(closest_point,0,1,obstacle.Adjacency.cols()).transpose();
    inds = (obstacle.A.block(0,0,obstacle.A.rows(),2) * x - obstacle.b).array() > -1e-2 && faces.array() > 0;

    int num_constraints_violated = inds.cast<int>().sum();
    if (num_constraints_violated == 1) {
        for (int j = 0; j < inds.size(); j++) {
            if (inds(j) > 0) {
                A_hyp = obstacle.A.block(j,0,1,2).transpose();
            }   
        }
    } else {
        A_hyp = x - obstacle.v.block(closest_point,0,1,2).transpose();
    }
    A_hyp = A_hyp / A_hyp.norm();
    b_hyp = (A_hyp.transpose() * obstacle.v.block(closest_point,0,1,2).transpose()).value();
}
