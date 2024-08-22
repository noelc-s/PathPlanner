#include "../inc/obstacle.h"

Obstacle::Obstacle() {
    Obs obstacle;
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
    obstacle.Adjacency << 1, 0, 1, 0,
        1, 0, 0, 1,
        0, 1, 0, 1,
        0, 1, 1, 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(-1, 1);
    std::uniform_real_distribution<> dis_y(-2, 2);
    std::uniform_real_distribution<> dist_freq(1, 3);

    double x_rand = dis_x(gen);
    double y_rand = dis_y(gen);
    obstacle.center << 0.5, 0;
    // obstacle.center << x_rand, y_rand;
    obstacle.b << 0.15 + obstacle.center(0), 0.15 - obstacle.center(0), 0.15 + obstacle.center(1), 0.15 - obstacle.center(1);
    obstacle.v << 0.15 + obstacle.center(0), 0.15 + obstacle.center(1),
        0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
        -0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
        -0.15 + obstacle.center(0), 0.15 + obstacle.center(1);
    obstacle.center << 0,0; // center should be difference in motion from IC
    b_obs.push_back(obstacle.b);
    obstacles.push_back(obstacle);
    freq.push_back(dist_freq(gen));

    x_rand = dis_x(gen);
    y_rand = dis_y(gen);
    // obstacle.center << x_rand, y_rand;
    obstacle.center << -0.5, 0;
    obstacle.b << 0.15 + obstacle.center(0), 0.15 - obstacle.center(0), 0.15 + obstacle.center(1), 0.15 - obstacle.center(1);
    obstacle.v << 0.15 + obstacle.center(0), 0.15 + obstacle.center(1),
        0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
        -0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
        -0.15 + obstacle.center(0), 0.15 + obstacle.center(1);
    obstacle.center << 0,0; // center should be difference in motion from IC
    b_obs.push_back(obstacle.b);
    obstacles.push_back(obstacle);
    freq.push_back(dist_freq(gen));

    // x_rand = dis_x(gen);
    // y_rand = dis_y(gen);
    // obstacle.center << 0.3, -0.3;
    // obstacle.b << 0.15 + obstacle.center(0), 0.15 - obstacle.center(0), 0.15 + obstacle.center(1), 0.15 - obstacle.center(1);
    // obstacle.v << 0.15 + obstacle.center(0), 0.15 + obstacle.center(1),
    //     0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
    //     -0.15 + obstacle.center(0), -0.15 + obstacle.center(1),
    //     -0.15 + obstacle.center(0), 0.15 + obstacle.center(1);
    // b_obs.push_back(obstacle.b);
    // obstacles.push_back(obstacle);
    freq.push_back(dist_freq(gen));
}

void Obstacle::updateObstaclePositions(double t) {
    for (int o = 0; o < obstacles.size(); o++)
    {
        double x_add = 0 * cos(2 * 3.14 * freq[o] * t);
        double y_add = 2 * sin(2 * 3.14 * freq[o] * t);
        obstacles[o].center << x_add, y_add;
        obstacles[o].b << b_obs[o](0) + obstacles[o].center(0), b_obs[o](1) - obstacles[o].center(0), b_obs[o](2) + obstacles[o].center(1), b_obs[o](3) - obstacles[o].center(1);
        obstacles[o].v << obstacles[o].b(0), obstacles[o].b(2),
            obstacles[o].b(0), -obstacles[o].b(3),
            -obstacles[o].b(1), -obstacles[o].b(3),
            -obstacles[o].b(1), obstacles[o].b(2);
    }
}

void Obstacle::updateObstaclePositions(int o, double x, double y) {

        obstacles[o].center << x, y;
        obstacles[o].b << b_obs[o](0) + obstacles[o].center(0), b_obs[o](1) - obstacles[o].center(0), b_obs[o](2) + obstacles[o].center(1), b_obs[o](3) - obstacles[o].center(1);
        obstacles[o].v << obstacles[o].b(0), obstacles[o].b(2),
            obstacles[o].b(0), -obstacles[o].b(3),
            -obstacles[o].b(1), -obstacles[o].b(3),
            -obstacles[o].b(1), obstacles[o].b(2);
}

