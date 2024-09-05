#pragma once
#include <chrono>
#include <iostream>
#include <fstream>
#include "yaml-cpp/yaml.h"
#include "Types.h"
#include "obstacle.h"

#define PRINT_TIMING false

class Timer
{
public: 
    Timer(bool print);
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    std::chrono::duration<scalar_t, std::nano> duration;

    void start();
    scalar_t time();
    scalar_t time(std::string info);

    bool print_;
};

struct Planner_Params
{
    int num_points;
    vector_t x_bounds;
    vector_t dx_bounds;
    bool log_edges;
    bool use_planner;
    scalar_t buffer;
};

struct Params
{
    int num_traj;
};

struct MPC_Params 
{
    int N;
    scalar_t dt;
    int SQP_iters;
    vector_t stateScaling;
    vector_t inputScaling;
    scalar_t terminalScaling;
    scalar_t tau_max;
    scalar_t vel_max;
    bool use_previous_reference;
    scalar_t buffer;
};

void loadPlannerParams(std::string filename, Params &p, MPC_Params &mpc_p, Planner_Params &p_p);

// ChatGPT with some crazy templating
template <typename Derived>
void log(std::vector<Derived> data, std::ofstream& file, std::string name,
    typename std::enable_if<std::is_base_of<Eigen::MatrixBase<Derived>, Derived>::value>::type* = 0)
{
    file << name << " = [" << std::endl;
    for (auto d : data)
    {
        file << d.transpose() << std::endl;
    }
    file << "];" << std::endl;
    std::cout << "Done Logging " << name << std::endl;
}

// ChatGPT with some crazy templating
template <typename T>
void log(std::vector<T> data, std::ofstream& file, std::string name,
    typename std::enable_if<!std::is_base_of<Eigen::MatrixBase<T>, T>::value>::type* = 0)
{
    file << name << " = [" << std::endl;
    for (auto d : data)
    {
        file << d << std::endl;
    }
    file << "];" << std::endl;
    std::cout << "Done Logging " << name << std::endl;
}

void log(matrix_t data, std::ofstream& file, std::string name);
void logEdges(Graph cut_graph, std::ofstream& file, std::string name);
void logObstacles(const std::vector<Obstacle> obstacles, std::ofstream& file);
void logObstaclePosition(const std::vector<Obstacle> obstacles, std::ofstream& file, const int i);
