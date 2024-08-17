#pragma once
#include <chrono>
#include <iostream>
#include <fstream>
#include "yaml-cpp/yaml.h"
#include "Types.h"

class Timer
{
public: 
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    std::chrono::duration<double, std::nano> duration;

    void start();
    void time();
    void time(std::string info);
};

struct Params
{
    int num_traj;
    int num_points;
    bool log_edges;
};

struct MPC_Params 
{
    int N;
    double dt;
    int SQP_iters;
    vector_t stateScaling;
    vector_t inputScaling;
    double terminalScaling;
    double tau_max;
    bool use_previous_reference;
};

void loadParams(std::string filename, Params &p, MPC_Params &mpc_p);

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