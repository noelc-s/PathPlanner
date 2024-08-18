#include "../inc/utils.h"


void Timer::start() {
    start_time = std::chrono::high_resolution_clock::now();
}

void Timer::time() {
    end_time = std::chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    std::cout << duration.count()*1e-6 << " ms" << std::endl;
    start_time = end_time;
}

void Timer::time(std::string info) {
    std::cout << info;
    time();
}

void loadParams(std::string filename, Params &p, MPC_Params &mpc_p, Planner_Params &p_p) {
    YAML::Node config = YAML::LoadFile("../config/params.yaml");
    p.num_traj = config["num_traj"].as<int>();
    p.log_edges = config["log_edges"].as<bool>();

    mpc_p.N = config["MPC"]["N"].as<int>();
    mpc_p.dt = config["MPC"]["dt"].as<double>();
    mpc_p.SQP_iters = config["MPC"]["SQP_iters"].as<int>();
    mpc_p.terminalScaling = config["MPC"]["terminalScaling"].as<double>();
    mpc_p.tau_max = config["MPC"]["tau_max"].as<double>();
    mpc_p.use_previous_reference = config["MPC"]["use_previous_reference"].as<bool>();
    mpc_p.stateScaling.resize(4);
    mpc_p.inputScaling.resize(2);
    auto tmp = config["MPC"]["stateScaling"].as<std::vector<double>>();
    for (int i = 0; i < mpc_p.stateScaling.size(); i++)
        mpc_p.stateScaling(i) = tmp[i];
    tmp = config["MPC"]["inputScaling"].as<std::vector<double>>();
    for (int i = 0; i < mpc_p.inputScaling.size(); i++)
        mpc_p.inputScaling(i) = tmp[i];

    p_p.x_bounds.resize(4);
    p_p.dx_bounds.resize(4);
    tmp = config["Planner"]["x_bounds"].as<std::vector<double>>();
    for (int i = 0; i < p_p.x_bounds.size(); i++)
        p_p.x_bounds(i) = tmp[i];
    tmp = config["Planner"]["dx_bounds"].as<std::vector<double>>();
    for (int i = 0; i < p_p.dx_bounds.size(); i++)
        p_p.dx_bounds(i) = tmp[i];
    p_p.num_points = config["Planner"]["num_points"].as<int>();

}

void log(matrix_t data, std::ofstream& file, std::string name)
{
    file << name << " = [" << std::endl;
    file << data << std::endl;
    file << "];" << std::endl;
    std::cout << "Done Logging " << name << std::endl;
}

void logEdges(Graph cut_graph, std::ofstream& file, std::string name) {
    boost::graph_traits<Graph>::edge_iterator ei, ei_end;
    file << name << " = [" << std::endl;
    for (boost::tie(ei, ei_end) = edges(cut_graph); ei != ei_end; ++ei)
    {
        file << source(*ei, cut_graph)
                    << "," << target(*ei, cut_graph) << std::endl;
    }
    file << "];" << std::endl;
    std::cout << "Done Logging " << name << std::endl;
}

void logObstacles(const std::vector<Obs> obstacles, std::ofstream& file) {
    for (int i = 0; i < obstacles.size(); i++) {
        file << "Obstacle_A{" << i+1 << "}=[" << std::endl;
        file << obstacles[i].A << std::endl;
        file << "];" << std::endl;
        file << "Obstacle_b{" << i+1 << "}=[" << std::endl;
        file << obstacles[i].b << std::endl;
        file << "];" << std::endl;
    }
    std::cout << "Done Logging obstacles" << std::endl;
}

void logObstaclePosition(const std::vector<Obs> obstacles, std::ofstream& file, const int i) {
    file << "Obs{" << i + 1 << "}=[" << std::endl;
    for (int o = 0; o < obstacles.size(); o++)
    {
        file << obstacles[o].center(0) << ", " << obstacles[o].center(1) << std::endl;
    }
    file << "];" << std::endl;
}
