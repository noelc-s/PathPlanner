#include <Eigen/Dense>
#include <fstream>

#include <vector>
#include <stdexcept>

#include "../inc/Types.h"
#include "../inc/graph.h"
#include "../inc/obstacle.h"
#include "../inc/mpc.h"
#include "../inc/utils.h"
#include "../inc/pathPlanner.h"

#include <random>

int main()
{
    // Parameters
    Params params;
    MPC_Params mpc_params;
    Planner_Params planner_params;
    loadParams("../config/params.yaml", params, mpc_params, planner_params);
    std::ofstream graph_file = open_log_file("../stored_graph.m");
    std::ofstream output_file = open_log_file("../output.m");

    const int state_size = 4;
    const int input_size = 2;

    vector_4t starting_loc;
    vector_4t ending_loc;
    starting_loc << -3, -3, 0, 0;
    ending_loc << 3, 3, 0, 0;

    Obstacle O = Obstacle();
    logObstacles(O.obstacles, graph_file);

    PathPlanner planner(state_size, input_size, mpc_params, planner_params);    
    planner.initialize(O);

    // Log graph
    log(planner.points, graph_file, "Points");
    log(planner.edges, graph_file, "EdgeControlPoints");

    for (int i = 0; i < params.num_traj; i++)
    {
        std::cout << "Trajectory percentage: " << (float)(i + 1) / params.num_traj << std::endl;
        O.updateObstaclePositions((float)i / params.num_traj);
        logObstaclePosition(O.obstacles, output_file, i);

        planner.cutGraph(O, output_file);

        std::vector<int> optimalInd;
        std::vector<vector_t> optimalPath;
        planner.findPath(starting_loc, ending_loc, optimalInd, optimalPath);
        log(optimalInd, output_file, "Path{" + (std::to_string)(i + 1) + "}");

        vector_t sol;
        planner.refineWithMPC(sol, O, optimalInd, optimalPath, starting_loc, ending_loc);       
        log(sol, output_file, "MPC{" + (std::to_string)(i + 1) + "}");
    }

    output_file.close();

    return 0;
}
