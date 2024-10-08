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

#include <thread>
#include <condition_variable>

#include <random>

int main()
{
    // Parameters
    Params params;
    MPC_Params mpc_params;
    Planner_Params planner_params;
    loadPlannerParams("../config/params.yaml", params, mpc_params, planner_params);
    std::ofstream graph_file = open_log_file("../stored_graph.m");
    std::ofstream output_file = open_log_file("../output.m");

    const int state_size = 4;
    const int input_size = 2;

    vector_t starting_loc;
    vector_t ending_loc;
    starting_loc.resize(4);
    ending_loc.resize(4);
    starting_loc <<-3.8, -3.8 , 0, 0;
    ending_loc << 3.4 ,3.4, 0, 0;

    ObstacleCollector O = ObstacleCollector();
    logObstacles(O.obstacles, graph_file);

    Timer timer(true);

    PathPlanner planner(state_size, input_size, mpc_params, planner_params);    
    timer.start();
    planner.initialize(O);
    timer.time("Initialization: ");

    // Log graph
    log(planner.points, graph_file, "Points");
    log(planner.edges, graph_file, "EdgeControlPoints");

    std::condition_variable cv;
    std::mutex m;

    std::cout << "Number of edges: " << planner.edges.size() << std::endl;

    // std::thread cutGraph(static_cast<void (PathPlanner::*)(ObstacleCollector&, std::ofstream&, std::condition_variable&, std::mutex&)>(&PathPlanner::cutGraph),
    //                     &planner, std::ref(O), std::ref(output_file), std::ref(cv), std::ref(m));

    //                     sleep(1);

    for (int i = 0; i < params.num_traj; i++)
    {
        std::cout << "Trajectory percentage: " << (float)(i + 1) / params.num_traj << std::endl;
        // O.updateObstaclePositions((float)i / params.num_traj);
        logObstaclePosition(O.obstacles, output_file, i);

        timer.start();
        planner.cutGraph(O, output_file, cv, m);
        timer.time("Cut: ");

        std::vector<int> optimalInd;
        std::vector<vector_t> optimalPath;
        timer.start();
        planner.findPath(O.obstacles, starting_loc, ending_loc, optimalInd, optimalPath, cv, m);
        timer.time("Find Path: ");
        log(optimalInd, output_file, "Path{" + (std::to_string)(i + 1) + "}");

        vector_t sol, graph_sol;
        graph_sol.resize(mpc_params.N * planner.mpc_->nx_);
        timer.start();
        planner.refineWithMPC(graph_sol, sol, O, optimalInd, optimalPath, starting_loc, ending_loc);       
        timer.time("Refine: ");
        log(sol, output_file, "MPC{" + (std::to_string)(i + 1) + "}");

        starting_loc << sol.segment(4,4);

        std::cout << std::endl;
    }

    output_file.close();

    return 0;
}
