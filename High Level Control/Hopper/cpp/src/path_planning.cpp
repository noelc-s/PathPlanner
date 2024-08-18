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

    PathPlanner planner(state_size, input_size, mpc_params, planner_params);    

    vector_4t starting_loc;
    vector_4t ending_loc;
    starting_loc << -3, -3, 0, 0;
    ending_loc << 3, 3, 0, 0;

    Obstacle O = Obstacle();
    logObstacles(O.obstacles, graph_file);

    planner.initialize(O);

    Timer timer;

    // Log graph
    log(planner.points, graph_file, "Points");
    log(planner.edges, graph_file, "EdgeControlPoints");

    for (int i = 0; i < params.num_traj; i++)
    {
        std::cout << "Trajectory percentage: " << (float)(i + 1) / params.num_traj << std::endl;
        O.updateObstaclePositions((float)i / params.num_traj);
        logObstaclePosition(O.obstacles, output_file, i);

        Graph cut_graph;
        boost::copy_graph(planner.graph, cut_graph);
        std::vector<vector_t> optimalSolutions;

        for (auto obstacle : O.obstacles)
        {
            planner.graphQP.updateConstraints(planner.graphSolver, obstacle, planner.edges);
            planner.graphQP.solveQP(planner.graphSolver);

            double optimal_objective = planner.graphSolver.objective_value();
            VectorXd optimal_solution = planner.graphSolver.primal_solution();

            cutEdges(cut_graph, planner.edges, planner.vertexInds, obstacle.b.size(), 
            optimal_solution);
            optimalSolutions.push_back(optimal_solution);
        }

        std::vector<double> d(num_vertices(cut_graph), std::numeric_limits<double>::max());
        std::vector<Vertex> p(num_vertices(cut_graph), graph_traits<Graph>::null_vertex()); // the predecessor array
        int starting_ind = -1;
        int ending_ind = -1;

        solveGraph(planner.points, starting_loc, ending_loc, starting_ind, ending_ind, cut_graph, d, p);

        if (params.log_edges) {
            logEdges(cut_graph, output_file, "Edges");
        }

        std::vector<int> optimalInd;
        std::vector<vector_t> optimalPath;
        for (Vertex v = vertex(ending_ind, cut_graph); v != vertex(starting_ind, cut_graph); v = p[v])
        {
            if (v < 0 || v > planner.points.size())
            {
                std::cout << "No optimal path found" << std::endl;
                return -1;
            }
            optimalInd.push_back(v);
            optimalPath.push_back(planner.points[v]);
        }
        optimalInd.push_back(vertex(starting_ind, cut_graph));
        optimalPath.push_back(planner.points[vertex(starting_ind, cut_graph)]);
        log(optimalInd, output_file, "Path{" + (std::to_string)(i + 1) + "}");

        vector_t sol = planner.mpc_->buildFromOptimalGraphSolve(O, optimalSolutions, optimalInd, optimalPath, ending_loc);
        if (i == 0)
        {
            planner.mpc_->initialize();
        }
        planner.mpc_->updateConstraints(starting_loc);
        planner.mpc_->updateCost();

        sol = planner.mpc_->solve(O, sol, starting_loc, ending_loc);

        log(sol, output_file, "MPC{" + (std::to_string)(i + 1) + "}");

        // starting_loc << sol.segment(4,4);
    }

    output_file.close();

    return 0;
}
