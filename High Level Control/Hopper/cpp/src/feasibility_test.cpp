#include <Eigen/Dense>
#include <fstream>

#include <vector>
#include <chrono>
#include <stdexcept>

#include "../inc/Types.h"
#include "../inc/graph.h"
#include "../inc/mpc.h"

int main()
{
    // Parameters
    const int num_pts = 500;
    const int state_size = 2;
    const int num_obstacle_faces = 4;

    vector_2t starting_loc(-3,-3);
    vector_2t ending_loc(3,3);

    std::ofstream graph_file = open_log_file("../stored_graph.m");
    std::ofstream output_file = open_log_file("../output.m");

    Obstacle obstacle;
    obstacle.A_obstacle.resize(num_obstacle_faces, state_size);
    obstacle.b_obstacle.resize(num_obstacle_faces);
    obstacle.A_obstacle << 1,  0,
                          -1,  0,
                           0,  1,
                           0, -1;
    obstacle.b_obstacle << 1.0, 1.0, 1.0, 1.0;
    
    std::vector<vector_2t> points = generateUniformPoints(num_pts, -3, 3, -3, 3);

    Graph graph = buildGraph(points, obstacle);
    std::vector<matrix_t> vertices = getReachableVertices(points);
    // Log graph
    log(points, graph_file, "Points");
    log(vertices, graph_file, "ReachableVertices");

    // Setup the OSQP instance
    OsqpInstance graphInstance;
    OsqpSolver graphSolver;
    OsqpSettings graphSettings;
    graphSettings.verbose = false;
    graphSettings.polish = false;

    const double num_traj = 100;
    for (int i = 0; i <= num_traj; i++)
    {
        obstacle.b_obstacle << 1.0 - 1.5 * sin(2*3.14 * i / num_traj), 1.0 + 1.5 * sin(2*3.14 * i / num_traj), 
                               1.0 + 1.5 * cos(2*3.14 * i / num_traj), 1.0 - 1.5 * cos(2*3.14 * i / num_traj);
    
        // auto start = std::chrono::high_resolution_clock::now();

    setupQP(graphInstance, vertices, obstacle);    
    initializeQP(graphSolver, graphInstance, graphSettings);
    solveQP(graphSolver);

    double optimal_objective = graphSolver.objective_value();
    VectorXd optimal_solution = graphSolver.primal_solution();

    Graph cut_graph;
    boost::copy_graph(graph, cut_graph);
    cutEdges(cut_graph, num_pts, obstacle.b_obstacle.size(), optimal_solution);

    std::vector<double> d(num_pts, std::numeric_limits<double>::max());
    std::vector<Vertex> p(num_vertices(cut_graph), graph_traits<Graph>::null_vertex()); //the predecessor array
    int starting_ind = -1;
    int ending_ind = -1;

    solveGraph(points, starting_loc, ending_loc, starting_ind, ending_ind, cut_graph, d, p);
    
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double, std::nano> duration = end - start;
    // std::cout << duration.count()*1e-6 << std::endl;

    // graph_traits<Graph>::edge_iterator ei, ei_end;
    // output_file << "Edges =[" << std::endl;
    // for (boost::tie(ei, ei_end) = edges(cut_graph); ei != ei_end; ++ei)
    // {
    //     output_file << source(*ei, cut_graph)
    //             << "," << target(*ei, cut_graph) << std::endl;
    // }
    // output_file << "];" << std::endl;

    // output_file << "Sol = [" << std::endl;
    // output_file << optimal_solution << std::endl;
    // output_file << "];" << std::endl;

    output_file << "Path{" << i+1 << "}=[" << std::endl;
    std::vector<int> optimalInd;
    std::vector<vector_t> optimalPath;
    for (Vertex v = vertex(ending_ind, cut_graph); v != vertex(starting_ind, cut_graph); v = p[v]) {
        output_file << v+1 << std::endl;
        optimalInd.push_back(v);
        optimalPath.push_back(points[v]);
    }
    output_file << "];" << std::endl;

    MPC::MPC_Params mpc_params = loadParams();
    MPC mpc(4, 2, mpc_params);
    mpc.buildDynamicEquality();
    mpc.buildCost();
    mpc.buildFromOptimalGraphSolve(obstacle, optimal_solution, optimalInd, optimalPath);
    
    vector_t x0(4);
    x0 << -3,-3, 0, 0;
    vector_t sol = mpc.solve(x0);
    output_file << "MPC{" << i+1 << "}=[" << std::endl;
    output_file << sol << std::endl;
    output_file << "];" << std::endl;

    }

    output_file.close();

    return 0;
}
