#include <Eigen/Dense>
#include <fstream>

#include <vector>
#include <stdexcept>

#include "../inc/Types.h"
#include "../inc/graph.h"
#include "../inc/mpc.h"
#include "../inc/utils.h"

#include <random>

int main()
{
    // Parameters
    const int num_pts = 2000;
    const int state_size = 4;
    const int num_obstacle_faces = 4;
    const int num_adjacent_pts = pow(2, 4);

    vector_4t starting_loc;
    vector_4t ending_loc;
    starting_loc << -3, -3, 0, 0;
    ending_loc << 3, 3, 0, 0;

    std::ofstream graph_file = open_log_file("../stored_graph.m");
    std::ofstream output_file = open_log_file("../output.m");

    std::vector<Obstacle> obstacles;
    Obstacle obstacle;
    std::vector<vector_t> b_obs;
    obstacle.A_obstacle.resize(num_obstacle_faces, state_size);
    obstacle.b_obstacle.resize(num_obstacle_faces);
    obstacle.A_obstacle << 1, 0, 0, 0,
        -1, 0, 0, 0,
        0, 1, 0, 0,
        0, -1, 0, 0;

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(-1, 1);
    std::uniform_real_distribution<> dis_y(-2, 2);

    double x_rand = dis_x(gen);
    double y_rand = dis_y(gen);
    obstacle.b_obstacle << 0.5 + x_rand, 0.5 - x_rand, 0.5 + y_rand, 0.5 - y_rand;
    b_obs.push_back(obstacle.b_obstacle);
    obstacles.push_back(obstacle);
    graph_file << "Obstacle_A{" << 1 << "}=[" << std::endl;
    graph_file << obstacle.A_obstacle << std::endl;
    graph_file << "];" << std::endl;
    graph_file << "Obstacle_b{" << 1 << "}=[" << std::endl;
    graph_file << obstacle.b_obstacle << std::endl;
    graph_file << "];" << std::endl;

    x_rand = dis_x(gen);
    y_rand = dis_y(gen);
    obstacle.b_obstacle << 0.5 + x_rand, 0.5 - x_rand, 0.5 + y_rand, 0.5 - y_rand;
    b_obs.push_back(obstacle.b_obstacle);
    obstacles.push_back(obstacle);
    graph_file << "Obstacle_A{" << 2 << "}=[" << std::endl;
    graph_file << obstacle.A_obstacle << std::endl;
    graph_file << "];" << std::endl;
    graph_file << "Obstacle_b{" << 2 << "}=[" << std::endl;
    graph_file << obstacle.b_obstacle << std::endl;
    graph_file << "];" << std::endl;

    x_rand = dis_x(gen);
    y_rand = dis_y(gen);
    obstacle.b_obstacle << 0.5 + x_rand, 0.5 - x_rand, 0.5 + y_rand, 0.5 - y_rand;
    b_obs.push_back(obstacle.b_obstacle);
    obstacles.push_back(obstacle);
    graph_file << "Obstacle_A{" << 3 << "}=[" << std::endl;
    graph_file << obstacle.A_obstacle << std::endl;
    graph_file << "];" << std::endl;
    graph_file << "Obstacle_b{" << 3 << "}=[" << std::endl;
    graph_file << obstacle.b_obstacle << std::endl;
    graph_file << "];" << std::endl;

    std::vector<vector_4t> points = generateUniformPoints(num_pts, -3, 3, -3, 3, 0, 0, 0, 0);

    Graph graph = buildGraph(points);
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

    Timer timer;
    GraphQP graphQP;

    timer.start();
    graphQP.setupQP(graphInstance, vertices, obstacles[0]);
    graphQP.initializeQP(graphSolver, graphInstance, graphSettings);

    MPC::MPC_Params mpc_params = loadParams();
    MPC mpc(4, 2, mpc_params);
    mpc.buildDynamicEquality();
    mpc.buildCost();

    const double num_traj = 100;
    for (int i = 0; i < num_traj; i++)
    {
        std::cout << "Trajectory percentage: " << (i+1)/num_traj << std::endl;
        // for (int o = 0; o < obstacles.size(); o++) {
        //     obstacles[o].b_obstacle << b_obs[o](0), b_obs[o](1),
        //         b_obs[o](2) + 1.0 * cos(2 * 3.14 * i / num_traj), b_obs[o](3) - 1.0 * cos(2 * 3.14 * i / num_traj);
        // }
        output_file << "Obs{" << i + 1 << "}=[" << std::endl;
        for (int o = 0; o < obstacles.size(); o++) {
            double x_add = sin(2 * 3.14 * i / num_traj);
            double y_add = cos(2 * 3.14 * i / num_traj);
            output_file << x_add << ", " << y_add << std::endl;
            obstacles[o].b_obstacle << b_obs[o](0) + x_add, b_obs[o](1) - x_add, b_obs[o](2) + y_add, b_obs[o](3) - y_add;
        }
        output_file << "];" << std::endl;

        
        Graph cut_graph;
        boost::copy_graph(graph, cut_graph);
        std::vector<vector_t> optimalSolutions;

        for (auto obstacle : obstacles)
        {
            graphQP.updateConstraints(graphSolver, obstacle, num_pts, num_obstacle_faces, vertices);
            graphQP.solveQP(graphSolver);

            double optimal_objective = graphSolver.objective_value();
            VectorXd optimal_solution = graphSolver.primal_solution();

            cutEdges(cut_graph, num_pts, num_adjacent_pts, obstacle.b_obstacle.size(), optimal_solution);
            optimalSolutions.push_back(optimal_solution);
        }

        std::vector<double> d(num_pts, std::numeric_limits<double>::max());
        std::vector<Vertex> p(num_vertices(cut_graph), graph_traits<Graph>::null_vertex()); // the predecessor array
        int starting_ind = -1;
        int ending_ind = -1;

        solveGraph(points, starting_loc, ending_loc, starting_ind, ending_ind, cut_graph, d, p);

        // Storing the edges is slow:

        // graph_traits<Graph>::edge_iterator ei, ei_end;
        // output_file << "Edges =[" << std::endl;
        // for (boost::tie(ei, ei_end) = edges(cut_graph); ei != ei_end; ++ei)
        // {
        //     output_file << source(*ei, cut_graph)
        //                 << "," << target(*ei, cut_graph) << std::endl;
        // }
        // output_file << "];" << std::endl;

        // output_file << "Sol = [" << std::endl;
        // output_file << optimal_solution << std::endl;
        // output_file << "];" << std::endl;

        output_file << "Path{" << i + 1 << "}=[" << std::endl;
        std::vector<int> optimalInd;
        std::vector<vector_t> optimalPath;
        for (Vertex v = vertex(ending_ind, cut_graph); v != vertex(starting_ind, cut_graph); v = p[v])
        {
            if (v < 0 || v > points.size())
            {
                output_file << "];" << std::endl;
                std::cout << "No optimal path found" << std::endl;
                return -1;
            }
            output_file << v + 1 << std::endl;
            optimalInd.push_back(v);
            optimalPath.push_back(points[v]);
        }
        output_file << "];" << std::endl;

        mpc.buildFromOptimalGraphSolve(obstacles, num_adjacent_pts,
                                       num_obstacle_faces, optimalSolutions, optimalInd, optimalPath);
        if (i == 0) {
            mpc.initialize();
        }
        mpc.updateConstraints(starting_loc);
        mpc.updateCost();

        vector_t sol;
        sol = mpc.solve(starting_loc);
        
        output_file << "MPC{" << i + 1 << "}=[" << std::endl;
        output_file << sol << std::endl;
        output_file << "];" << std::endl;
    }

    output_file.close();

    return 0;
}
