#include <Eigen/Dense>
#include <fstream>

#include <vector>
#include <stdexcept>

#include "../inc/Types.h"
#include "../inc/graph.h"
#include "../inc/mpc.h"
#include "../inc/utils.h"

#include "Bezier.h"

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
    std::uniform_real_distribution<> dist_freq(1,3);
    std::vector<double> freq;

    double x_rand = dis_x(gen);
    double y_rand = dis_y(gen);
    obstacle.b_obstacle << 0.5 + x_rand, 0.5 - x_rand, 0.5 + y_rand, 0.5 - y_rand;
    b_obs.push_back(obstacle.b_obstacle);
    obstacles.push_back(obstacle);
    freq.push_back(dist_freq(gen));
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
    freq.push_back(dist_freq(gen));
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
    freq.push_back(dist_freq(gen));
    graph_file << "Obstacle_A{" << 3 << "}=[" << std::endl;
    graph_file << obstacle.A_obstacle << std::endl;
    graph_file << "];" << std::endl;
    graph_file << "Obstacle_b{" << 3 << "}=[" << std::endl;
    graph_file << obstacle.b_obstacle << std::endl;
    graph_file << "];" << std::endl;

    std::vector<vector_4t> points = generateUniformPoints(num_pts, -3, 3, -3, 3, 0, 0, 0, 0);

    // Setup the OSQP instance
    OsqpInstance graphInstance;
    OsqpSolver graphSolver;
    OsqpSettings graphSettings;
    graphSettings.verbose = false;
    graphSettings.polish = false;

    Timer timer;
    GraphQP graphQP;

    MPC::MPC_Params mpc_params = loadParams();
    MPC mpc(4, 2, mpc_params);
    mpc.buildDynamicEquality();
    mpc.buildCost();

    int order = 3;
    int gamma = 2;
    double tau = mpc.p.dt;
    int m = 2;
    Bezier B = Bezier(order,gamma,tau);
    matrix_t H_vec = B.H_vec(B.H, m, order, gamma, gamma-1);
    matrix_t D_nT_vec = B.inv_DT_vec(m, order, gamma);
    matrix_t Bez = H_vec * D_nT_vec;
    mpc.setBez(Bez);

    // Right now this is done independently for x and y!
    matrix_t A_x, H, xbar, f_xbar, g_xbar;
    vector_t b_x, K;
    std::vector<matrix_t> Q;
    int m_ind = 1;
    scalar_t Lf = 0;
    scalar_t Lg = 1;
    scalar_t e_bar = 0;
    scalar_t u_max = 5;
    A_x.resize(4,2);
    A_x << 1,     0,
            -1,     0,
            0 ,    1,
            0 ,   -1;
    b_x.resize(4);
    // TODO: change this when the state space changes
    b_x << 3,3,3,3;
    H = B.H_matrix(order);
    xbar.resize(2,1);
    xbar << 0,0;
    f_xbar.resize(1,1);
    f_xbar << 0;
    g_xbar.resize(1,1);
    g_xbar << 1;
    Q.push_back(matrix_t::Identity(4,4));
    K.resize(2);
    K << -1,-1;
    matrix_t F, G;
    B.F_G(A_x, b_x, H, m_ind, xbar, f_xbar, g_xbar, gamma, Q, 
                    Lg, Lf, e_bar, K, u_max,
                    F, G);
    auto F_G = [&B, A_x, b_x, H, m_ind,gamma, Q, Lg, Lf, e_bar, K, u_max]
                (matrix_t xbar, matrix_t f_xbar, matrix_t g_xbar, matrix_t &F, matrix_t &G)
                {B.F_G(A_x, b_x, H, m_ind, xbar, f_xbar, g_xbar, gamma, Q, Lg, Lf, e_bar, K, u_max, F, G);};

    // Graph graph = buildGraph(points);
    Graph graph = buildGraph(points, F_G);
    std::vector<matrix_t> vertices = getVerticesOfBezPoly(points);

    // Log graph
    log(points, graph_file, "Points");
    log(vertices, graph_file, "ReachableVertices");

    graphQP.setupQP(graphInstance, vertices, obstacles[0]);
    graphQP.initializeQP(graphSolver, graphInstance, graphSettings);

    const double num_traj = 100;
    for (int i = 0; i < num_traj; i++)
    {
        std::cout << "Trajectory percentage: " << (i+1)/num_traj << std::endl;
        output_file << "Obs{" << i + 1 << "}=[" << std::endl;
        for (int o = 0; o < obstacles.size(); o++) {
            double x_add = 0*cos(2 * 3.14 * freq[o] * i / num_traj);
            double y_add = 2*sin(2 * 3.14 * freq[o] * i / num_traj);
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

            cutEdges(cut_graph, num_pts, vertices, obstacle.b_obstacle.size(), optimal_solution);
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
        // output_file << optimalSolutions[0] << std::endl;
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

        // starting_loc << sol.segment(4,4);
    }

    output_file.close();

    return 0;
}
