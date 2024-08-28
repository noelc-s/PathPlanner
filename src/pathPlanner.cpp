#include "../inc/pathPlanner.h"
#include "../inc/kernel.hpp"

PathPlanner::PathPlanner(int state_size, int input_size, const MPC_Params mpc_params, const Planner_Params p_p)
    : state_size_(state_size), params_(p_p)
{

    int order = 3;
    B_p.gamma = 2;
    scalar_t tau = mpc_params.dt;
    B = std::make_unique<Bezier>(order, B_p.gamma, tau);
    matrix_t H_vec = B->H_vec(B->H, input_size, order, B_p.gamma, B_p.gamma - 1);
    matrix_t D_nT_vec = B->inv_DT_vec(input_size, order, B_p.gamma);
    Bez_ = H_vec * D_nT_vec;

    // Right now this is done independently for x and y!

    B_p.m_ind = 1;
    B_p.Lf = 0;
    B_p.Lg = 1;
    B_p.e_bar = 0;
    B_p.u_max = mpc_params.tau_max;

    B_p.A_x.resize(4, 2);
    B_p.A_x << 1, 0,
        -1, 0,
        0, 1,
        0, -1;
    B_p.b_x.resize(4);
    // TODO: change this when the state space changes
    B_p.b_x << 3, 3, 3, 3;
    B_p.H = B->H_matrix(order);
    B_p.Q.push_back(matrix_t::Identity(4, 4));
    B_p.K.resize(2);
    B_p.K << -1, -1;

    // xbar.resize(2, 1);
    // xbar << 0, 0;
    // f_xbar.resize(1, 1);
    // f_xbar << 0;
    // g_xbar.resize(1, 1);
    // g_xbar << 1;
    // matrix_t F, G;

    mpc_ = std::make_unique<MPC>(state_size, input_size, mpc_params, Bez_);
    mpc_->buildDynamicEquality();
    mpc_->buildCost();

    D_nT_vec_ = B->inv_DT_vec(B_p.m_ind, order, B_p.gamma);
}

void PathPlanner::F_G(const matrix_t &xbar, const matrix_t &f_xbar, const matrix_t &g_xbar, matrix_t &F, matrix_t &G)
{
    B->F_G(B_p.A_x, B_p.b_x, B_p.H, B_p.m_ind, xbar, f_xbar, g_xbar, B_p.gamma, B_p.Q, B_p.Lg, B_p.Lf, B_p.e_bar, B_p.K, B_p.u_max, F, G);
}

void PathPlanner::initialize(const ObstacleCollector O)
{

    graphSettings.verbose = false;
    graphSettings.polish = false;

    points = generateUniformPoints(
        params_.num_points,
        params_.x_bounds[0], params_.x_bounds[1],
        params_.x_bounds[2], params_.x_bounds[3],
        params_.dx_bounds[0], params_.dx_bounds[1],
        params_.dx_bounds[2], params_.dx_bounds[3]);
    
    // points = generateGridPoints(
    //     params_.num_points,
    //     params_.x_bounds[0], params_.x_bounds[1],
    //     params_.x_bounds[2], params_.x_bounds[3],
    //     params_.dx_bounds[0], params_.dx_bounds[1],
    //     params_.dx_bounds[2], params_.dx_bounds[3]);

    auto F_G_bound = std::bind(&PathPlanner::F_G, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);

    graph = buildGraph(points, F_G_bound, D_nT_vec_, Bez_);

    edges = getBezEdges(graph, vertexInds);

    graphQP.setupQP(graphInstance, edges, O.obstacles[0]);
    graphQP.initializeQP(graphSolver, graphInstance, graphSettings);
}

std::vector<vector_4t> generateUniformPoints(int n, scalar_t min_x, scalar_t max_x, scalar_t min_y, scalar_t max_y,
                                             scalar_t min_dx, scalar_t max_dx, scalar_t min_dy, scalar_t max_dy)
{
    std::vector<vector_4t> points;
    points.reserve(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(min_x, max_x);
    std::uniform_real_distribution<> dis_y(min_y, max_y);
    std::uniform_real_distribution<> dis_dx(min_dx, max_dx);
    std::uniform_real_distribution<> dis_dy(min_dy, max_dy);

    for (int i = 0; i < n; ++i)
    {
        vector_4t point;
        point << dis_x(gen), dis_y(gen), dis_dx(gen), dis_dy(gen);
        points.push_back(point);
    }

    return points;
}

std::vector<vector_4t> generateGridPoints(int n, scalar_t min_x, scalar_t max_x, scalar_t min_y, scalar_t max_y,
                                             scalar_t min_dx, scalar_t max_dx, scalar_t min_dy, scalar_t max_dy)
{
    std::vector<vector_4t> points;
    points.reserve(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(min_x, max_x);
    std::uniform_real_distribution<> dis_y(min_y, max_y);
    std::uniform_real_distribution<> dis_dx(min_dx, max_dx);
    std::uniform_real_distribution<> dis_dy(min_dy, max_dy);

    for (int i = 0; i < n; ++i)
    {
        vector_4t point;
        point << i%10 * (max_x - min_x) / 10 + min_x, i/10 * (max_y - min_y) / 10 + min_y, 0, 0;
        points.push_back(point);
    }

    return points;
}

void PathPlanner::cutGraph(ObstacleCollector O)
{
    cut_graph.clear();
    boost::copy_graph(graph, cut_graph);
    Timer timer(PRINT_TIMING);
    int_vector_t membership(edges.size());

    for (auto obstacle : O.obstacles)
    {
        timer.start();
        // graphQP.ObstacleMembershipHeuristic(obstacle, edges, membership);
        Kernel::GraphQP_ObstacleMembershipHeuristic(obstacle, edges, membership);
        timer.time("    Heuristic check: ");
        //////////////////////////////////////////////
        VectorXd optimal_solution;
        if ((membership.array() == 2).any()) {
            std::vector<matrix_t> uncertain_edges;
            for (size_t i = 0; i < membership.size(); ++i) {
                if (membership[i] == 2) {
                    uncertain_edges.push_back(edges[i]);
                }
            }
            timer.time("    Copy Edges: ");
            graphQP.setupQP(graphInstance, uncertain_edges, obstacle);
            graphQP.initializeQP(graphSolver, graphInstance, graphSettings);

            // graphQP.updateConstraints(graphSolver, obstacle, edges, membership);
            // graphQP.updateConstraints(graphSolver, obstacle, edges);
            //////////////////////////////////////////////
            timer.time("    Update constraints: ");
            graphQP.solveQP(graphSolver);
            timer.time("    Solve: ");

            optimal_solution = graphSolver.primal_solution();
            percentageEdgesRemovedWithHeuristic = 1. - (float)uncertain_edges.size() / edges.size();
        } else {
            optimal_solution.resize(0);
            percentageEdgesRemovedWithHeuristic = 1;
        }
        

        timer.start();
        //////////////////////////////////////////////
        // cutGraphEdges(cut_graph, edges, vertexInds, obstacle.b.size(), optimal_solution);
        cutGraphEdges(cut_graph, edges, vertexInds, obstacle.b.size(), optimal_solution, membership);
        //////////////////////////////////////////////
        timer.time("    Cut graph edges: ");
         
        // std::cout << "Percentage of nodes removed with heuristic: " <<  << std::endl;
    }
}

void PathPlanner::cutGraph(ObstacleCollector O, std::ofstream &output_file)
{
    cutGraph(O);
    if (params_.log_edges)
    {
        logEdges(cut_graph, output_file, "Edges");
    }
}

void PathPlanner::findPath(vector_t starting_location, vector_t ending_location, std::vector<int> &optimalInd, std::vector<vector_t> &optimalPath)
{
    std::vector<scalar_t> d(num_vertices(cut_graph), std::numeric_limits<scalar_t>::max());
    std::vector<Vertex> p(num_vertices(cut_graph), graph_traits<Graph>::null_vertex()); // the predecessor array
    int starting_ind = -1;
    int ending_ind = -1;

    solveGraph(points, starting_location, ending_location, starting_ind, ending_ind, cut_graph, d, p);

    for (Vertex v = vertex(ending_ind, cut_graph); v != vertex(starting_ind, cut_graph); v = p[v])
    {
        if (v < 0 || v > points.size())
        {
            optimalPathFound = 0;
            return;
        }
        optimalInd.push_back(v);
        optimalPath.push_back(points[v]);
    }
    optimalPathFound = 1;
    optimalInd.push_back(vertex(starting_ind, cut_graph));
    optimalPath.push_back(points[vertex(starting_ind, cut_graph)]);
}

void PathPlanner::refineWithMPC(vector_t &sol, ObstacleCollector O, std::vector<int> optimalInd, std::vector<vector_t> optimalPath, vector_t starting_loc, vector_t ending_loc) {
        Timer timer(PRINT_TIMING);
        timer.start();
        sol = mpc_->buildFromOptimalGraphSolve(O, optimalInd, optimalPath, ending_loc);
        timer.time("    Build Graph from Optimal Solve: ");
        mpc_->updateConstraints(starting_loc);
        timer.time("    Update Constraints: ");
        mpc_->updateCost();
        timer.time("    Update Cost: ");

        sol = mpc_->solve(O, sol, starting_loc, ending_loc);
        timer.time("    Solve: ");
}