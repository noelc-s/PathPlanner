#include "../inc/pathPlanner.h"


PathPlanner::PathPlanner(int state_size, int input_size, const MPC_Params mpc_params, const Planner_Params p_p)
         : state_size_(state_size), params_(p_p) {
    
    int order = 3;
    B_p.gamma = 2;
    double tau = mpc_params.dt;
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

void PathPlanner::F_G(const matrix_t& xbar, const matrix_t& f_xbar, const matrix_t& g_xbar, matrix_t& F, matrix_t& G) {
    B->F_G(B_p.A_x, B_p.b_x, B_p.H, B_p.m_ind, xbar, f_xbar, g_xbar, B_p.gamma, B_p.Q, B_p.Lg, B_p.Lf, B_p.e_bar, B_p.K, B_p.u_max, F, G);
}

void PathPlanner::initialize(const Obstacle O) {

    graphSettings.verbose = false;
    graphSettings.polish = false;

    points = generateUniformPoints(
            params_.num_points,
            params_.x_bounds[0], params_.x_bounds[1],
            params_.x_bounds[2], params_.x_bounds[3],
            params_.dx_bounds[0], params_.dx_bounds[1],
            params_.dx_bounds[2], params_.dx_bounds[3]);

    auto F_G_bound = std::bind(&PathPlanner::F_G, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4, std::placeholders::_5);

    graph = buildGraph(points, F_G_bound, D_nT_vec_, Bez_);

    edges = getBezEdges(graph, vertexInds);

    graphQP.setupQP(graphInstance, edges, O.obstacles[0]);
    graphQP.initializeQP(graphSolver, graphInstance, graphSettings);
}

std::vector<vector_4t> generateUniformPoints(int n, double min_x, double max_x, double min_y, double max_y,
                                                    double min_dx, double max_dx, double min_dy, double max_dy)
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