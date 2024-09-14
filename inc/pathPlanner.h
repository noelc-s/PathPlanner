#include "utils.h"
#include "mpc.h"
#include "Bezier.h"
#include "graph.h"
#include <random>
#include <condition_variable>

using FunctionType = std::function<void(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, int,
                                         const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                                         int, const std::vector<Eigen::MatrixXd>&, scalar_t, scalar_t, scalar_t,
                                         const Eigen::MatrixXd&, scalar_t, Eigen::MatrixXd&, Eigen::MatrixXd&)>;

struct BezierParams
{
    matrix_t A_x, H, xbar, f_xbar, g_xbar;
    vector_t b_x, K;
    std::vector<matrix_t> Q;
    int m_ind;
    scalar_t Lf;
    scalar_t Lg;
    scalar_t e_bar;
    scalar_t u_max;
    int gamma;
};

class PathPlanner {
public:
    PathPlanner(int state_size, int input_size, const MPC_Params mpc_params, const Planner_Params p_p);
    void initialize(const ObstacleCollector O);

    Planner_Params params_;

    int state_size_;
    std::unique_ptr<MPC> mpc_;
    Graph graph;
    Graph cut_graph;
    GraphQP graphQP;
    std::vector<matrix_t> edges;
    std::vector<std::pair<int,int>> vertexInds;
    std::vector<vector_4t> points;

    FunctionType F_G_;
    matrix_t D_nT_vec_;
    matrix_t Bez_;
    std::unique_ptr<Bezier> B;

    // Setup the OSQP instance
    OsqpInstance graphInstance;
    OsqpSolver graphSolver;
    OsqpSettings graphSettings;

    BezierParams B_p;

    scalar_t percentageEdgesRemovedWithHeuristic = 0;
    bool optimalPathFound = false;

    void F_G(const matrix_t& xbar, const matrix_t& f_xbar, const matrix_t& g_xbar, matrix_t& F, matrix_t& G);

    void cutGraph(ObstacleCollector &O, std::condition_variable &cv, std::mutex &m);
    void cutGraph(ObstacleCollector &O, std::ofstream &output_file, std::condition_variable &cv, std::mutex &m);
    void cutGraphLoop(ObstacleCollector &O, std::ofstream &output_file, double &timing, std::condition_variable &cv, std::mutex &m);
    void findPath(const std::vector<Obstacle> obstacles, vector_t starting_location, vector_t &ending_location, std::vector<int> &optimalInd, std::vector<vector_t> &optimalPath, std::condition_variable &cv, std::mutex &m);
    void refineWithMPC(vector_t &graph_sol, vector_t &sol, ObstacleCollector O, std::vector<int> optimalInd, std::vector<vector_t> optimalPath, vector_t starting_loc, vector_t ending_loc);
};

std::vector<vector_4t> generateUniformPoints(int n, scalar_t min_x, scalar_t max_x, scalar_t min_y, scalar_t max_y,
                                                    scalar_t min_dx, scalar_t max_dx, scalar_t min_dy, scalar_t max_dy);

std::vector<vector_4t> generateGridPoints(int n, scalar_t min_x, scalar_t max_x, scalar_t min_y, scalar_t max_y,
                                             scalar_t min_dx, scalar_t max_dx, scalar_t min_dy, scalar_t max_dy);

