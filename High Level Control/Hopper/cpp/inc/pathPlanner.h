#include "utils.h"
#include "mpc.h"
#include "Bezier.h"
#include "graph.h"
#include <random>

using FunctionType = std::function<void(const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&, int,
                                         const Eigen::MatrixXd&, const Eigen::MatrixXd&, const Eigen::MatrixXd&,
                                         int, const std::vector<Eigen::MatrixXd>&, double, double, double,
                                         const Eigen::MatrixXd&, double, Eigen::MatrixXd&, Eigen::MatrixXd&)>;

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
    void initialize(const Obstacle O);

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

    void F_G(const matrix_t& xbar, const matrix_t& f_xbar, const matrix_t& g_xbar, matrix_t& F, matrix_t& G);

    void cutGraph(Obstacle O, std::ofstream &output_file);
    void findPath(vector_t starting_location, vector_t ending_location, std::vector<int> &optimalInd, std::vector<vector_t> &optimalPath);
    void refineWithMPC(vector_t &sol, Obstacle O, std::vector<int> optimalInd, std::vector<vector_t> optimalPath, vector_t starting_loc, vector_t ending_loc);
};

std::vector<vector_4t> generateUniformPoints(int n, double min_x, double max_x, double min_y, double max_y,
                                                    double min_dx, double max_dx, double min_dy, double max_dy);

