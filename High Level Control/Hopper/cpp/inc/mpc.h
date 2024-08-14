
#include "Types.h"
#include "osqp++.h"

#include <Eigen/Dense>
#include <numeric>

using namespace osqp;
using ::Eigen::SparseMatrix;
using ::Eigen::Triplet;

class MPC {
public:

    // Parameters for the MPC program
    struct MPC_Params {
        int N;
        double dt;
	    int SQP_iter;
        vector_t stateScaling;
        vector_t inputScaling;
        double terminalScaling;
        double tau_max;
    } p;

    SparseMatrix<double> constraints;
    vector_t lb, ub;

    MPC(const int nx, const int nu, const MPC_Params loaded_p);
    void setBez(matrix_t Bez);
    void buildDynamicEquality();
    void buildConstraintInequality(const std::vector<matrix_t> A_constraint, const std::vector<vector_t> b_constraint);
    void buildCost();
    void buildFromOptimalGraphSolve(const std::vector<Obstacle> obstacles, const int num_adjacent_pts,
                         const int num_obstacle_faces, const std::vector<vector_t> optimal_solution,
                         const std::vector<int> optimalInd, const std::vector<vector_t> optimalPath);
    void initialize();
    void updateConstraints(const vector_t& x0);
    void updateCost();
    vector_t solve(const vector_t& x0); 
    void reset();

private:
    int nx_, nu_, nvar_;
    matrix_t Ad_;         // State transition matrix
    matrix_t Bd_;         // Control input matrix
    matrix_t Bez_;
    SparseMatrix<double> SparseIdentity;
    SparseMatrix<double> dynamics_A;
    matrix_t constraint_A;
    vector_t dynamics_b_lb, dynamics_b_ub, constraint_b;
    SparseMatrix<double> H;
    vector_t f, full_ref;

    OsqpSolver solver;
    OsqpInstance instance;
    OsqpSettings settings;
};   

MPC::MPC_Params loadParams();