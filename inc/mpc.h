#pragma once
#include "Types.h"
#include "utils.h"
#include "obstacle.h"
#include "osqp++.h"

#include <Eigen/Dense>
#include <numeric>

using namespace osqp;
using ::Eigen::SparseMatrix;
using ::Eigen::Triplet;

class MPC {
public:

    MPC_Params mpc_params_;

    SparseMatrix<double> constraints;
    Eigen::VectorXd lb, ub;
    bool isInitialized = false;
    int nx_, nu_, nvar_;

    MPC(const int nx, const int nu, const MPC_Params loaded_p, const matrix_t &Bez);
    void buildDynamicEquality();
    void buildConstraintInequality(const std::vector<matrix_t> A_constraint, const std::vector<vector_t> b_constraint);
    void buildCost();
    void initialize();
    void updateConstraints(const vector_t& x0);
    void updateConstraintsSQP(ObstacleCollector O, vector_t sol, const vector_t& xg);
    void updateCost();
    vector_t solve(ObstacleCollector O, vector_t sol, const vector_t& x0, const vector_t& xg); 
    void reset();

private:
    matrix_t Ad_;         // State transition matrix
    matrix_t Bd_;         // Control input matrix
    matrix_t Bez_;
    SparseMatrix<double> SparseIdentity;
    SparseMatrix<double> dynamics_A;
    matrix_t constraint_A;
    Eigen::VectorXd dynamics_b_lb, dynamics_b_ub, constraint_b;
    SparseMatrix<double> H;
    matrix_t H_cost;
    Eigen::VectorXd f, full_ref;

    OsqpSolver solver;
    OsqpInstance instance;
    OsqpSettings settings;
};   