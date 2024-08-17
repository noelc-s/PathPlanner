
#include "../inc/mpc.h"

// Constructor: Initialize the MPC parameters
MPC::MPC(const int nx, const int nu, const MPC_Params loaded_p) : nx_(nx), nu_(nu), mpc_params_(loaded_p) {
    nvar_ = nx_*mpc_params_.N+nu_*(mpc_params_.N-1);

    dynamics_A.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1),(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)));
    SparseIdentity.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1),(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)));
    dynamics_b_lb.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1)); // dynamics and torque bounds
    dynamics_b_ub.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1)); // dynamics and torque bounds
    dynamics_b_lb.setZero();
    dynamics_b_ub.setZero();

    f.resize(nvar_);
	full_ref.resize(nvar_);
	H.resize(nvar_,nvar_);

    Ad_.resize(nx_,nx_);
    Bd_.resize(nx_,nu_);
    Ad_ << 1, 0, mpc_params_.dt, 0,
          0, 1, 0, mpc_params_.dt,
          0, 0, 1, 0,
          0, 0, 0, 1;
    Bd_ <<  pow(mpc_params_.dt, 2)/2, 0,
            0, pow(mpc_params_.dt, 2)/2,
            mpc_params_.dt, 0,
            0, mpc_params_.dt;

    settings.verbose = false;
    settings.polish = true;
}

void MPC::buildDynamicEquality() {
    int offset = nx_;
    for (int i = 0; i < mpc_params_.N-1; i++) {
        for (int j = 0; j < nx_; j++) {
            for (int k = 0; k < nx_; k++) {
                dynamics_A.insert(offset+i * nx_ + j, i * nx_ + k) = Ad_(j,k);
            }
            for (int k = 0; k < nu_; k++) {
                dynamics_A.insert(offset+i*nx_+j,nx_*mpc_params_.N+i*nu_+k) = Bd_(j, k);
            }
        }
    }
    for (int i = offset; i < nx_*mpc_params_.N; i++) {
        SparseIdentity.insert(i,i) = 1;
    }
    dynamics_A = SparseIdentity - dynamics_A;
    // for initial condition
    for (int i = 0; i < offset; i++) {
        dynamics_A.insert(i,i) = 1;
    }
    // torque_bounds
    for (int i = 0; i < mpc_params_.N-1; i++) {
      for (int j = 0; j < nu_; j++) {
        dynamics_A.insert(nx_*mpc_params_.N+i*nu_+j,nx_*mpc_params_.N+i*nu_+j) = 1;
      }
    }
    // set torque limits
    for (int i = 0; i < mpc_params_.N-1; i++) {
      for (int j = 0; j < nu_; j++) {
        dynamics_b_lb(nx_*mpc_params_.N+i*nu_+j) = -mpc_params_.tau_max;
        dynamics_b_ub(nx_*mpc_params_.N+i*nu_+j) = mpc_params_.tau_max;
      }
    }
}

void MPC::buildConstraintInequality(const std::vector<matrix_t> A_constraint, const std::vector<vector_t> b_constraint) {
    const int num_constraints = A_constraint.size();

    int total_rows = std::accumulate(A_constraint.begin(), A_constraint.end(), 0, 
                                 [](int sum, const matrix_t& mat) { return sum + mat.rows(); });
    constraint_A.resize(4*total_rows, nx_*mpc_params_.N); // constrain x_k and x_kp1 except for last point
    constraint_b.resize(4*total_rows);
    constraint_A.setZero();
    constraint_b.setZero();
    int row = 0;
    for (int i = 0; i < num_constraints-1; i++) {
        // constrain over all control points
        // v_0 = x_k
        constraint_A.block(row, i*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(0,0,nx_,nx_);
        constraint_A.block(row, (i+1)*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(0,nx_,nx_,nx_);
        constraint_b.segment(row,A_constraint[i].rows()) = b_constraint[i];
        row += A_constraint[i].rows();

        // v_1
        constraint_A.block(row, i*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(nx_,0,nx_,nx_);
        constraint_A.block(row, (i+1)*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(nx_,nx_,nx_,nx_);
        constraint_b.segment(row,A_constraint[i].rows()) = b_constraint[i];
        row += A_constraint[i].rows();

        // v_2
        constraint_A.block(row, i*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(2*nx_,0,nx_,nx_);
        constraint_A.block(row, (i+1)*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(2*nx_,nx_,nx_,nx_);
        constraint_b.segment(row,A_constraint[i].rows()) = b_constraint[i];
        row += A_constraint[i].rows();

        // v_3 = x_kp1
        constraint_A.block(row, i*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(3*nx_,0,nx_,nx_);
        constraint_A.block(row, (i+1)*nx_, A_constraint[i].rows(), nx_) = A_constraint[i]*Bez_.block(3*nx_,nx_,nx_,nx_);
        constraint_b.segment(row,A_constraint[i].rows()) = b_constraint[i];
        row += A_constraint[i].rows();
    }
}

void MPC::setBez(matrix_t Bez) {
    Bez_ = Bez;
}

void MPC::buildCost()
{
    for (int i = 0; i < mpc_params_.N; i++) {
        for (int j = 0; j < nx_; j++) {
	  if (i == mpc_params_.N-1) {
            H.insert(i*nx_+j,i*nx_+j) = mpc_params_.terminalScaling*mpc_params_.stateScaling(j);
	  } else {
            H.insert(i*nx_+j,i*nx_+j) = mpc_params_.stateScaling(j);
	  }
        }
	if (i < mpc_params_.N-1) {
          for (int j = 0; j < nu_; j++) {
              H.insert(mpc_params_.N*nx_+i*nu_+j,mpc_params_.N*nx_+i*nu_+j) = mpc_params_.inputScaling(j);
          }
	}
    }
    f.setZero();
}

vector_t MPC::buildFromOptimalGraphSolve(const std::vector<Obstacle> obstacles,
                    const int num_adjacent_pts,
                    const int num_obstacle_faces,
                    const std::vector<vector_t> optimalSolutions, const std::vector<int> optimalInd,
                    const std::vector<vector_t> optimalPath,
                    const vector_t& xg)
{
    std::vector<matrix_t> A_constraint;
    std::vector<vector_t> b_constraint;

    vector_t sol(mpc_params_.N * nx_);
    for (int i = 0; i < std::min((size_t)mpc_params_.N,optimalInd.size()); i++) {
        // Have to traverse the list backwards
        const int vertex_ind = optimalInd[optimalInd.size()-1-i];
        vector_t x = optimalPath[optimalInd.size()-1-i];
        sol.segment(i*nx_,nx_) << x;
    }
    if (optimalInd.size() < mpc_params_.N) {
        for (int i = optimalInd.size(); i < mpc_params_.N; i++) 
        {    
            sol.segment(i*nx_,nx_) << xg;
        }
    }

    // for (int i = 0; i < optimalInd.size(); i++) {
    //     // Have to traverse the list backwards
    //     const int vertex_ind = optimalInd[optimalInd.size()-1-i];

    //     double min_violation = 1e5;
    //     matrix_t A_obs(1,4);
    //     vector_t b_obs(1);
    //     for (int i = 0; i < optimalSolutions.size(); i++) {
    //         vector_t optimal_solution = optimalSolutions[i];
    //         Obstacle obstacle = obstacles[i];
    //         vector_t violation = optimal_solution.segment(vertex_ind*(num_adjacent_pts+num_obstacle_faces)+num_adjacent_pts,num_obstacle_faces);
    //         Eigen::Index maxIndex;
    //         double violation_mag = violation.maxCoeff(&maxIndex);
    //         if (violation_mag < min_violation)
    //         {
    //             A_obs << obstacle.A.row(maxIndex);
    //             b_obs << obstacle.b.segment(maxIndex,1);
    //             min_violation = violation_mag;
    //         }
    //     }
    //     A_constraint.push_back(A_obs);
    //     b_constraint.push_back(b_obs);
    //     f_ref.segment(i*nx_,nx_) << optimalPath[optimalInd.size()-1-i];
    //     // f_ref.segment(i*nx_,nx_) << optimalPath.front();
    // }
    // matrix_t A_term(1,4);
    // vector_t b_term(1);
    // A_term << 0,0,0,0;
    // b_term << -1;
    // vector_t x_terminal;
    // x_terminal = optimalPath.front();
    // for (int i = optimalInd.size(); i < mpc_params_.N; i++) 
    // {    
    //     A_constraint.push_back(A_term);
    //     b_constraint.push_back(b_term);
    //     f_ref.segment(i*nx_,nx_) << x_terminal;
    // }
    // buildConstraintInequality(A_constraint, b_constraint);

    updateConstraintsSQP(obstacles, sol, xg);
    return sol;
}

void MPC::initialize()
{
    instance.objective_matrix = H;
    instance.objective_vector.resize(nvar_);
    instance.objective_vector << f;
    constraints.resize(dynamics_A.rows() + constraint_A.rows(), nvar_);
    lb.resize(dynamics_A.rows() + constraint_A.rows());
    ub.resize(dynamics_A.rows() + constraint_A.rows());
    lb.setZero();
    ub.setZero();
    lb << dynamics_b_lb, constraint_b;
    ub << dynamics_b_ub, std::numeric_limits<double>::infinity() * vector_t::Ones(constraint_b.size());

    for (int i = 0; i < dynamics_A.rows(); i++) {
        for (int j = 0; j < nvar_; j++) {
            constraints.insert(i, j) = dynamics_A.coeff(i, j);
        }
    }
    for (int i = 0; i < constraint_A.rows(); i++) {
        for (int j = 0; j < nx_*mpc_params_.N; j++) {
            constraints.insert(i + dynamics_A.rows(), j) = constraint_A(i, j);
        }
    }

    instance.constraint_matrix = constraints;
    instance.lower_bounds.resize(lb.size());
    instance.upper_bounds.resize(lb.size());
    instance.lower_bounds.setZero();
    instance.lower_bounds.setZero();
    instance.lower_bounds << lb;
    instance.upper_bounds << ub;

    auto status = solver.Init(instance, settings);

}

void MPC::updateConstraints(const vector_t& x0)
{
    dynamics_b_lb.segment(0,nx_) = x0;
    dynamics_b_ub.segment(0,nx_) = x0;

    lb.segment(0,nx_) = x0;
    ub.segment(0,nx_) = x0;
    lb.segment(dynamics_b_lb.size(), constraint_b.size()) << constraint_b;

    for (int i = 0; i < constraint_A.rows(); i++) {
        for (int j = 0; j < nx_*mpc_params_.N; j++) {
            constraints.coeffRef(i + dynamics_A.rows(), j) = constraint_A(i, j);
        }
    }
}

void MPC::updateConstraintsSQP(std::vector<Obstacle> obstacles, vector_t sol, const vector_t& xg) {
    std::vector<matrix_t> A_constraint;
    std::vector<vector_t> b_constraint;
    vector_t f_ref(nvar_);
    f_ref.setZero();
    
    for (int i = 0; i < mpc_params_.N; i++) {
        if (mpc_params_.use_previous_reference) {
            f_ref.segment(i*nx_,nx_) << sol.segment(i*nx_,nx_);
        } else {
            f_ref.segment(i*nx_,nx_) << xg;
        }
        matrix_t A(1,4);
        vector_t b(1);
        double max_viol = 1e3;
        for (auto obstacle : obstacles) {
            vector_t x(2);
            x << sol.segment(nx_*i,2);
            int closest_point = -1;
            double closest_dist = 1e3;
            for (int j = 0; j < obstacle.v.rows(); j++) {
                double dist_to_point = (x - obstacle.v.block(j,0,1,2).transpose()).norm();
                if (dist_to_point < closest_dist) {
                    closest_point = j;
                    closest_dist = dist_to_point;
                }
            }
            vector_t faces = obstacle.Adjacency.block(closest_point,0,1,obstacle.Adjacency.cols()).transpose();
            Eigen::Array<bool, Eigen::Dynamic, 1> inds = (obstacle.A.block(0,0,obstacle.A.rows(),2) * x - obstacle.b).array() > -1e-2 && faces.array() > 0;
            matrix_t A_hyp(1,2);
            for (int j = 0; j < inds.size(); j++) {
                if (inds(j) > 0) {
                    A_hyp += obstacle.A.block(j,0,1,2);
                }   
            }
            A_hyp = A_hyp / A_hyp.norm();
            vector_t b_hyp = A_hyp * obstacle.v.block(closest_point,0,1,2).transpose();

            double v = (A_hyp * x - b_hyp).maxCoeff();
            if (v < max_viol) {
                A << A_hyp,0,0;
                b << b_hyp;
                max_viol = v;
            }
        }
        A_constraint.push_back(A);
        b_constraint.push_back(b);
    }
    buildConstraintInequality(A_constraint, b_constraint);
    f = -H*f_ref;
}

void MPC::updateCost()
{
    auto status = solver.SetObjectiveVector(f);
}

// Solve the MPC problem
vector_t MPC::solve(std::vector<Obstacle> obstacles, vector_t sol, const vector_t& x0, const vector_t& xg) {

    vector_t mpc_sol = sol;
    for (int i = 0; i < mpc_params_.SQP_iters; i++) {
        updateConstraintsSQP(obstacles, mpc_sol, xg);
        updateCost();

        updateConstraints(x0);
        auto status1 = solver.UpdateConstraintMatrix(constraints);
        auto status2 = solver.SetBounds(lb, ub);
        
        OsqpExitCode exit_code = solver.Solve();
        mpc_sol = solver.primal_solution();
    }
    return mpc_sol;
    
}

void MPC::reset() {
    dynamics_A.setZero();
    dynamics_b_lb.setZero();
    dynamics_b_ub.setZero();
    H.setZero();
    f.setZero();
}