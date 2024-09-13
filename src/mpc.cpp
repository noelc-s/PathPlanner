
#include "../inc/mpc.h"

// Constructor: Initialize the MPC parameters
MPC::MPC(const int nx, const int nu, const MPC_Params loaded_p, const matrix_t &Bez)
     : nx_(nx), nu_(nu), mpc_params_(loaded_p), Bez_(Bez) {
    nvar_ = nx_*mpc_params_.N+nu_*(mpc_params_.N-1);

    dynamics_A.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1) + 2*mpc_params_.N,(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)));
    SparseIdentity.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1) + 2*mpc_params_.N,(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)));
    dynamics_b_lb.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1) + 2*mpc_params_.N); // dynamics and torque bounds and vel_bound
    dynamics_b_ub.resize(nx_*mpc_params_.N + nu_*(mpc_params_.N-1) + 2*mpc_params_.N); // dynamics and torque bounds and vel_bound
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

    settings.verbose = PRINT_TIMING;
    settings.polish = true;
    settings.max_iter = 200;
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
    // vel_bounds
    for (int i = 0; i < mpc_params_.N; i++) {
      for (int j = 0; j < 2; j++) {
        dynamics_A.insert(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)+i*2+j,i*4+j+2) = 1;
      }
    }
    // set torque limits
    for (int i = 0; i < mpc_params_.N-1; i++) {
      for (int j = 0; j < nu_; j++) {
        dynamics_b_lb(nx_*mpc_params_.N+i*nu_+j) = -mpc_params_.tau_max;
        dynamics_b_ub(nx_*mpc_params_.N+i*nu_+j) = mpc_params_.tau_max;
      }
    }
    // set vel limits
    for (int i = 0; i < mpc_params_.N-1; i++) {
      for (int j = 0; j < 2; j++) {
        dynamics_b_lb(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)+i*2+j) = -mpc_params_.vel_max;
        dynamics_b_ub(nx_*mpc_params_.N+nu_*(mpc_params_.N-1)+i*2+j) = mpc_params_.vel_max;
      }
    }
}

void MPC::buildConstraintInequality(const std::vector<matrix_t> A_constraint, const std::vector<vector_t> b_constraint) {
    // const int num_constraints = A_constraint.size();
    // int total_rows = std::accumulate(A_constraint.begin(), A_constraint.end(), 0, 
                            //  [](int sum, const matrix_t& mat) { return sum + mat.rows(); });
    // General case ^ not neded here
    const int num_constraints = mpc_params_.N;
    int total_rows = mpc_params_.N;

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

Eigen::MatrixXd createPathLengthMatrix(int n) {
    Eigen::MatrixXd mat = Eigen::MatrixXd::Zero(n, n);

    mat(0, 0) = 1;
    mat(n - 1, n - 1) = 1;

    for (int i = 1; i < n - 1; ++i) {
        mat(i, i) = 2;
        mat(i, i - 1) = -1;
        mat(i, i + 1) = -1;
    }

    // Handle the first and last row separately
    mat(0, 1) = -1;
    mat(n - 1, n - 2) = -1;

    return mat;
}

void MPC::buildCost()
{
    H_cost.resize(mpc_params_.N*nx_ + (mpc_params_.N-1) * nu_, mpc_params_.N*nx_ + (mpc_params_.N-1) * nu_);
    H_cost.setZero();
    for (int i = 0; i < mpc_params_.N; i++) {
        for (int j = 0; j < nx_; j++) {
            if (i == mpc_params_.N-1) {
                    H.insert(i*nx_+j,i*nx_+j) = mpc_params_.terminalScaling*mpc_params_.stateScaling(j);
                    H_cost(i*nx_+j,i*nx_+j) = mpc_params_.terminalScaling*mpc_params_.stateScaling(j);
            } else {
                    H.insert(i*nx_+j,i*nx_+j) = mpc_params_.stateScaling(j);
                    H_cost(i*nx_+j,i*nx_+j) = mpc_params_.stateScaling(j);
            }
        }
        if (i < mpc_params_.N-1) {
            for (int j = 0; j < nu_; j++) {
                H.insert(mpc_params_.N*nx_+i*nu_+j,mpc_params_.N*nx_+i*nu_+j) = mpc_params_.inputScaling(j);
                H_cost(mpc_params_.N*nx_+i*nu_+j,mpc_params_.N*nx_+i*nu_+j) = mpc_params_.inputScaling(j);
            }
        }
    }

    matrix_t pathLengthCost = createPathLengthMatrix(mpc_params_.N);
    for (int i = 0; i < mpc_params_.N; i++) {
        for (int j = 0; j < nx_; j++) {
            H.coeffRef(i*nx_+j,i*nx_+j) += pathLengthCost(i,i) * mpc_params_.stateScaling(j) * mpc_params_.path_length_cost;
            if (i > 0) {
                H.insert(i*nx_+j,(i-1)*nx_+j) = pathLengthCost(i,(i-1)) * mpc_params_.stateScaling(j) * mpc_params_.path_length_cost;
            }
            if (i < mpc_params_.N - 1) {
                H.insert(i*nx_+j,(i+1)*nx_+j) = pathLengthCost(i,(i+1)) * mpc_params_.stateScaling(j) * mpc_params_.path_length_cost;
            }
        }
    }

    f.setZero();
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
            if (dynamics_A.coeff(i, j) != 0)
                constraints.insert(i, j) = dynamics_A.coeff(i, j);
        }
    }
    for (int i = 0; i < (mpc_params_.N-1); i++) {
        for (int j = 0; j < nx_; j++) {
            constraints.insert(4*i   + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i  , i*nx_ + j);
            constraints.insert(4*i   + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i  , (i+1)*nx_ + j);
            constraints.insert(4*i+1 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+1, i*nx_ + j);
            constraints.insert(4*i+1 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+1, (i+1)*nx_ + j);
            constraints.insert(4*i+2 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+2, i*nx_ + j);
            constraints.insert(4*i+2 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+2, (i+1)*nx_ + j);
            constraints.insert(4*i+3 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+3, i*nx_ + j);
            constraints.insert(4*i+3 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+3, (i+1)*nx_ + j);
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
    isInitialized = true;

}

void MPC::updateConstraints(const vector_t& x0)
{
    dynamics_b_lb.segment(0,nx_) = x0;
    dynamics_b_ub.segment(0,nx_) = x0;

    lb.segment(0,nx_) = x0;
    ub.segment(0,nx_) = x0;
    lb.segment(dynamics_b_lb.size(), constraint_b.size()) << constraint_b;

    for (int i = 0; i < (mpc_params_.N-1); i++) {
        for (int j = 0; j < nx_; j++) {
            constraints.coeffRef(4*i   + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i  , i*nx_ + j);
            constraints.coeffRef(4*i   + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i  , (i+1)*nx_ + j);
            constraints.coeffRef(4*i+1 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+1, i*nx_ + j);
            constraints.coeffRef(4*i+1 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+1, (i+1)*nx_ + j);
            constraints.coeffRef(4*i+2 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+2, i*nx_ + j);
            constraints.coeffRef(4*i+2 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+2, (i+1)*nx_ + j);
            constraints.coeffRef(4*i+3 + dynamics_A.rows(), i*nx_ + j)     = constraint_A(4*i+3, i*nx_ + j);
            constraints.coeffRef(4*i+3 + dynamics_A.rows(), (i+1)*nx_ + j) = constraint_A(4*i+3, (i+1)*nx_ + j);
        }
    }
}

void MPC::updateConstraintsSQP(ObstacleCollector O, vector_t sol, const vector_t& xg) {
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
        A.setZero();
        b << 1;
        double max_viol = 1e3;
        for (auto obstacle : O.obstacles) {
            vector_t x(2);
            x << sol.segment(nx_*i,2);
            vector_t A_hyp(2);
            double b_hyp;
            A_hyp.setZero();
            b_hyp = 0;
            getSeparatingHyperplane(obstacle, x, A_hyp, b_hyp);

            double v = (A_hyp.transpose() * x - b_hyp);
            if (v < max_viol) {
                A << A_hyp.transpose(),0,0;
                b << b_hyp;
                b += mpc_params_.buffer*vector_t::Ones(1);
                max_viol = v;
            }
        }
        A_constraint.push_back(A);
        b_constraint.push_back(b);
    }
    buildConstraintInequality(A_constraint, b_constraint);
    f = -H_cost*f_ref;
}

void MPC::updateCost()
{
    auto status = solver.SetObjectiveVector(f);
}

// Solve the MPC problem
vector_t MPC::solve(ObstacleCollector O, vector_t sol, const vector_t& x0, const vector_t& xg) {
    Timer timer(PRINT_TIMING);
    timer.start();
    vector_t mpc_sol = sol;
    for (int i = 0; i < mpc_params_.SQP_iters; i++) {
        updateConstraintsSQP(O, mpc_sol, xg);
        timer.time("        Update Constraints: ");
        updateCost();
        timer.time("        Update Cost: ");

        updateConstraints(x0);
        timer.time("        Update Constraints: ");
        auto status1 = solver.UpdateConstraintMatrix(constraints);
        timer.time("        Update UpdateConstraintMatrix: ");
        auto status2 = solver.SetBounds(lb, ub);
        timer.time("        Update SetBounds: ");
        
        OsqpExitCode exit_code = solver.Solve();
        timer.time("        Solve: ");
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
