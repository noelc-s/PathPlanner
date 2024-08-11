#include "osqp++.h"
#include <Eigen/Dense>
#include <random>
#include <vector>
#include <fstream>
// #include "matplotlibcpp.h"

// namespace plt = matplotlibcpp;

using namespace osqp;
using ::Eigen::SparseMatrix;
using ::Eigen::Triplet;
using ::Eigen::VectorXd;

using vector_2t = Eigen::Matrix<double, 2, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

std::vector<vector_2t> generateUniformPoints(int n, double min_x, double max_x, double min_y, double max_y) {
    std::vector<vector_2t> points;
    points.reserve(n);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis_x(min_x, max_x);
    std::uniform_real_distribution<> dis_y(min_y, max_y);

    for (int i = 0; i < n; ++i) {
        vector_2t point;
        point << dis_x(gen), dis_y(gen);
        points.push_back(point);
    }

    return points;
}

int main()
{
    const int state_size = 2;
    const int num_adjacent_pts = 4;
    const int num_obstacle_faces = 4;
    const double kInfinity = std::numeric_limits<double>::infinity();

    matrix_t A_obstacle(num_obstacle_faces,state_size);
    vector_t b_obstacle(num_obstacle_faces);

    A_obstacle << 1, 0,
                 -1, 0,
                  0, 1,
                  0, -1;
    b_obstacle << 1.0, 1.0, 1.0, 1.0;

    const int num_pts = 5000;
    std::vector<vector_2t> points = generateUniformPoints(num_pts,-3,3,-3,3);
    std::vector<matrix_t> vertices;
    const double dim_of_R = .3;

    std::ofstream output_file("../output.m");
    // Check if the file is open
    if (!output_file.is_open()) {
        std::cerr << "Error: Could not open the file!" << std::endl;
        return 1;
    }
    output_file << "O = [" << std::endl;

    for (int i = 0; i < num_pts; i++)
    {
        matrix_t v;
        v.resize(state_size,num_adjacent_pts);

        v << points[i][0]+dim_of_R, points[i][0]+dim_of_R, points[i][0]-dim_of_R,points[i][0]-dim_of_R,
                    points[i][1]+dim_of_R, points[i][1]-dim_of_R, points[i][1]-dim_of_R,points[i][1]+dim_of_R;
        vertices.push_back(v);
        output_file << points[i].transpose() << std::endl << v.transpose() << std::endl << std::endl;
    }
    output_file << "];" << std::endl;

    // Setup the OSQP instance
    OsqpInstance instance;
    
    // decision variables are lambda_i and slack_i for each adj pt. [l1 l2 ... s1 s2 ...]
    const int num_dec_per_pt = num_adjacent_pts+num_obstacle_faces;
    const int num_const_per_pt = 1+num_adjacent_pts+num_obstacle_faces;
    SparseMatrix<double> objective_matrix(num_pts*num_dec_per_pt, num_pts*num_dec_per_pt);
    objective_matrix.setIdentity(); // cost does not matter right now
    instance.objective_matrix = 0.001*objective_matrix;
    instance.objective_vector.resize(num_pts*num_dec_per_pt);
    instance.objective_vector.setZero(); 

    // Constraint matrix (A)
    // simplex constraint -> sum to 1, elementwise positive, obstacle for each adj
    instance.lower_bounds.resize(num_pts*num_const_per_pt);
    instance.upper_bounds.resize(num_pts*num_const_per_pt);
    SparseMatrix<double> constraint_matrix(num_pts*num_const_per_pt, num_pts*num_dec_per_pt);
    std::vector<Triplet<double>> tripletsA;

    for (int p = 0; p < num_pts; p++)
    {
        for (int j = 0; j < num_adjacent_pts; ++j) {
            tripletsA.emplace_back(p*num_const_per_pt,p*num_dec_per_pt+j,1); // lambda_i sum to 1
        }
        instance.lower_bounds[p*num_const_per_pt] = 1; // lambda_i sum to 1
        instance.upper_bounds[p*num_const_per_pt] = 1; // lambda_i sum to 1
        for (int i = 0; i < num_adjacent_pts; ++i) {
            tripletsA.emplace_back(p*num_const_per_pt+1+i, p*num_dec_per_pt+i, 1); // elementwise_positive
        }
            for (int i = 0; i < num_adjacent_pts; ++i) {
            instance.lower_bounds[p*num_const_per_pt+1+i] = 0; // elementwise_positive
            instance.upper_bounds[p*num_const_per_pt+1+i] = 1; // elementwise_positive
        }

        for (int i = 0; i < num_adjacent_pts; ++i) {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << A_obstacle * vertices[p].block(0,i,state_size,1);
            for (int j = 0; j < num_obstacle_faces; j++) {
                tripletsA.emplace_back(p*num_const_per_pt+1+num_adjacent_pts+j, p*num_dec_per_pt+i, constraint(j));
            }
            tripletsA.emplace_back(p*num_const_per_pt+1+num_adjacent_pts+i, p*num_dec_per_pt+num_adjacent_pts+i, -1);
        }

        for (int j = 0; j < num_obstacle_faces; ++j) {
            instance.lower_bounds[p*num_const_per_pt+1+num_adjacent_pts+j] = -kInfinity;
            instance.upper_bounds[p*num_const_per_pt+1+num_adjacent_pts+j] = b_obstacle[j];
        }
    }

    constraint_matrix.setFromTriplets(tripletsA.begin(), tripletsA.end());
    instance.constraint_matrix = constraint_matrix;

    // Debug check
    // std::cout << constraint_matrix << std::endl;
    // std::cout << instance.upper_bounds << std::endl;

    // Create and solve the problem 
    OsqpSolver solver;
    OsqpSettings settings;
    settings.verbose = true;
    settings.polish = true;

    auto status = solver.Init(instance, settings);
    if (!status.ok()) {
        std::cout << status << std::endl;
        std::cerr << "Solver initialization failed!" << std::endl;
        return -1;
    }

    OsqpExitCode exit_code = solver.Solve();
    if (exit_code != OsqpExitCode::kOptimal) {
        std::cerr << "Solver did not find an optimal solution!" << std::endl;
        return -1;
    }

    double optimal_objective = solver.objective_value();
    VectorXd optimal_solution = solver.primal_solution();

    output_file << "V = [" << std::endl;

    // Print the optimal solution and objective value
    output_file << optimal_solution << std::endl;
    output_file << "];" << std::endl;

    return 0;
}
