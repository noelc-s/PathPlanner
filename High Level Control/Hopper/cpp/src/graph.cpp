#include "../inc/graph.h"

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

bool adjacent(vector_4t p1, vector_4t p2)
{
    return (p1 - p2).norm() < distance_tol;
}

double get_weight(vector_4t p1, vector_4t p2)
{
    return (p1 - p2).norm();
}

std::ofstream open_log_file(std::string filename)
{
    std::ofstream output_file(filename);
    if (!output_file.is_open())
    {
        throw std::runtime_error( "Error: Could not open the log file!" );
    }
    return output_file;
}

Graph buildGraph(std::vector<vector_4t> points)
{
    const int num_pts = points.size();
    Graph g(num_pts);
    for (int i = 0; i < num_pts; ++i)
    {
        for (int j = 0; j < num_pts; ++j)
        {
            if (i != j && adjacent(points[i], points[j]))
            {
                add_edge(i, j, get_weight(points[i], points[j]), g);
            }
        }
    }
    return g;
}

std::vector<matrix_t> getReachableVertices(const std::vector<vector_4t> points)
{
    const int num_pts = points.size();
    const int num_adjacent_pts = pow(2,4);

    std::vector<matrix_t> vertices;
    const double dim_of_R = .03;

    for (int i = 0; i < num_pts; i++)
    {
        matrix_t v;
        v.resize(points[0].size(), num_adjacent_pts);
        v << points[i][0] + dim_of_R, points[i][0] + dim_of_R, points[i][0] + dim_of_R, points[i][0] + dim_of_R,
                points[i][0] + dim_of_R, points[i][0] + dim_of_R, points[i][0] + dim_of_R, points[i][0] + dim_of_R,
                points[i][0] - dim_of_R, points[i][0] - dim_of_R, points[i][0] - dim_of_R, points[i][0] - dim_of_R,
                points[i][0] - dim_of_R, points[i][0] - dim_of_R, points[i][0] - dim_of_R, points[i][0] - dim_of_R,
            points[i][1] + dim_of_R, points[i][1] + dim_of_R, points[i][1] + dim_of_R, points[i][1] + dim_of_R,
                points[i][1] - dim_of_R, points[i][1] - dim_of_R, points[i][1] - dim_of_R, points[i][1] - dim_of_R,
                points[i][1] - dim_of_R, points[i][1] - dim_of_R, points[i][1] - dim_of_R, points[i][1] - dim_of_R,
                points[i][1] + dim_of_R, points[i][1] + dim_of_R, points[i][1] + dim_of_R, points[i][1] + dim_of_R,
            points[i][2] + dim_of_R, points[i][2] + dim_of_R, points[i][2] - dim_of_R, points[i][2] - dim_of_R,
                points[i][2] + dim_of_R, points[i][2] + dim_of_R, points[i][2] - dim_of_R, points[i][2] - dim_of_R,
                points[i][2] + dim_of_R, points[i][2] + dim_of_R, points[i][2] - dim_of_R, points[i][2] - dim_of_R,
                points[i][2] + dim_of_R, points[i][2] + dim_of_R, points[i][2] - dim_of_R, points[i][2] - dim_of_R,
            points[i][3] + dim_of_R, points[i][3] - dim_of_R, points[i][3] - dim_of_R, points[i][3] + dim_of_R,
                points[i][3] + dim_of_R, points[i][3] - dim_of_R, points[i][3] - dim_of_R, points[i][3] + dim_of_R,
                points[i][3] + dim_of_R, points[i][3] - dim_of_R, points[i][3] - dim_of_R, points[i][3] + dim_of_R,
                points[i][3] + dim_of_R, points[i][3] - dim_of_R, points[i][3] - dim_of_R, points[i][3] + dim_of_R,
        vertices.push_back(v);
        
    }
    return vertices;
}



void setupQP(OsqpInstance& instance, const std::vector<matrix_t> vertices, const Obstacle obstacle)
{
    const int num_adjacent_pts = pow(2,4); // TODO: make this stored via size of Reachable sets
    // decision variables are lambda_i and slack_i for each adj pt. [l1 l2 ... s1 s2 ...]
    const int num_obstacle_faces = obstacle.b_obstacle.size();
    const int num_pts = vertices.size();
    const int num_dec_per_pt = num_adjacent_pts + num_obstacle_faces;
    const int num_const_per_pt = 1 + num_adjacent_pts + num_obstacle_faces;
    SparseMatrix<double> objective_matrix(num_pts * num_dec_per_pt, num_pts * num_dec_per_pt);
    objective_matrix.setIdentity(); // cost does not matter right now
    instance.objective_matrix = objective_matrix;
    instance.objective_vector.resize(num_pts * num_dec_per_pt);
    instance.objective_vector.setZero();

    // Constraint matrix (A)
    // simplex constraint -> sum to 1, elementwise positive, obstacle for each adj
    instance.lower_bounds.resize(num_pts * num_const_per_pt);
    instance.upper_bounds.resize(num_pts * num_const_per_pt);
    SparseMatrix<double> constraint_matrix(num_pts * num_const_per_pt, num_pts * num_dec_per_pt);
    std::vector<Triplet<double>> tripletsA;

    for (int p = 0; p < num_pts; p++)
    {
        for (int j = 0; j < num_adjacent_pts; ++j)
        {
            tripletsA.emplace_back(p * num_const_per_pt, p * num_dec_per_pt + j, 1); // lambda_i sum to 1
        }
        instance.lower_bounds[p * num_const_per_pt] = 1; // lambda_i sum to 1
        instance.upper_bounds[p * num_const_per_pt] = 1; // lambda_i sum to 1
        for (int i = 0; i < num_adjacent_pts; ++i)
        {
            tripletsA.emplace_back(p * num_const_per_pt + 1 + i, p * num_dec_per_pt + i, 1); // elementwise_positive
        }
        for (int i = 0; i < num_adjacent_pts; ++i)
        {
            instance.lower_bounds[p * num_const_per_pt + 1 + i] = 0; // elementwise_positive
            instance.upper_bounds[p * num_const_per_pt + 1 + i] = 1; // elementwise_positive
        }

        for (int i = 0; i < num_adjacent_pts; ++i)
        {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << obstacle.A_obstacle * vertices[p].block(0, i, vertices[p].rows(), 1);
            for (int j = 0; j < num_obstacle_faces; j++)
            {
                tripletsA.emplace_back(p * num_const_per_pt + 1 + num_adjacent_pts + j, p * num_dec_per_pt + i, constraint(j));
            }
        }
        for (int i = 0; i < num_obstacle_faces; i++)
        {
            tripletsA.emplace_back(p * num_const_per_pt + 1 + num_adjacent_pts + i, p * num_dec_per_pt + num_adjacent_pts + i, -1);
        }

        for (int j = 0; j < num_obstacle_faces; ++j)
        {
            instance.lower_bounds[p * num_const_per_pt + 1 + num_adjacent_pts + j] = -kInfinity;
            instance.upper_bounds[p * num_const_per_pt + 1 + num_adjacent_pts + j] = obstacle.b_obstacle[j];
        }
    }

    constraint_matrix.setFromTriplets(tripletsA.begin(), tripletsA.end());
    instance.constraint_matrix = constraint_matrix;

    // Debug check
    // std::cout << constraint_matrix << std::endl;
    // std::cout << instance.upper_bounds << std::endl;
}

int initializeQP(OsqpSolver &solver, OsqpInstance instance, OsqpSettings settings)
{
    auto status = solver.Init(instance, settings);
    if (!status.ok())
    {
        std::cout << status << std::endl;
        std::cerr << "Solver initialization failed!" << std::endl;
        return -1;
    }
    return 0;
}

int solveQP(OsqpSolver &solver)
{
    OsqpExitCode exit_code = solver.Solve();
    if (exit_code != OsqpExitCode::kOptimal)
    {
        std::cerr << "Solver did not find an optimal solution!" << std::endl;
        return -1;
    }
    return 0;
}

void cutEdges(Graph &g, const int num_pts, const int num_adjacent_pts, const int num_obstacle_faces, VectorXd optimal_solution)
{
    for (int i = 0; i < num_pts; i++) {
        if (optimal_solution.segment(i*(num_adjacent_pts+num_obstacle_faces)+num_adjacent_pts,num_obstacle_faces).norm() < viol_tol)
        {
            clear_vertex(i, g);
        }
    }
}

void solveGraph(std::vector<vector_4t> points, vector_4t starting_loc, vector_4t ending_loc,
                int &starting_ind, int& ending_ind, Graph g, std::vector<double> d, std::vector<Vertex>& p)
{
    double closest_starting_dist = 1e5;
    double closest_ending_dist = 1e5;
    const int num_pts = points.size();
    for (int i = 0; i < num_pts; i++) {
        double s_dist = (points[i] - starting_loc).norm();
        if (s_dist < closest_starting_dist) {
            closest_starting_dist = s_dist;
            starting_ind = i;
        }
        double e_dist = (points[i] - ending_loc).norm();
        if (e_dist < closest_ending_dist) {
            closest_ending_dist = e_dist;
            ending_ind = i;
        }
    }

    dijkstra_shortest_paths(g, vertex(starting_ind, g),
                distance_map(&d[0]).visitor(make_predecessor_recorder(&p[0])));
}