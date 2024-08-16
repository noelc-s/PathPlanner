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

// does there exist and edge from p1 to p2.
bool adjacent(vector_4t p1, vector_4t p2, const matrix_t &Fx, const matrix_t & Gx,const matrix_t & Fy, const matrix_t & Gy, const matrix_t D_nT) {
    bool connected = true; //innocent until proven guilty

    matrix_t x(4,1);
    matrix_t y(4,1);
    x << p1(0), p1(2), p2(0), p2(2);
    y << p1(1), p1(3), p2(1), p2(3);

    if (((Fx * D_nT * x - Gx).array() > 0).any()) {
        connected = false;
    }
    
    if (((Fy * D_nT * y - Gy).array() > 0).any()) {
        connected = false;
    }
    return connected;
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

std::vector<matrix_t> getVerticesOfBezPoly(const std::vector<vector_4t> points)
{
    const int num_pts = points.size();
    const int num_reachable_pts = pow(2,4);

    std::vector<matrix_t> vertices;

    for (int i = 0; i < num_pts; i++)
    {
        matrix_t v;
        v.resize(points[0].size(), num_reachable_pts);
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

GraphQP::GraphQP()
{

}

void GraphQP::setupQP(OsqpInstance& instance, const std::vector<matrix_t> vertices, const Obstacle obstacle)
{
    // decision variables are lambda_i and slack_i for each adj pt. [l1 l2 ... s1 s2 ...]
    // decision variables per point = number of vertices in reachable set + number of obstacle faces;

    const int total_adjacent_pts = std::accumulate(vertices.begin(), vertices.end(), 0, 
                                 [](int sum, const matrix_t& mat) { return sum + mat.cols(); });    

    const int num_obstacle_faces = obstacle.b.size();
    const int num_pts = vertices.size();
    
    SparseMatrix<double> objective_matrix(total_adjacent_pts + num_pts * num_obstacle_faces, 
                                total_adjacent_pts + num_pts * num_obstacle_faces);
    objective_matrix.setIdentity(); // cost does not matter right now
    instance.objective_matrix = objective_matrix;
    instance.objective_vector.resize(total_adjacent_pts + num_pts * num_obstacle_faces);
    instance.objective_vector.setZero();

    // Constraint matrix (A)
    // simplex constraint -> sum to 1, elementwise positive, obstacle for each adj
    const int num_constaints = num_pts + num_pts * num_obstacle_faces + total_adjacent_pts;

    constraint_matrix.resize(num_constaints, total_adjacent_pts + num_pts * num_obstacle_faces);
    lb.resize(num_constaints);
    ub.resize(num_constaints);
    buildConstraintMatrix(obstacle, num_pts, num_obstacle_faces, vertices);

    instance.lower_bounds = lb;
    instance.upper_bounds = ub;
    instance.constraint_matrix = constraint_matrix;

    // Debug check
    // std::cout << constraint_matrix << std::endl;
    // std::cout << instance.upper_bounds << std::endl;
}

void GraphQP::buildConstraintMatrix(Obstacle obstacle,
                        const int num_pts, const int num_obstacle_faces, const std::vector<matrix_t> vertices)
{
    std::vector<Triplet<double>> tripletsA;

    int variable_index = 0;
    int constraint_index = 0;
    int num_reachable_pts = -1;

    for (int p = 0; p < num_pts; p++)
    {
        num_reachable_pts = vertices[p].cols();
        for (int j = 0; j < num_reachable_pts; ++j)
        {
            tripletsA.emplace_back(constraint_index, variable_index + j, 1); // lambda_i sum to 1
        }
        lb[constraint_index] = 1; // lambda_i sum to 1
        ub[constraint_index] = 1; // lambda_i sum to 1
        for (int i = 0; i < num_reachable_pts; ++i)
        {
            tripletsA.emplace_back(constraint_index + 1 + i, variable_index + i, 1); // elementwise_positive
        }
        for (int i = 0; i < num_reachable_pts; ++i)
        {
            lb[constraint_index + 1 + i] = 0; // elementwise_positive
            ub[constraint_index + 1 + i] = 1; // elementwise_positive
        }

        for (int i = 0; i < num_reachable_pts; ++i)
        {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << obstacle.A * vertices[p].block(0, i, vertices[p].rows(), 1);
            for (int j = 0; j < num_obstacle_faces; j++)
            {
                tripletsA.emplace_back(constraint_index + 1 + num_reachable_pts + j, variable_index + i, constraint(j));
            }
        }
        for (int i = 0; i < num_obstacle_faces; i++)
        {
            tripletsA.emplace_back(constraint_index + 1 + num_reachable_pts + i, variable_index + num_reachable_pts + i, -1);
        }

        for (int j = 0; j < num_obstacle_faces; ++j)
        {
            lb[constraint_index + 1 + num_reachable_pts + j] = -kInfinity;
            ub[constraint_index + 1 + num_reachable_pts + j] = obstacle.b[j];
        }
        constraint_index += 1 + num_reachable_pts + num_obstacle_faces;
        variable_index += num_reachable_pts + num_obstacle_faces;
    }

    constraint_matrix.setFromTriplets(tripletsA.begin(), tripletsA.end());
}

void GraphQP::updateConstraints(OsqpSolver &solver, Obstacle obstacle, 
                    const int num_pts, const int num_obstacle_faces, const std::vector<matrix_t> vertices)
{
    const int total_adjacent_pts = std::accumulate(vertices.begin(), vertices.end(), 0, 
                                 [](int sum, const matrix_t& mat) { return sum + mat.cols(); }); 
    const int num_constraints = num_pts + num_pts * num_obstacle_faces + total_adjacent_pts;

    int variable_index = 0;
    int constraint_index = 0;
    int num_reachable_pts = -1;

    for (int p = 0; p < num_pts; p++)
    {
        num_reachable_pts = vertices[p].cols();
        
        for (int i = 0; i < num_reachable_pts; ++i)
        {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << obstacle.A * vertices[p].block(0, i, vertices[p].rows(), 1);
            for (int j = 0; j < num_obstacle_faces; j++)
            {
                constraint_matrix.coeffRef(constraint_index + 1 + num_reachable_pts + j, variable_index + i) = constraint(j);
            }
        }

        for (int j = 0; j < num_obstacle_faces; ++j)
        {
            lb[constraint_index + 1 + num_reachable_pts + j] = -kInfinity;
            ub[constraint_index + 1 + num_reachable_pts + j] = obstacle.b[j];
        }
        constraint_index += 1 + num_reachable_pts + num_obstacle_faces;
        variable_index += num_reachable_pts + num_obstacle_faces;
    }

    auto status1 = solver.UpdateConstraintMatrix(constraint_matrix);
    auto status2 = solver.SetBounds(lb, ub);
}

int GraphQP::initializeQP(OsqpSolver &solver, OsqpInstance instance, OsqpSettings settings)
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

int GraphQP::solveQP(OsqpSolver &solver)
{
    OsqpExitCode exit_code = solver.Solve();
    if (exit_code != OsqpExitCode::kOptimal)
    {
        std::cerr << "Solver did not find an optimal solution!" << std::endl;
        return -1;
    }
    return 0;
}

void cutEdges(Graph &g, const int num_pts, const std::vector<matrix_t> vertices, const int num_obstacle_faces, VectorXd optimal_solution)
{
    int dec_var_ind = 0;
    for (int i = 0; i < num_pts; i++) {
        if (optimal_solution.segment(dec_var_ind + vertices[i].cols(),num_obstacle_faces).norm() < viol_tol)
        {
            clear_vertex(i, g);
        }
        dec_var_ind += vertices[i].cols() + num_obstacle_faces;
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