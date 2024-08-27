#include "../inc/graph.h"

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
    EdgeProperties ep;
    vector_t x1_x2(2*points[0].size());
    const int num_pts = points.size();
    Graph g(num_pts);
    for (int i = 0; i < num_pts; ++i)
    {
        for (int j = 0; j < num_pts; ++j)
        {
            if (i != j && adjacent(points[i], points[j]))
            {
                ep.weight = get_weight(points[i], points[j]);
                x1_x2 << points[i], points[j];
                ep.controlPoints = x1_x2;
                ep.source_vertex_ind = i;
                ep.target_vertex_ind = j;
                add_edge(i, j, ep, g);
            }
        }
    }
    return g;
}

std::vector<matrix_t> getBezEdges(const Graph graph, std::vector<std::pair<int,int>> &vertexInds)
{
    std::vector<matrix_t> edgePolytopeVertices;
    std::pair<EdgeIterator, EdgeIterator> edge_pair = edges(graph);
    for (EdgeIterator ei = edge_pair.first; ei != edge_pair.second; ++ei) {
        edgePolytopeVertices.push_back(graph[*ei].controlPoints);
        vertexInds.push_back(std::pair(graph[*ei].source_vertex_ind, graph[*ei].target_vertex_ind));
    }
    return edgePolytopeVertices;
}

GraphQP::GraphQP()
{

}

void GraphQP::setupQP(OsqpInstance& instance, const std::vector<matrix_t> edges, const Obstacle obstacle)
{
    // decision variables are lambda_i and slack_i for each adj pt. [l1 l2 ... s1 s2 ...]
    // decision variables per point = number of vertices in edge control point polytope + number of obstacle faces;

    const int total_adjacent_pts = std::accumulate(edges.begin(), edges.end(), 0, 
                                 [](int sum, const matrix_t& mat) { return sum + mat.cols(); });    

    const int num_obstacle_faces = obstacle.b.size();
    const int num_edges = edges.size();
    
    SparseMatrix<double> objective_matrix(total_adjacent_pts + num_edges * num_obstacle_faces, 
                                total_adjacent_pts + num_edges * num_obstacle_faces);
    objective_matrix.setIdentity(); // cost does not matter right now
    instance.objective_matrix = objective_matrix;
    instance.objective_vector.resize(total_adjacent_pts + num_edges * num_obstacle_faces);
    instance.objective_vector.setZero();

    // Constraint matrix (A)
    // simplex constraint -> sum to 1, elementwise positive, obstacle for each adj
    const int num_constaints = num_edges + num_edges * num_obstacle_faces + total_adjacent_pts;

    constraint_matrix.resize(num_constaints, total_adjacent_pts + num_edges * num_obstacle_faces);
    lb.resize(num_constaints);
    ub.resize(num_constaints);
    buildConstraintMatrix(obstacle, edges);

    instance.lower_bounds = lb;
    instance.upper_bounds = ub;
    instance.constraint_matrix = constraint_matrix;

    // Debug check
    // std::cout << constraint_matrix << std::endl;
    // std::cout << instance.upper_bounds << std::endl;
}

void GraphQP::buildConstraintMatrix(Obstacle obstacle, const std::vector<matrix_t> edges)
{
    const int num_obstacle_faces = obstacle.b.size();
    std::vector<Triplet<double>> tripletsA;

    int variable_index = 0;
    int constraint_index = 0;
    int num_reachable_pts = -1;

    const int num_edges = edges.size();

    for (int p = 0; p < num_edges; p++)
    {
        num_reachable_pts = edges[p].cols();
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
            constraint << obstacle.A * edges[p].block(0, i, edges[p].rows(), 1);
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

void GraphQP::ObstacleMembershipHeuristic(Obstacle obstacle, const std::vector<matrix_t> edges, int_vector_t &member)
{
    // 0 if out, 1 if in, 2 if uncertain
    #pragma omp parallel for
    for (int i = 0; i < edges.size(); i++) {
        vector_t A_hyp_(2), A_hyp(2);
        double b_hyp_, b_hyp;
        matrix_t coll = (obstacle.A * edges[i]).colwise() - obstacle.b;
        if (((coll.array() <= 0).colwise().all()).any()) {
            member[i] = 1;
        } else {
            A_hyp.setZero();
            b_hyp = 0;
            A_hyp_.setZero();
            b_hyp_ = 0;
            for (int j = 0; j < edges[i].cols(); j++) {
                getSeparatingHyperplane(obstacle, edges[i].block(0,j,2, 1), A_hyp_, b_hyp_);
                A_hyp += A_hyp_;
                b_hyp += b_hyp_;
            }
            A_hyp /= edges[i].cols();
            b_hyp /= edges[i].cols();
            
            bool safe = (((A_hyp.transpose() * edges[i].block(0, 0, 2, edges[i].cols())).array() - b_hyp).array() >= 0).all();
            if (safe == 1) {
                member[i] = 0;
            } else {
                member[i] = 2;
            }
        }
    }
}

void GraphQP::updateConstraints(OsqpSolver &solver, Obstacle obstacle, const std::vector<matrix_t> edges)
{
    const int total_adjacent_pts = std::accumulate(edges.begin(), edges.end(), 0, 
                                 [](int sum, const matrix_t& mat) { return sum + mat.cols(); }); 
    const int num_edges = edges.size();
    const int num_obstacle_faces = obstacle.b.size();
    const int num_constraints = num_edges + num_edges * num_obstacle_faces + total_adjacent_pts;

    int variable_index = 0;
    int constraint_index = 0;
    int num_reachable_pts = -1;

    Timer timer(PRINT_TIMING);
    timer.start();
    for (int p = 0; p < num_edges; p++)
    {
        num_reachable_pts = edges[p].cols();
        
        for (int i = 0; i < num_reachable_pts; ++i)
        {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << obstacle.A * edges[p].block(0, i, edges[p].rows(), 1);
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
    timer.time("        Modification:");

    auto status1 = solver.UpdateConstraintMatrix(constraint_matrix);
    timer.time("        update constraint matrix:");
    auto status2 = solver.SetBounds(lb, ub);
    timer.time("        set bounds:");
}

void GraphQP::updateConstraints(OsqpSolver &solver, Obstacle obstacle, const std::vector<matrix_t> edges, const int_vector_t &member) {
    // constraint_matrix.setZero();
    const int num_edges = edges.size();
    const int num_obstacle_faces = obstacle.b.size();

    int variable_index = 0;
    int constraint_index = 0;
    int num_reachable_pts = -1;

    Timer timer(PRINT_TIMING);
    timer.start();
    for (int p = 0; p < num_edges; p++)
    {
        num_reachable_pts = edges[p].cols();
    
        for (int i = 0; i < num_reachable_pts; ++i)
        {
            matrix_t constraint(num_obstacle_faces, 1);
            constraint << obstacle.A * edges[p].block(0, i, edges[p].rows(), 1);
            for (int j = 0; j < num_obstacle_faces; j++)
            {
                if (member(p) == 2) {
                    constraint_matrix.coeffRef(constraint_index + 1 + num_reachable_pts + j, variable_index + i) = constraint(j);
                } else {
                    constraint_matrix.coeffRef(constraint_index + 1 + num_reachable_pts + j, variable_index + i) = 0;
                }
            }
        }

        for (int j = 0; j < num_obstacle_faces; ++j)
        {
            if (member(p) == 2) {
                lb[constraint_index + 1 + num_reachable_pts + j] = -kInfinity;
                ub[constraint_index + 1 + num_reachable_pts + j] = obstacle.b[j];
            } else {
                lb[constraint_index + 1 + num_reachable_pts + j] = -kInfinity;
                ub[constraint_index + 1 + num_reachable_pts + j] = 1;
            }
        }
        constraint_index += 1 + num_reachable_pts + num_obstacle_faces;
        variable_index += num_reachable_pts + num_obstacle_faces;
    }
    timer.time("        Modification:");

    auto status1 = solver.UpdateConstraintMatrix(constraint_matrix);
    timer.time("        update constraint matrix:");
    auto status2 = solver.SetBounds(lb, ub);
    timer.time("        set bounds:");
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

void cutGraphEdges(Graph &g, const std::vector<matrix_t> edges, std::vector<std::pair<int,int>> vertexInds, 
                    const int num_obstacle_faces, VectorXd optimal_solution)
{
    int dec_var_ind = 0;
    for (int i = 0; i < edges.size(); i++) {
        if (optimal_solution.segment(dec_var_ind + edges[i].cols(),num_obstacle_faces).norm() < viol_tol)
        {
            remove_edge(vertexInds[i].first, vertexInds[i].second, g);
        }
        dec_var_ind += edges[i].cols() + num_obstacle_faces;
    }
}

void cutGraphEdges(Graph &g, const std::vector<matrix_t> edges, std::vector<std::pair<int,int>> vertexInds, 
                    const int num_obstacle_faces, VectorXd optimal_solution, int_vector_t membership)
{
    int dec_var_ind = 0;
    for (int i = 0; i < edges.size(); i++) {
        if (membership[i] == 1) {
            remove_edge(vertexInds[i].first, vertexInds[i].second, g);
        } else if (membership[i] == 2) {
            if (optimal_solution.segment(dec_var_ind + edges[i].cols(),num_obstacle_faces).norm() < viol_tol)
            {
                remove_edge(vertexInds[i].first, vertexInds[i].second, g);
            }
            dec_var_ind += edges[i].cols() + num_obstacle_faces;
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
        int num_edges = out_degree(i, g);
        if (num_edges > 0 && s_dist < closest_starting_dist) {
            closest_starting_dist = s_dist;
            starting_ind = i;
        }
        double e_dist = (points[i] - ending_loc).norm();
        if (e_dist < closest_ending_dist) {
            closest_ending_dist = e_dist;
            ending_ind = i;
        }
    }

    // Create a property map that extracts the weight from the custom EdgeProperties
    auto weight_map = boost::weight_map(boost::get(&EdgeProperties::weight, g));

    dijkstra_shortest_paths(g, vertex(starting_ind, g),
                weight_map.distance_map(&d[0]).visitor(make_predecessor_recorder(&p[0])));
}
