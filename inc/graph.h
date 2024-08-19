#pragma once 
#include "Types.h"
#include "obstacle.h"

#include <numeric> 
#include "osqp++.h"
#include <fstream>

const double distance_tol = 0.5;
const double viol_tol = 0.3;

using namespace osqp;
using namespace boost;
using ::Eigen::SparseMatrix;
using ::Eigen::Triplet;
using ::Eigen::VectorXd;

class GraphQP
{
public:
    GraphQP();

    SparseMatrix<double> constraint_matrix;
    vector_t lb, ub;

    void setupQP(OsqpInstance& instance, const std::vector<matrix_t> edges, const Obs obstacle);
    int initializeQP(OsqpSolver &solver, OsqpInstance instance, OsqpSettings settings);
    int solveQP(OsqpSolver &solver);

    void buildConstraintMatrix(Obs obstacle, const std::vector<matrix_t> edges);
    void updateConstraints(OsqpSolver &solver, Obs obstacle, const std::vector<matrix_t> edges);
};

const double kInfinity = std::numeric_limits<double>::infinity();

template <class PredecessorMap>
class record_predecessorsC : public dijkstra_visitor<>
{
public:
record_predecessorsC(PredecessorMap p)
    : m_predecessor(p) { }

template <class Edge, class Graph>
void edge_relaxed(Edge e, Graph& g) {
    // set the parent of the target(e) to source(e)
    put(m_predecessor, target(e, g), source(e, g));
}
protected:
PredecessorMap m_predecessor;
};

template <class PredecessorMap>
record_predecessorsC<PredecessorMap>
make_predecessor_recorder(PredecessorMap p) {
    return record_predecessorsC<PredecessorMap>(p);
}

bool adjacent(vector_4t p1, vector_4t p2);
double get_weight(vector_4t p1, vector_4t p2);

std::ofstream open_log_file(std::string filename);
Graph buildGraph(std::vector<vector_4t> points);
void cutGraphEdges(Graph &g, const std::vector<matrix_t> edges, std::vector<std::pair<int,int>> vertexInds,
                const int num_obstacle_faces, VectorXd optimal_solution);
void solveGraph(std::vector<vector_4t> points, vector_4t starting_loc, vector_4t ending_loc,
    int &starting_ind, int& ending_ind, Graph g, std::vector<double> d, std::vector<Vertex>& p);
std::vector<matrix_t> getBezEdges(const Graph graph, std::vector<std::pair<int,int>> &vertexInds);

// does there exist and edge from p1 to p2.
bool adjacent(vector_4t p1, vector_4t p2, const matrix_t &Fx, const matrix_t & Gx,const matrix_t & Fy, const matrix_t & Gy, const matrix_t D_nT);

template<typename Func>
Graph buildGraph(std::vector<vector_4t> points, Func&& F_G, const matrix_t& D_nT, const matrix_t& Bez)
{
    EdgeProperties ep;
    matrix_t controlPoints(4,4);
    vector_t x1_x2(2*points[0].size());
    const int num_pts = points.size();
    Graph g(num_pts);
    for (int i = 0; i < num_pts; ++i)
    {
        matrix_t xbar_x(2,1);
        matrix_t xbar_y(2,1);
        xbar_x << points[i](0), points[i](2);
        xbar_y << points[i](1), points[i](3);
        matrix_t f_xbar_x(1,1);
        matrix_t f_xbar_y(1,1);
        f_xbar_x.setZero();
        f_xbar_y.setZero(); // double integrator has no nonlinearities.
        matrix_t g_xbar_x(1,1);
        matrix_t g_xbar_y(1,1);
        g_xbar_x << 1;
        g_xbar_y << 1;
        // TODO: This function in Bezier tubes does a sketchy resize so we have to re-
        // initialize here. FIX this for general use!
        matrix_t Fx, Fy, Gx, Gy;
        F_G(xbar_x, f_xbar_x, g_xbar_x, Fx, Gx);
        F_G(xbar_y, f_xbar_y, g_xbar_y, Fy, Gy);
        for (int j = 0; j < num_pts; ++j)
        {
            if (i != j && adjacent(points[i], points[j], Fx, Gx, Fy, Gy, D_nT))
            {
                ep.weight = get_weight(points[i], points[j]);
                x1_x2 << points[i], points[j];
                matrix_t mul = Bez*x1_x2;
                controlPoints << Eigen::Map<Eigen::MatrixXd>(mul.data(),4,4);
                ep.controlPoints = controlPoints;
                ep.source_vertex_ind = i;
                ep.target_vertex_ind = j;
                add_edge(i, j, ep, g);
            }
        }
    }
    return g;
}

