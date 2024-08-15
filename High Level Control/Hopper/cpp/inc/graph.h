#pragma once 
#include "Types.h"

#include "osqp++.h"
#include <random>
#include <fstream>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/copy.hpp>

const double distance_tol = 0.5;
const double viol_tol = 0.25;
const double dim_of_R = .3;

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

    void setupQP(OsqpInstance& instance, const std::vector<matrix_t> vertices, const Obstacle obstacle);
    int initializeQP(OsqpSolver &solver, OsqpInstance instance, OsqpSettings settings);
    int solveQP(OsqpSolver &solver);

    void buildConstraintMatrix(Obstacle obstacle,
                            const int num_pts, const int num_obstacle_faces, const std::vector<matrix_t> vertices);
    void updateConstraints(OsqpSolver &solver, Obstacle obstacle,
                            const int num_pts, const int num_obstacle_faces, const std::vector<matrix_t> vertices);
};

const double kInfinity = std::numeric_limits<double>::infinity();

typedef adjacency_list<vecS, vecS, directedS, no_property, property<edge_weight_t, double>> Graph;
typedef graph_traits<Graph>::vertex_descriptor Vertex;

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

template <typename Derived>
void log(std::vector<Derived> data, std::ofstream& file, std::string name)
{
    file << name << " = [" << std::endl;
    for (auto d : data)
    {
        file << d.transpose() << std::endl;
    }
    file << "];" << std::endl;
}

std::vector<vector_4t> generateUniformPoints(int n, double min_x, double max_x, double min_y, double max_y,
                                                    double min_dx, double max_dx, double min_dy, double max_dy);
bool adjacent(vector_4t p1, vector_4t p2);
double get_weight(vector_4t p1, vector_4t p2);

std::ofstream open_log_file(std::string filename);
Graph buildGraph(std::vector<vector_4t> points);
void cutEdges(Graph &g, const int num_pts, const std::vector<matrix_t> vertices, const int num_obstacle_faces, VectorXd optimal_solution);
void solveGraph(std::vector<vector_4t> points, vector_4t starting_loc, vector_4t ending_loc,
    int &starting_ind, int& ending_ind, Graph g, std::vector<double> d, std::vector<Vertex>& p);
std::vector<matrix_t> getVerticesOfBezPoly(const std::vector<vector_4t> points);

template<typename Func>
Graph buildGraph(std::vector<vector_4t> points, Func&& F_G)
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

// does there exist and edge from p1 to p2.
template<typename Func>
bool adjacent(vector_4t p1, vector_4t p2, Func&& F_G) {
    bool connected = true; //innocent until proven guilty

    matrix_t xbar_x(2,1);
    matrix_t xbar_y(2,1);
    matrix_t x(4,1);
    matrix_t y(4,1);
    matrix_t f_xbar_x(1,1);
    matrix_t f_xbar_y(1,1);
    f_xbar_x.setZero();
    f_xbar_y.setZero(); // double integrator has no nonlinearities.
    matrix_t g_xbar_x(1,1);
    matrix_t g_xbar_y(1,1);
    g_xbar_x << 1;
    g_xbar_y << 1;
    matrix_t F, G;

    x << p1(0), p1(2), p2(0), p2(2);
    y << p1(1), p1(3), p2(1), p2(3);
    xbar_x << p1(0), p1(2);
    xbar_y << p1(1), p1(3);

    F_G(xbar_x, f_xbar_x, g_xbar_x, F, G);
    if (((F * x - G).array() > 0).any()) {
        connected = false;
    }
    F_G(xbar_y, f_xbar_y, g_xbar_y, F, G);
    if (((F * y - G).array() > 0).any()) {
        connected = false;
    }
    return connected;

}

