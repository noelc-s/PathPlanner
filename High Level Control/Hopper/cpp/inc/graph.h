#pragma once 
#include "Types.h"

#include "osqp++.h"
#include <random>
#include <fstream>

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

// does there exist and edge from p1 to p2.
bool adjacent(vector_4t p1, vector_4t p2, const matrix_t &Fx, const matrix_t & Gx,const matrix_t & Fy, const matrix_t & Gy, const matrix_t D_nT);

template<typename Func>
Graph buildGraph(std::vector<vector_4t> points, Func&& F_G, const matrix_t D_nT)
{
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
                add_edge(i, j, get_weight(points[i], points[j]), g);
            }
        }
    }
    return g;
}

