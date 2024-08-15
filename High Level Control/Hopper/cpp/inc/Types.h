# pragma once
#include <Eigen/Dense>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/copy.hpp>

using vector_2t = Eigen::Matrix<double, 2, 1>;
using vector_4t = Eigen::Matrix<double, 4, 1>;
using matrix_t = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using vector_t = Eigen::Matrix<double, Eigen::Dynamic, 1>;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::property<boost::edge_weight_t, double>> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;

struct Obstacle
{
    matrix_t A;
    vector_t b;
    matrix_t v; // vertex representation
    matrix_t Adjacency; // points to faces (1 at (i,j) if point i touches face j)
};