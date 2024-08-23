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
using int_vector_t = Eigen::Matrix<int, Eigen::Dynamic, 1>;

// Define a property that includes a matrix_t
struct EdgeProperties {
    double weight;   // standard weight
    matrix_t controlPoints; // custom matrix property
    int source_vertex_ind;
    int target_vertex_ind;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::no_property, EdgeProperties> Graph;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::edge_iterator EdgeIterator;
