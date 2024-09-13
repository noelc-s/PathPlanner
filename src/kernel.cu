#include <kernel.hpp>
#include <Eigen/Core>

#include <iostream>
#include <stdio.h>

#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

// CUDA Version
namespace Kernel
{
    __device__ void getSeparatingHyperplane(double *obstacle_A, double *obstacle_b, double *obstacle_Adj, double *obstacle_v, const double *x, double *A_hyp, double &b_hyp)
    {
        int closest_point = -1;
        double closest_dist = 1e3;
        double dist_to_point = 0.0;
        double temp[2] = {0.0};

        for (int j = 0; j < 4; j++) {
            dist_to_point = 0.0;
            for (int k = 0; k < 2; k++) {
                temp[k] = x[k] - obstacle_v[j + k*4];
                dist_to_point += temp[k] * temp[k];
            }
            if (dist_to_point < closest_dist) {
                closest_point = j;
                closest_dist = dist_to_point;
            }
        }

        const int num_faces = 4;
        for (int j = 0; j < num_faces; j++) {
            if ((obstacle_Adj[closest_point + j * 4] > 0) && (obstacle_A[j] * x[0] + obstacle_A[j + 4] * x[1] - obstacle_b[j] > -1e-2)) {
                for (int k = 0; k < 2; k++) {
                    A_hyp[k] += obstacle_A[j + 4*k];
                }
            }
        }

        double norm = sqrt(A_hyp[0] * A_hyp[0] + A_hyp[1] * A_hyp[1]);
        A_hyp[0] /= norm;
        A_hyp[1] /= norm;
        b_hyp = A_hyp[0] * obstacle_v[closest_point] + A_hyp[1] * obstacle_v[closest_point + 4];
    }

    __global__ void obstacleMembershipHeuristic(double *obstacle_A, double *obstacle_b, double *obstacle_Adj, double *obstacle_v, const double *edges, int *member, const int num_edges, const int num_obstacles)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= num_edges * num_obstacles)
            return;

        const int edge_number = i % num_edges;
        const int obstacle_number = i / num_edges;

        const double *edge = &edges[edge_number * 16];
        double *A = &obstacle_A[obstacle_number *  16];
        double *b = &obstacle_b[obstacle_number * 4];
        double *Adj = &obstacle_Adj[obstacle_number * 16];
        double *v = &obstacle_v[obstacle_number * 8];

        bool in_obstacle = false;
        for (int k = 0; k < 4; k++) {
            bool result[4] = {false};
            for (int j = 0; j < 4; j++)
            {
                result[j] = (A[j] * edge[4*k] + A[j + 4] * edge[4*k+1] + A[j + 8] * edge[4*k+2] + A[j + 12] * edge[4*k+3] - b[j]) <= 0;
            }
            if (result[0] & result[1] & result[2] & result[3])
            {
                in_obstacle = true;
                break;
            }
        }
        if (in_obstacle)
        {
            member[i] = 1;
        }
        else
        {
            // If not in obstacle, then we have to see if we have a separating hyperplane
            double A_hyp[2] = {0.0};
            double b_hyp = 0.0;

            for (int j = 0; j < 4; j++) {
                double A_hyp_[2] {0.0};
                double b_hyp_ = 0.0;
                getSeparatingHyperplane(A, b, Adj, v, &edge[j * 4], A_hyp_, b_hyp_);
                A_hyp[0] += A_hyp_[0];
                A_hyp[1] += A_hyp_[1];
                b_hyp += b_hyp_;
            }

            A_hyp[0] /= 4;
            A_hyp[1] /= 4;
            b_hyp /= 4;

            bool safe = true;
            for (int j = 0; j < 4; j++) {
                double result = A_hyp[0] * edge[j * 4] + A_hyp[1] * edge[j * 4 + 1] - b_hyp;
                if (result < 0) {
                    safe = false;
                    break;
                }
            }

            member[i] = safe ? 0 : 2;

        }
    }

    void GraphQP_ObstacleMembershipHeuristic(std::vector<Obstacle> obstacles, const std::vector<matrix_t> &edges, int_vector_t &member)
    {
        int num_edges = edges.size();
        int num_obstacles = obstacles.size();
        int member_size = num_obstacles * num_edges * sizeof(int);

        // Memory sizes
        size_t obstacle_A_size = num_obstacles * 16 * sizeof(double);        // 4x4 matrix
        size_t obstacle_b_size = num_obstacles * 4 * sizeof(double);         // 4x1 vector
        size_t obstacle_Adj_size = num_obstacles * 16 * sizeof(double);         // 4x1 vector
        size_t obstacle_v_size = num_obstacles * 8 * sizeof(double);         // 4x1 vector
        size_t edges_size = num_edges * 16 * sizeof(double); // num_edges x 4x4 matrix

        // Allocate memory on the device
        double *d_obstacle_A;
        double *d_obstacle_b;
        double *d_obstacle_Adj;
        double *d_obstacle_v;
        double *d_edges;
        int *d_member;

        cudaMalloc((void **)&d_obstacle_A, obstacle_A_size);
        cudaMalloc((void **)&d_obstacle_b, obstacle_b_size);
        cudaMalloc((void **)&d_obstacle_Adj, obstacle_Adj_size);
        cudaMalloc((void **)&d_obstacle_v, obstacle_v_size);
        cudaMalloc((void **)&d_edges, edges_size);
        cudaMalloc((void **)&d_member, member_size);

        // Prepare edge data
        double obstacle_A_flat[num_obstacles * 16];
        double obstacle_b_flat[num_obstacles * 4];
        double obstacle_Adj_flat[num_obstacles * 16];
        double obstacle_v_flat[num_obstacles * 8];
        double edges_flat[num_edges * 16];

        for (int o = 0; o < num_obstacles; o++) {
            // copy obstacle A
            for (int col = 0; col < 4; col++)
            {
                for (int row = 0; row < 4; row++)
                {
                    obstacle_A_flat[o * 16 + row + col*4] = obstacles[o].A(row, col);
                    obstacle_Adj_flat[o * 16 + row + col*4] = obstacles[o].Adjacency(row, col);
                }
                obstacle_b_flat[o * 4 + col] = obstacles[o].b(col);
            }
            for (int col = 0; col < 2; col++)
            {
                for (int row = 0; row < 4; row++)
                {
                    obstacle_v_flat[o * 8 + row + col*4] = obstacles[o].v(row, col);
                }
            }
        }
        // copy edges
        for (int i = 0; i < num_edges; i++)
        {
            Eigen::MatrixXd mat = edges[i];
            for (int col = 0; col < 4; col++)
            {
                for (int row = 0; row < 4; row++)
                {
                    edges_flat[i * 16 + row + col*4] = mat(row, col);
                }
            }
        }

        // cudaEvent_t start, stop;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);

        // Copy data to device
        cudaMemcpy(d_obstacle_A, obstacle_A_flat, obstacle_A_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_b, obstacle_b_flat, obstacle_b_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_Adj, obstacle_Adj_flat, obstacle_Adj_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_v, obstacle_v_flat, obstacle_v_size, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_obstacle_A, obstacles[0].A.data(), obstacle_A_size, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_obstacle_b, obstacles[0].b.data(), obstacle_b_size, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_obstacle_Adj, obstacles[0].Adjacency.data(), obstacle_Adj_size, cudaMemcpyHostToDevice);
        // cudaMemcpy(d_obstacle_v, obstacles[0].v.data(), obstacle_v_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_edges, edges_flat, edges_size, cudaMemcpyHostToDevice);

        // Launch the kernel
        int blockSize = 256;
        int gridSize = (num_obstacles * num_edges + blockSize - 1) / blockSize;

        // cudaEventRecord(start);
        obstacleMembershipHeuristic<<<gridSize, blockSize>>>(d_obstacle_A, d_obstacle_b, d_obstacle_Adj, d_obstacle_v, d_edges, d_member, num_edges, num_obstacles);
        // cudaEventRecord(stop);

        // Copy the result back to the host
        cudaMemcpy(member.data(), d_member, member_size, cudaMemcpyDeviceToHost);

        // cudaEventSynchronize(stop);
        // float milliseconds = 0;
        // cudaEventElapsedTime(&milliseconds, start, stop);
        // printf("That took: %f ms\n", milliseconds);

        // Free device memory
        cudaFree(d_obstacle_A);
        cudaFree(d_obstacle_b);
        cudaFree(d_obstacle_Adj);
        cudaFree(d_obstacle_v);
        cudaFree(d_edges);
        cudaFree(d_member);
    }
}

// void GraphQP::ObstacleMembershipHeuristic(Obstacle obstacle, const std::vector<matrix_t> edges, int_vector_t &member)
// {
//     // 0 if out, 1 if in, 2 if uncertain
//     // #pragma omp parallel for
//     for (int i = 0; i < edges.size(); i++) {
//         vector_t A_hyp_(2), A_hyp(2);
//         scalar_t b_hyp_, b_hyp;
//         matrix_t coll = (obstacle.A * edges[i]).colwise() - obstacle.b;
//         if (((coll.array() <= 0).colwise().all()).any()) {
//             member[i] = 1;
//         } else {
//             A_hyp.setZero();
//             b_hyp = 0;
//             A_hyp_.setZero();
//             b_hyp_ = 0;
//             for (int j = 0; j < edges[i].cols(); j++) {
//                 getSeparatingHyperplane(obstacle, edges[i].block(0,j,2, 1), A_hyp_, b_hyp_);
//                 A_hyp += A_hyp_;
//                 b_hyp += b_hyp_;
//             }
//             A_hyp /= edges[i].cols();
//             b_hyp /= edges[i].cols();

//             bool safe = (((A_hyp.transpose() * edges[i].block(0, 0, 2, edges[i].cols())).array() - b_hyp).array() >= 0).all();
//             if (safe == 1) {
//                 member[i] = 0;
//             } else {
//                 member[i] = 2;
//             }
//         }
//     }
// }

// void getSeparatingHyperplane(Obstacle obstacle, vector_t x, vector_t &A_hyp, scalar_t &b_hyp)
// {
//     int closest_point = -1;
//     scalar_t closest_dist = 1e3;
//     scalar_t dist_to_point;
//     Eigen::Array<bool, Eigen::Dynamic, 1> inds;

//     // (obstacle.v.block(0,0,obstacle.v.rows(),2).rowwise() - x.transpose()).rowwise().squaredNorm().minCoeff(&closest_point);

//     for (int j = 0; j < obstacle.v.rows(); j++) {
//         dist_to_point = (x - obstacle.v.block(j,0,1,2).transpose()).squaredNorm();
//         if (dist_to_point < closest_dist) {
//             closest_point = j;
//             closest_dist = dist_to_point;
//         }
//     }
//     vector_t faces = obstacle.Adjacency.block(closest_point,0,1,obstacle.Adjacency.cols()).transpose();
//     inds = (obstacle.A.block(0,0,obstacle.A.rows(),2) * x - obstacle.b).array() > -1e-2 && faces.array() > 0;
//     for (int j = 0; j < inds.size(); j++) {
//         if (inds(j) > 0) {
//             A_hyp += obstacle.A.block(j,0,1,2).transpose();
//         }
//     }
//     A_hyp = A_hyp / A_hyp.norm();
//     b_hyp = (A_hyp.transpose() * obstacle.v.block(closest_point,0,1,2).transpose()).value();
// }
