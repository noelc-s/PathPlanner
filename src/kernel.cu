#include <kernel.hpp>
#include <Eigen/Core>

#include <iostream>
#include <stdio.h>

#define HANDLE_ERROR(err) (HandleError(err, __FILE__, __LINE__))

// CUDA Version
namespace Kernel
{
    // __device__ void getSeparatingHyperplane(double *obstacle_A, double *obstacle_b, double *obstacle_Adj, double *obstacle_v, const double *x, double *A_hyp, double &b_hyp, double &dist)
    // {
    //     int closest_point = -1;
    //     double closest_dist = 1e3;
    //     double dist_to_point = 0.0;
    //     double temp[2] = {0.0};

    //     for (int j = 0; j < 4; j++) {
    //         dist_to_point = 0.0;
    //         for (int k = 0; k < 2; k++) {
    //             temp[k] = x[k] - obstacle_v[j + k*4];
    //             dist_to_point += temp[k] * temp[k];
    //         }
    //         if (dist_to_point < closest_dist) {
    //             closest_point = j;
    //             closest_dist = dist_to_point;
    //         }
    //     }
    //     dist = closest_dist;

    //     const int num_faces = 4;
    //     for (int j = 0; j < num_faces; j++) {
    //         if ((obstacle_Adj[closest_point + j * 4] > 0) && (obstacle_A[j] * x[0] + obstacle_A[j + 4] * x[1] - obstacle_b[j] > -1e-2)) {
    //             for (int k = 0; k < 2; k++) {
    //                 A_hyp[k] += obstacle_A[j + 4*k];
    //             }
    //         }
    //     }

    //     double norm = sqrt(A_hyp[0] * A_hyp[0] + A_hyp[1] * A_hyp[1]);
    //     A_hyp[0] /= norm;
    //     A_hyp[1] /= norm;
    //     b_hyp = A_hyp[0] * obstacle_v[closest_point] + A_hyp[1] * obstacle_v[closest_point + 4];
    // }

        __device__ void getSeparatingHyperplane(double *obstacle_A, double *obstacle_b, double *obstacle_Adj, double *obstacle_v, const double *x, double *A_hyp, double &b_hyp, double &dist)
    {
        int closest_point = -1;
        double closest_dist = 1e3;
        double dist_to_point = 0;
        for (int j = 0; j < 4; j++) {
            dist_to_point = (x[0] - obstacle_v[j]) * (x[0] - obstacle_v[j]) +  (x[1] - obstacle_v[j + 4])*(x[1] - obstacle_v[j + 4]);
            if (dist_to_point < closest_dist) {
                closest_point = j;
                closest_dist = dist_to_point;
            }
        }
        dist = closest_dist;

        int num_constraint_violated = 0;
        const int num_faces = 4;
        for (int j = 0; j < num_faces; j++) {
            if ((obstacle_Adj[closest_point + j * 4] > 0) && (obstacle_A[j] * x[0] + obstacle_A[j + 4] * x[1] - obstacle_b[j] > -1e-2)) {
                A_hyp[0] = obstacle_A[j + 4*0];
                A_hyp[1] = obstacle_A[j + 4*1];
                num_constraint_violated++;
            }
        }
        if (num_constraint_violated == 0) {
            // printf("Constraint Violated . . . . . ");
            A_hyp[0] = 0;
            A_hyp[1] = 0;
            b_hyp = 1;
        } else {
            if (num_constraint_violated > 1) {
                A_hyp[0] = x[0] - obstacle_v[closest_point];
                A_hyp[1] = x[1] - obstacle_v[closest_point + 4];
            }
            double norm = sqrt(A_hyp[0] * A_hyp[0] + A_hyp[1] * A_hyp[1]);
            A_hyp[0] /= norm;
            A_hyp[1] /= norm;
            b_hyp = A_hyp[0] * obstacle_v[closest_point] + A_hyp[1] * obstacle_v[closest_point + 4];   
        }
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
                double dist;
                getSeparatingHyperplane(A, b, Adj, v, &edge[j * 4], A_hyp_, b_hyp_, dist);
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
        cudaMemcpy(d_edges, edges_flat, edges_size, cudaMemcpyHostToDevice);

        // Launch the kernel
        // CAUTION: THIS CANNOT BE MORE THAN YOUR TENSOR CORE COUNT
        int blockSize = 128;
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

    __global__ void getAllHyperplanes(double *obstacle_A, double *obstacle_b, double *obstacle_Adj, double *obstacle_v, const double *sol, double *A1, double *A2, double *b, double *dist, const int N, const int num_obstacles)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N * num_obstacles)
            return;

        const int sol_number = i % N;
        const int obstacle_number = i / N;

        const double *sol_k = &sol[sol_number * 4];
        double *A_obst = &obstacle_A[obstacle_number *  16];
        double *b_obst = &obstacle_b[obstacle_number * 4];
        double *Adj_obst = &obstacle_Adj[obstacle_number * 16];
        double *v_obst = &obstacle_v[obstacle_number * 8];
        
        double A_hyp[2] = {0.0};
        double b_hyp = 0.0;
        double the_dist = 1e3;

        getSeparatingHyperplane(A_obst, b_obst, Adj_obst, v_obst, sol_k, A_hyp, b_hyp, the_dist);
        A1[i] = A_hyp[0];
        A2[i] = A_hyp[1];
        b[i] = b_hyp;
        dist[i] = the_dist;

    }

    void MPC_GetActiveConstraints(std::vector<Obstacle> obstacles, const vector_t &sol, vector_t &A1, vector_t &A2, vector_t &b, vector_t &dist)
    {
        int num_obstacles = obstacles.size();
        int N = sol.size() / 4;

        // Memory sizes
        size_t obstacle_A_size = num_obstacles * 16 * sizeof(double);        // 4x4 matrix
        size_t obstacle_b_size = num_obstacles * 4 * sizeof(double);         // 4x1 vector
        size_t obstacle_Adj_size = num_obstacles * 16 * sizeof(double);         // 4x4 matrix
        size_t obstacle_v_size = num_obstacles * 8 * sizeof(double);         // 4x2 matrix
        size_t sol_size = sol.size() * sizeof(double);
        size_t A1_size = A1.size() * sizeof(double); 
        size_t A2_size = A2.size() * sizeof(double); 
        size_t b_size = b.size() * sizeof(double); 
        size_t dist_size = dist.size() * sizeof(double); 

        // Allocate memory on the device
        double *d_obstacle_A;
        double *d_obstacle_b;
        double *d_obstacle_Adj;
        double *d_obstacle_v;
        double *d_sol;
        double *d_A1;
        double *d_A2;
        double *d_b;
        double *d_dist;

        cudaMalloc((void **)&d_obstacle_A, obstacle_A_size);
        cudaMalloc((void **)&d_obstacle_b, obstacle_b_size);
        cudaMalloc((void **)&d_obstacle_Adj, obstacle_Adj_size);
        cudaMalloc((void **)&d_obstacle_v, obstacle_v_size);
        cudaMalloc((void **)&d_sol, sol_size);
        cudaMalloc((void **)&d_A1, A1_size);
        cudaMalloc((void **)&d_A2, A2_size);
        cudaMalloc((void **)&d_b, b_size);
        cudaMalloc((void **)&d_dist, dist_size);

        // Prepare edge data
        double obstacle_A_flat[num_obstacles * 16];
        double obstacle_b_flat[num_obstacles * 4];
        double obstacle_Adj_flat[num_obstacles * 16];
        double obstacle_v_flat[num_obstacles * 8];

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

        // cudaEvent_t start, stop;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);

        // Copy data to device
        cudaMemcpy(d_obstacle_A, obstacle_A_flat, obstacle_A_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_b, obstacle_b_flat, obstacle_b_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_Adj, obstacle_Adj_flat, obstacle_Adj_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_v, obstacle_v_flat, obstacle_v_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_sol, sol.data(), sol_size, cudaMemcpyHostToDevice);

        // Launch the kernel
        // CAUTION: THIS CANNOT BE MORE THAN YOUR TENSOR CORE COUNT
        int blockSize = 128;
        int gridSize = (num_obstacles * N + blockSize - 1) / blockSize;

        // cudaEventRecord(start);
        getAllHyperplanes<<<gridSize, blockSize>>>(d_obstacle_A, d_obstacle_b, d_obstacle_Adj, d_obstacle_v, d_sol, d_A1, d_A2, d_b, d_dist, N, num_obstacles);
        // cudaEventRecord(stop);

        // Copy the result back to the host
        cudaMemcpy(A1.data(), d_A1, A1_size, cudaMemcpyDeviceToHost);
        cudaMemcpy(A2.data(), d_A2, A2_size, cudaMemcpyDeviceToHost);
        cudaMemcpy(b.data(), d_b, b_size, cudaMemcpyDeviceToHost);
        cudaMemcpy(dist.data(), d_dist, dist_size, cudaMemcpyDeviceToHost);

        // cudaEventSynchronize(stop);
        // float milliseconds = 0;
        // cudaEventElapsedTime(&milliseconds, start, stop);
        // printf("That took: %f ms\n", milliseconds);

        // Free device memory
        cudaFree(d_obstacle_A);
        cudaFree(d_obstacle_b);
        cudaFree(d_obstacle_Adj);
        cudaFree(d_obstacle_v);
        cudaFree(d_sol);
        cudaFree(d_A1);
        cudaFree(d_A2);
        cudaFree(d_b);
        cudaFree(d_dist);
    }

}
