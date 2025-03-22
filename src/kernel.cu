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

        __device__ void getSeparatingHyperplane(double *obstacle_A, double *obstacle_b, double *obstacle_v, const double *x, double *A_hyp, double &b_hyp, double &dist, const int num_faces)
    {
        int closest_point = -1;
        double closest_dist = 1e3;
        double dist_to_point = 0;
        for (int j = 0; j < num_faces; j++) {
            dist_to_point = (x[0] - obstacle_v[j]) * (x[0] - obstacle_v[j]) +  (x[1] - obstacle_v[j + num_faces])*(x[1] - obstacle_v[j + num_faces]);
            if (dist_to_point < closest_dist) {
                closest_point = j;
                closest_dist = dist_to_point;
            }
        }
        dist = closest_dist;

        // Assuming the ordering that edge j connects vertex j to j+1
        int plane_1_ind = closest_point;
        int plane_2_ind = (closest_point + 1) % num_faces;

        int num_constraint_violated = 0;

        if (obstacle_A[plane_1_ind] * x[0] + obstacle_A[plane_1_ind + num_faces] * x[1] - obstacle_b[plane_1_ind] > -1e-2) {
            A_hyp[0] = obstacle_A[plane_1_ind + num_faces*0];
            A_hyp[1] = obstacle_A[plane_1_ind + num_faces*1];
            num_constraint_violated++;
        }
        if (obstacle_A[plane_2_ind] * x[0] + obstacle_A[plane_2_ind + num_faces] * x[1] - obstacle_b[plane_2_ind] > -1e-2) {
            A_hyp[0] = obstacle_A[plane_2_ind + num_faces*0];
            A_hyp[1] = obstacle_A[plane_2_ind + num_faces*1];
            num_constraint_violated++;
        }

        if (num_constraint_violated == 0) {
            // printf("Constraint Violated . . . . . ");
            A_hyp[0] = 0;
            A_hyp[1] = 0;
            b_hyp = 1;
        } else {
            if (num_constraint_violated > 1) {
                A_hyp[0] = x[0] - obstacle_v[closest_point];
                A_hyp[1] = x[1] - obstacle_v[closest_point + num_faces];
            }
            double norm = sqrt(A_hyp[0] * A_hyp[0] + A_hyp[1] * A_hyp[1]);
            A_hyp[0] /= norm;
            A_hyp[1] /= norm;
            b_hyp = A_hyp[0] * obstacle_v[closest_point] + A_hyp[1] * obstacle_v[closest_point + num_faces];   
        }
    }

    __global__ void obstacleMembershipHeuristic(double *obstacle_A, double *obstacle_b, double *obstacle_v, const double *edges, int *member, const int* num_constaint_vec, const int* constraint_vec_ind, const int* obstacle_type, const int num_edges, const int num_obstacles)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= num_edges * num_obstacles)
            return;

        const int edge_number = i % num_edges;
        const int obstacle_number = i / num_edges;

        int num_faces = num_constaint_vec[obstacle_number];
        int constaint_index = constraint_vec_ind[obstacle_number];

        const double *edge = &edges[edge_number * 16];
        double *A = &obstacle_A[constaint_index * 4];
        double *b = &obstacle_b[constaint_index];
        double *v = &obstacle_v[obstacle_number * num_faces * 2];

        bool in_obstacle = false;
        bool in_freespace = true;
        for (int k = 0; k < 4; k++) {
            bool result_and = true;
            for (int j = 0; j < num_faces; j++)
            {
                bool result = (A[j] * edge[4*k] + A[j + num_faces] * edge[4*k+1] + A[j + 2*num_faces] * edge[4*k+2] + A[j + 3*num_faces] * edge[4*k+3] - b[j]) <= 0;
                if (!result) {
                    result_and = false;
                }
            }
            // if (result[0] & result[1] & result[2] & result[3])
            if (obstacle_type[obstacle_number] == 0 && result_and)
            {
                in_obstacle = true;
                break;
            }
            if (obstacle_type[obstacle_number] == 1 && !result_and)
            {
                in_freespace = false;
                break;
            }
        }
        if (in_obstacle && obstacle_type[obstacle_number] == 0)
        {
            member[i] = 1;
        }
        else if (in_freespace && obstacle_type[obstacle_number] == 1)
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
                getSeparatingHyperplane(A, b, v, &edge[j * 4], A_hyp_, b_hyp_, dist, num_faces);
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

        int num_constaint_vec[num_obstacles];
        int constraint_vec_ind[num_obstacles];
        int obstacle_type[num_obstacles];
        num_constaint_vec[0] = obstacles[0].b.rows();
        constraint_vec_ind[0] = 0;
        obstacle_type[0] = obstacles[0].occType == FREE ? 1 : 0;
        for (int i = 1; i < num_obstacles; i++) {
            num_constaint_vec[i] = obstacles[i].b.rows();
            constraint_vec_ind[i] = constraint_vec_ind[i - 1] + num_constaint_vec[i - 1];
            obstacle_type[i] = obstacles[i].occType == FREE ? 1 : 0;
        }
        int num_total_constraints = constraint_vec_ind[num_obstacles - 1] +  num_constaint_vec[num_obstacles - 1];

        // Memory sizes
        size_t obstacle_A_size = num_total_constraints * 4 * sizeof(double);        // 4x4 matrix
        size_t obstacle_b_size = num_total_constraints * sizeof(double);         // 4x1 vector
        size_t obstacle_v_size = num_total_constraints * 2 * sizeof(double);         // 4x1 vector
        size_t edges_size = num_edges * 16 * sizeof(double); // num_edges x 4x4 matrix
        size_t num_constaints = num_obstacles * sizeof(int);

        // Allocate memory on the device
        double *d_obstacle_A;
        double *d_obstacle_b;
        double *d_obstacle_v;
        double *d_edges;
        int *d_member;
        int *d_num_constaint_vec;
        int *d_constraint_vec_ind;
        int* d_obstacle_type;

        cudaMalloc((void **)&d_obstacle_A, obstacle_A_size);
        cudaMalloc((void **)&d_obstacle_b, obstacle_b_size);
        cudaMalloc((void **)&d_obstacle_v, obstacle_v_size);
        cudaMalloc((void **)&d_edges, edges_size);
        cudaMalloc((void **)&d_member, member_size);
        cudaMalloc((void **)&d_num_constaint_vec, num_constaints);
        cudaMalloc((void **)&d_constraint_vec_ind, num_constaints);
        cudaMalloc((void **)&d_obstacle_type, num_constaints);

        // Prepare edge data
        double obstacle_A_flat[num_total_constraints * 4];
        double obstacle_b_flat[num_total_constraints];
        double obstacle_v_flat[num_total_constraints * 2];
        double edges_flat[num_edges * 16];

        for (int o = 0; o < num_obstacles; o++) {
            // copy obstacle A
            for (int col = 0; col < 4; col++)
            {
                for (int row = 0; row < num_constaint_vec[o]; row++)
                {
                    obstacle_A_flat[constraint_vec_ind[o] * 4 + row + col*num_constaint_vec[o]] = obstacles[o].A(row, col);
                }
            }
            for (int row = 0; row < num_constaint_vec[o]; row++) {
                obstacle_b_flat[constraint_vec_ind[o] + row] = obstacles[o].b(row);
            }
            for (int col = 0; col < 2; col++)
            {
                for (int row = 0; row < num_constaint_vec[o]; row++)
                {
                    obstacle_v_flat[constraint_vec_ind[o] * 2 + row + col*num_constaint_vec[o]] = obstacles[o].v(row, col);
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
        cudaMemcpy(d_obstacle_v, obstacle_v_flat, obstacle_v_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_edges, edges_flat, edges_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_num_constaint_vec, num_constaint_vec, num_constaints, cudaMemcpyHostToDevice);
        cudaMemcpy(d_constraint_vec_ind, constraint_vec_ind, num_constaints, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_type, obstacle_type, num_constaints, cudaMemcpyHostToDevice);

        // Launch the kernel
        // CAUTION: THIS CANNOT BE MORE THAN YOUR TENSOR CORE COUNT
        int blockSize = 128;
        int gridSize = (num_obstacles * num_edges + blockSize - 1) / blockSize;

        // cudaEventRecord(start);
        obstacleMembershipHeuristic<<<gridSize, blockSize>>>(d_obstacle_A, d_obstacle_b, d_obstacle_v, d_edges, d_member, d_num_constaint_vec, d_constraint_vec_ind, d_obstacle_type, num_edges, num_obstacles);
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
        cudaFree(d_obstacle_v);
        cudaFree(d_edges);
        cudaFree(d_member);
        cudaFree(d_num_constaint_vec);
        cudaFree(d_constraint_vec_ind);
        cudaFree(d_obstacle_type);
    }

    __global__ void getAllHyperplanes(double *obstacle_A, double *obstacle_b, double *obstacle_v, const double *sol, double *A1, double *A2, double *b, double *dist, const int N, const int num_obstacles, const int* num_constaint_vec, const int* constraint_vec_ind)
    {
        int i = blockIdx.x * blockDim.x + threadIdx.x;
        if (i >= N * num_obstacles)
            return;

        const int sol_number = i % N;
        const int obstacle_number = i / N;

        const double *sol_k = &sol[sol_number * 4];
        double *A_obst = &obstacle_A[constraint_vec_ind[obstacle_number] * 4];
        double *b_obst = &obstacle_b[constraint_vec_ind[obstacle_number]];
        double *v_obst = &obstacle_v[constraint_vec_ind[obstacle_number] * 2];
        
        double A_hyp[2] = {0.0};
        double b_hyp = 0.0;
        double the_dist = 1e3;

        getSeparatingHyperplane(A_obst, b_obst, v_obst, sol_k, A_hyp, b_hyp, the_dist, num_constaint_vec[obstacle_number]);
        A1[i] = A_hyp[0];
        A2[i] = A_hyp[1];
        b[i] = b_hyp;
        dist[i] = the_dist;

    }

    void MPC_GetActiveConstraints(std::vector<Obstacle> obstacles, const vector_t &sol, vector_t &A1, vector_t &A2, vector_t &b, vector_t &dist)
    {
        int num_obstacles = obstacles.size();
        int N = sol.size() / 4;

        int num_constaint_vec[num_obstacles];
        int constraint_vec_ind[num_obstacles];
        num_constaint_vec[0] = obstacles[0].b.rows();
        constraint_vec_ind[0] = 0;
        for (int i = 1; i < num_obstacles; i++) {
            num_constaint_vec[i] = obstacles[i].b.rows();
            constraint_vec_ind[i] = constraint_vec_ind[i - 1] + num_constaint_vec[i - 1];
        }
        int num_total_constraints = constraint_vec_ind[num_obstacles - 1] +  num_constaint_vec[num_obstacles - 1];

        // Memory sizes
        size_t obstacle_A_size = num_total_constraints * 4 * sizeof(double);        // 4x4 matrix
        size_t obstacle_b_size = num_total_constraints * sizeof(double);         // 4x1 vector
        size_t obstacle_v_size = num_total_constraints * 2 * sizeof(double);         // 4x2 matrix
        size_t sol_size = sol.size() * sizeof(double);
        size_t A1_size = A1.size() * sizeof(double); 
        size_t A2_size = A2.size() * sizeof(double); 
        size_t b_size = b.size() * sizeof(double); 
        size_t dist_size = dist.size() * sizeof(double); 
        size_t num_constaints = num_obstacles * sizeof(int);

        // Allocate memory on the device
        double *d_obstacle_A;
        double *d_obstacle_b;
        double *d_obstacle_v;
        double *d_sol;
        double *d_A1;
        double *d_A2;
        double *d_b;
        double *d_dist;
        int *d_num_constaint_vec;
        int *d_constraint_vec_ind;

        cudaMalloc((void **)&d_obstacle_A, obstacle_A_size);
        cudaMalloc((void **)&d_obstacle_b, obstacle_b_size);
        cudaMalloc((void **)&d_obstacle_v, obstacle_v_size);
        cudaMalloc((void **)&d_sol, sol_size);
        cudaMalloc((void **)&d_A1, A1_size);
        cudaMalloc((void **)&d_A2, A2_size);
        cudaMalloc((void **)&d_b, b_size);
        cudaMalloc((void **)&d_dist, dist_size);
        cudaMalloc((void **)&d_num_constaint_vec, num_constaints);
        cudaMalloc((void **)&d_constraint_vec_ind, num_constaints);

        // Prepare edge data
        double obstacle_A_flat[num_total_constraints * 4];
        double obstacle_b_flat[num_total_constraints];
        double obstacle_v_flat[num_total_constraints * 2];

        for (int o = 0; o < num_obstacles; o++) {
            // copy obstacle A
            for (int col = 0; col < 4; col++)
            {
                for (int row = 0; row < num_constaint_vec[o]; row++)
                {
                    obstacle_A_flat[constraint_vec_ind[o] * 4 + row + col*num_constaint_vec[o]] = obstacles[o].A(row, col);
                }
            }
            for (int row = 0; row < num_constaint_vec[o]; row++) {
                obstacle_b_flat[constraint_vec_ind[o] + row] = obstacles[o].b(row);
            }
            for (int col = 0; col < 2; col++)
            {
                for (int row = 0; row < num_constaint_vec[o]; row++)
                {
                    obstacle_v_flat[constraint_vec_ind[o] * 2 + row + col*num_constaint_vec[o]] = obstacles[o].v(row, col);
                }
            }
        }

        // cudaEvent_t start, stop;
        // cudaEventCreate(&start);
        // cudaEventCreate(&stop);

        // Copy data to device
        cudaMemcpy(d_obstacle_A, obstacle_A_flat, obstacle_A_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_b, obstacle_b_flat, obstacle_b_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_obstacle_v, obstacle_v_flat, obstacle_v_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_sol, sol.data(), sol_size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_num_constaint_vec, num_constaint_vec, num_constaints, cudaMemcpyHostToDevice);
        cudaMemcpy(d_constraint_vec_ind, constraint_vec_ind, num_constaints, cudaMemcpyHostToDevice);

        // Launch the kernel
        // CAUTION: THIS CANNOT BE MORE THAN YOUR TENSOR CORE COUNT
        int blockSize = 128;
        int gridSize = (num_obstacles * N + blockSize - 1) / blockSize;

        // cudaEventRecord(start);
        getAllHyperplanes<<<gridSize, blockSize>>>(d_obstacle_A, d_obstacle_b, d_obstacle_v, d_sol, d_A1, d_A2, d_b, d_dist, N, num_obstacles, d_num_constaint_vec, d_constraint_vec_ind);
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
        cudaFree(d_obstacle_v);
        cudaFree(d_sol);
        cudaFree(d_A1);
        cudaFree(d_A2);
        cudaFree(d_b);
        cudaFree(d_dist);
        cudaFree(d_num_constaint_vec);
        cudaFree(d_constraint_vec_ind);
    }

}
