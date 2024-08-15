#include <iostream>
#include <random>

__global__ void addMatrices(float *P_x, float *P_y, float *C, float *A_obs, float *b_obs, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < N)
    {
        float A_i__x;
        float b_i;
        float in = 1;

        for (int i = 0; i < 4; i++)
        {
            A_i__x = A_obs[i * 2] * P_x[idx] + A_obs[i * 2 + 1] * P_y[idx];
            b_i = b_obs[i];
            if (A_i__x > b_i)
            {
                in = 0;
            }
        }

        C[idx] = in;
    }
}

cudaEvent_t start, stop;

int main()
{

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(-2, 2);

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    const int N = 300;
    const int size = N * sizeof(float);

    // Host matrices
    float h_x[N], h_y[N], h_C[N];
    float A_O[8] = {-0.2577, 0.9662,
                    0.9545, -0.2983,
                    -0.9545, 0.2983,
                    0.2577, -0.9662};
    float b_O[4] = {0.7343, 0.5608, 0.7517, 0.6828};

    // Initialize matrices A and B
    for (int i = 0; i < N; ++i)
    {
        h_x[i] = static_cast<float>(dis(gen));
        h_y[i] = static_cast<float>(dis(gen));
    }

    // Device matrices
    float *d_A, *d_B, *d_C, *A_obs, *b_obs;
    cudaMalloc((void **)&d_A, size);
    cudaMalloc((void **)&d_B, size);
    cudaMalloc((void **)&d_C, size);
    cudaMalloc((void **)&A_obs, 8 * sizeof(float));
    cudaMalloc((void **)&b_obs, 4 * sizeof(float));

    // Copy matrices A and B to the device
    cudaMemcpy(d_A, h_x, size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_B, h_y, size, cudaMemcpyHostToDevice);
    cudaMemcpy(A_obs, A_O, 8 * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(b_obs, b_O, 4 * sizeof(float), cudaMemcpyHostToDevice);

    // Define the number of threads per block and the number of blocks
    int threadsPerBlock = 16;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    // Launch the kernel
    cudaEventRecord(start);
    addMatrices<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, A_obs, b_obs, N);
    cudaEventRecord(stop);

    // Copy the result matrix C back to the host
    cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);
    cudaEventSynchronize(stop);
    float milliseconds = 0;
    cudaEventElapsedTime(&milliseconds, start, stop);

    // Print the result
    for (int i = 0; i < N; ++i) {
        std::cout << h_x[i] << ", " << h_y[i] << ", " << h_C[i] << std::endl;
    }
    std::cout << std::endl;
    std::cout << "This took: " << milliseconds << "milliseconds" << std::endl;
    // Free device memory
    cudaFree(d_A);
    cudaFree(d_B);
    cudaFree(d_C);

    return 0;
}
