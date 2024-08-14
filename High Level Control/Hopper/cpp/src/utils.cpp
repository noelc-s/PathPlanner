#include "../inc/utils.h"


void Timer::start() {
    start_time = std::chrono::high_resolution_clock::now();
}

void Timer::time() {
    end_time = std::chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    std::cout << duration.count()*1e-6 << " ms" << std::endl;
    start_time = end_time;
}

void Timer::time(std::string info) {
    std::cout << info;
    time();
}