#include <chrono>
#include <iostream>

class Timer
{
public: 
    std::chrono::high_resolution_clock::time_point start_time;
    std::chrono::high_resolution_clock::time_point end_time;
    std::chrono::duration<double, std::nano> duration;

    void start();
    void time();
    void time(std::string info);
};