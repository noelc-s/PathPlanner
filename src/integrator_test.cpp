
#include "../inc/Types.h"
#include <chrono>
#include <iostream>

#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>

using namespace boost::numeric::odeint;

scalar_t sigma = 1;
scalar_t R = 1;
scalar_t b = 1;

void lorenz_eigen( const vector_t &x , vector_t &dxdt , scalar_t t )
{
    dxdt.resize(x.size());
    dxdt.setZero();
    for (int i = 0; i < x.size(); i++)
    {
        dxdt(i) = x(i);
    }
}


int main() {
     // simple example with eigen types
 vector_t x0_eigen(3000);
 x0_eigen.setIdentity();
 auto start = std::chrono::high_resolution_clock::now();
 typedef runge_kutta4<vector_t, scalar_t, vector_t, scalar_t, vector_space_algebra> stepper_lorenz;
 stepper_lorenz* step = new stepper_lorenz();
 for (int i = 0; i < 100; i++) {
   //integrate_const( stepper(), lorenz_eigen , x0_eigen , 0.0 , 0.1 , 0.1);
   step->do_step(lorenz_eigen, x0_eigen, 0.1, 0.1);
 }
 auto finish = std::chrono::high_resolution_clock::now();
 std::cout << "integration took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()*1e-3 << "micro s\n";
 std::cout << "per step integration took: " << std::chrono::duration_cast<std::chrono::nanoseconds>(finish-start).count()*1e-3/100 << "micro s\n";
}