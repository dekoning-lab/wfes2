#pragma once

#include <cassert>
#include <cstdlib> 
#include <cmath>
#include <cfloat>
#include <cstring>

#include <iostream>
#include <fstream>
#include <limits>

#include <deque>
#include <vector>
#include <sstream>
#include <chrono>
#include <stdexcept>

#include <mkl.h>

#include <Eigen/Core>

#define DPF "%.10e" // Double print format
#define LPF "%lld" // Long long int print format
#define IPF "%d" // Int print format

typedef long long int llong;

typedef std::chrono::time_point<std::chrono::system_clock> time_point;
typedef std::chrono::duration<double> time_diff;

template<typename T>
void print_buffer(T* buffer, size_t size, std::ostream& os = std::cout, bool newline = true) {
    for(size_t i = 0; i < size; i++) {
        os << buffer[i] << "\t";
    }
    if (newline) std::cout << std::endl;
}

typedef Eigen::VectorXd dvec;
typedef Eigen::Matrix<llong, Eigen::Dynamic, 1> lvec;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> dmat;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> diagmat;
