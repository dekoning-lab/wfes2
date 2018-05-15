#pragma once

#include <cassert>
#include <cstdlib> 
#include <cmath>
#include <cfloat>
#include <cstring>

#include <iostream>
#include <limits>

// Fot arg parsing:
#include <deque>
#include <vector>
#include <sstream>

#include <mkl.h>

#include "MoveOnly.hpp"

#define DPF "%.10e" // Double print format
#define LPF "%lld" // Long long int print format
#define IPF "%d" // Int print format

typedef long long int llong;

template<typename T>
void print_buffer(T* buffer, size_t size, std::ostream& os = std::cout, bool newline = true) {
    for(size_t i = 0; i < size; i++) {
        os << buffer[i] << "\t";
    }
    if (newline) std::cout << std::endl;
}
