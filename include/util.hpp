#pragma once

#include "common.hpp"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

void write_matrix_to_file(const dmat& A, std::string path);
void write_vector_to_file(const dvec& A, std::string path);

void print_vector(const dvec& src, const char* prefix = "", const char* postfix = "", const char* delim = ",");
void print_vector(const lvec& src, const char* prefix = "", const char* postfix = "", const char* delim = ",");

bool approx_eq(const dvec& a, const dvec& b, double tol = 1e-10);
bool approx_eq(const dmat& a, const dmat& b, double tol = 1e-10);