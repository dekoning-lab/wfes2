#pragma once

#include "common.hpp"

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");
const static Eigen::IOFormat CSVRowFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", ", ");

void write_matrix_to_file(const dmat& A, std::string path, bool append = false);
void write_vector_map_to_file(const std::map<llong,dvec>& A, std::string path, bool append = false);
void write_vector_to_file(const dvec& A, std::string path, bool append = false);

void print_vector(const dvec& src, const char* prefix = "", const char* postfix = "", const char* delim = ",");
void print_vector(const lvec& src, const char* prefix = "", const char* postfix = "", const char* delim = ",");

llong positive_min(llong a, llong b);
bool approx_eq(const dvec& a, const dvec& b, double tol = 1e-10);
bool approx_eq(const dmat& a, const dmat& b, double tol = 1e-10);
double total_diff(const dmat& a, const dmat& b);

lvec start_indeces(lvec n);
lvec range_step(llong a, llong b, llong s);
