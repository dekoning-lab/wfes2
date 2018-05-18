#pragma once

#include "common.hpp"

// Tempalte sto<xx> functions:
template<typename T>
T from_string(std::string const& str, size_t* pos = 0);

struct TokenReader
{
    void split_string(const std::string &s, char delim, std::back_insert_iterator<std::deque<std::string>> result);
    std::deque<std::string> split(const std::string &s, char delim = ',');
};

template<typename T>
struct NumericVectorReader : TokenReader
{
    void operator()(const std::string& name, const std::string& value, Eigen::Matrix<T, Eigen::Dynamic, 1>& destination, char delim = ',') {
        std::deque<std::string> tokens = split(value, delim);
        destination.resize(tokens.size());
        for(size_t i = 0; i < tokens.size(); i++) {
            destination[i] = from_string<T>(tokens[i]);
        }
    }
};

template<typename T>
struct NumericMatrixReader : TokenReader
{
    void operator()(const std::string& name, const std::string& value, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& destination, char delim = ',') {
        std::deque<std::string> rows = split(value, ';');
        size_t n_row = rows.size();
        std::deque<std::deque<std::string>> m(n_row);
        for(size_t i = 0; i < n_row; i++) {
            m[i] = split(rows[i], delim);
        }

        destination.resize(n_row, n_row);
        for(size_t i = 0; i < n_row; i++) {
            for(size_t j = 0; j < n_row; j++) {
                destination(i, j) = from_string<T>(m[i][j]);
            }
        }
    }
};

template<typename T>
const char* string_format();

template<typename T>
void print_vector(const Eigen::Matrix<T, Eigen::Dynamic, 1>& src, const char* prefix = "", const char* postfix = "", const char* delim = ",") {
    size_t size = src.size();
    const char* f = string_format<T>();
    printf("%s", prefix);
    for(size_t i = 0; i < size - 1; i++) {
        printf(f, src(i));
        printf("%s ", delim);
    }
    printf(f, src(size-1));
    printf("%s", postfix);
}

const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

template<typename T>
void write_matrix_to_file(const Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& A, std::string path) {
    std::ofstream file(path);
    if (file.is_open()) file << A.format(CSVFormat) << std::endl;
}


template<typename T>
void write_vector_to_file(const Eigen::Matrix<T, Eigen::Dynamic, 1>& A, std::string path) {
    std::ofstream file(path);
    if (file.is_open()) file << A.format(CSVFormat) << std::endl;
}
