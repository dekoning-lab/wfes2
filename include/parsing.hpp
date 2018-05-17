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
const char* string_format();