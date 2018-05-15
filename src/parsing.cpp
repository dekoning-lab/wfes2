#include "parsing.hpp"

template<>
llong from_string<llong>(std::string const& str, size_t* pos) {
    return stoll(str, pos, 10);
}

template<>
int from_string<int>(std::string const& str, size_t* pos) {
    return stoi(str, pos, 10);
}

template<>
double from_string<double>(std::string const& str, size_t* pos) {
    return stod(str, pos);
}

void TokenReader::split_string(const std::string &s, char delim, std::back_insert_iterator<std::deque<std::string>> result) {
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::deque<std::string> TokenReader::split(const std::string &s, char delim) {
    std::deque<std::string> elems;
    split_string(s, delim, back_inserter(elems));
    return elems;
}

template <>
const char* string_format<double>() { return DPF; }

template <>
const char* string_format<llong>() { return LPF; }