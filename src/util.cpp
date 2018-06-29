#include "util.hpp"

bool approx_eq(const dvec& a, const dvec& b, double tol) {
	return ((a - b).array().abs() < tol).all();
}

bool approx_eq(const dmat& a, const dmat& b, double tol) {
	return ((a - b).array().abs() < tol).all();
}

double total_diff(const dmat& a, const dmat& b) {
    return (a - b).array().abs().sum();
}

void write_matrix_to_file(const dmat& A, std::string path, bool append) {
	if (path == "stdout") {
		std::cout << A.format(CSVFormat) << std::endl;
	} else {
		std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;
		std::ofstream file(path, mode);
    	if (file.is_open()) file << A.format(CSVFormat) << std::endl;	
	}
}

void write_vector_to_file(const dvec& A, std::string path, bool append) {
	if (path == "stdout") {
		std::cout << A.format(CSVRowFormat) << std::endl;
	} else {
		std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;
		std::ofstream file(path, mode);
		if (file.is_open()) file << A.format(CSVRowFormat) << std::endl;
	}
}

void print_vector(const dvec& src, const char* prefix, const char* postfix, const char* delim) {
    size_t size = src.size();
    printf("%s", prefix);
    for(size_t i = 0; i < size - 1; i++) {
        printf(DPF, src(i));
        printf("%s ", delim);
    }
    printf(DPF, src(size-1));
    printf("%s", postfix);
}

void print_vector(const lvec& src, const char* prefix, const char* postfix, const char* delim) {
    size_t size = src.size();
    printf("%s", prefix);
    for(size_t i = 0; i < size - 1; i++) {
        printf(LPF, src(i));
        printf("%s ", delim);
    }
    printf(LPF, src(size-1));
    printf("%s", postfix);
}