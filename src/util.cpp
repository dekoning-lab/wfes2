#include "util.hpp"

bool approx_eq(const dvec& a, const dvec& b, double tol) {
	return ((a - b).array().abs() < tol).all();
}

bool approx_eq(const dmat& a, const dmat& b, double tol) {
	return ((a - b).array().abs() < tol).all();
}

void write_matrix_to_file(const dmat& A, std::string path) {
	if (path == "stdout") {
		std::cout << A.format(CSVFormat) << std::endl;
	} else {
		std::ofstream file(path);
    	if (file.is_open()) file << A.format(CSVFormat) << std::endl;	
	}
}

void write_vector_to_file(const dvec& A, std::string path) {
	if (path == "stdout") {
		std::cout << A.format(CSVFormat) << std::endl;
	} else {
		std::ofstream file(path);
		if (file.is_open()) file << A.format(CSVFormat) << std::endl;
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