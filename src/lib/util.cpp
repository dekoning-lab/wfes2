#include "util.hpp"

llong positive_min(llong a, llong b)
{
    if (a == 0 && b == 0) return 0;
    if (a >= 0 && b <  0) return a;
    if (a <  0 && b >= 0) return b;
    if (a <= b) return a;
    else return b;
}
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

void write_vector_map_to_file(const std::map<llong, dvec>& A, std::string path, bool append) {
    if (path == "stdout") {
	for(auto const &item : A) {
	    llong key = item.first;
	    dvec value = item.second;
	    std::cout << key << ", " << value.format(CSVRowFormat) << std::endl;
	}
    } else {
	std::ios_base::openmode mode = append ? std::ios_base::app : std::ios_base::out;
	std::ofstream file(path, mode);
	if (file.is_open()) {
	    for(auto const &item : A) {
		llong key = item.first;
		dvec value = item.second;
		file << key << ", " << value.format(CSVRowFormat) << std::endl;
	    }
	}
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

// Pick start indeces for a given vector of population sizes
// 100, 200 -> 0, 199
lvec start_indeces(lvec n) {
    lvec si = lvec::Zero(n.size());
    for (llong i = 1; i < n.size(); i++) {
	si[i] = si[i-1] + n[i-1];
    }
    return si;
}

// Like python's strided slice array[a:b:s]
lvec range_step(llong a, llong b, llong s) {
    lvec r = lvec::Zero(ceil((b-a)/double(s)));
    for(llong v = a, i = 0; v < b; v += s, i++) {
	r[i] = v;
    }
    return r;
}
