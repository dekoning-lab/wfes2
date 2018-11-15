#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "WrightFisherElementary.hpp"
#include "WrightFisher.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
namespace WFE = WrightFisherElementary;

using namespace std;

TEST_CASE("Positive min", "[util]")
{
    REQUIRE(positive_min(-1,0) == 0);
    REQUIRE(positive_min(0,-1) == 0);
    REQUIRE(positive_min(3,-1) == 3);
    REQUIRE(positive_min(-1,7) == 7);
    REQUIRE(positive_min(2,4)  == 2);
    REQUIRE(positive_min(8,9)  == 8);
}

TEST_CASE("Sparse Matrix Column accessor", "[sparse]") {
	WF::Matrix W = WF::Single(100, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, true, 0);
	REQUIRE(approx_eq(W.Q.dense().col(10), W.Q.col_copy(10)));
	REQUIRE(approx_eq(W.Q.dense().col(0), W.Q.col_copy(0)));
	REQUIRE(approx_eq(W.Q.dense().col(198), W.Q.col_copy(198)));
}

TEST_CASE("Sparse Matrix subtract_identity is inverse of itself", "[sparse]") {
	WF::Matrix W1 = WF::Single(100, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, true, 0);
	WF::Matrix W2 = WF::Single(100, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, true, 0);

	W1.Q.subtract_identity();
	REQUIRE(!W1.Q.approx_eq(W2.Q));	

	W1.Q.subtract_identity();
	REQUIRE(W1.Q.approx_eq(W2.Q));	
}

TEST_CASE("Sparse Matrix vector multiplication", "[sparse]") {
    llong N = 100;
	WF::Matrix W = WF::Single(N, N, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, true, 0);

    dvec a = dvec::Ones((2*N)-1).array().abs();
    dvec b = a;

    W.Q.multiply_inplace_rep(a, 100, true);

    for(llong i = 0; i < 100; i++) {
        b = W.Q.multiply(b, true);
    }
    REQUIRE(approx_eq(a, b));
}

TEST_CASE("Sparse Matrix in COO format converts to CSR", "[sparse]") 
{
    llong k = 4;
    SparseMatrix A(k, k);
    dmat B(k, k);
    for(llong i = 0; i < k; i++) {
        dvec range = dvec::LinSpaced(k, 1, 100);
        A.append_row(range, 0, k-1);
        B.row(i) = range;
    }

    REQUIRE(approx_eq(A.dense(), B));
}


SCENARIO("Runs with large magnitude selection", "[selection]") 
{
    llong N = 100;
    dvec sel(4); sel << -0.99, -0.9, 5.7, 10.99;
    for (int i = 0; i < 4; i++) {
        GIVEN("Selection coefficient is " + to_string(sel(i))) {
            WF::Matrix W = WF::Single(N, N, WF::NON_ABSORBING, sel(i));

            WHEN("Building a Single model matrix") {
                W.Q.get_diag_copy();
                REQUIRE(W.Q.get_diag_copy().size() == (2 * N) + 1);
            }

            WHEN("Building an equilibrium matrix") {
                WF::Matrix E = WF::Equilibrium(N, sel(i));
                REQUIRE(E.Q.get_diag_copy().size() == (2 * N) + 1);
            }
        }
    }
}

TEST_CASE("Psi Calculations", "[binom]") {
	REQUIRE(WF::psi_diploid(0, 100) == WFE::psi_diploid(0, 100));
}

TEST_CASE("Binomial row", "[binom]") {
	dvec a = WFE::binom_row(100, 0.5);
	dvec b = WF::binom_row(100, 0.5, 0).Q;
	REQUIRE(approx_eq(a, b));
}

// SINGLE 

SCENARIO("Single Square matrix are built correctly", "[single]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(100, 100, a_t, 0, 0.5, 1e-9, 1e-9, true, 0);
			std::pair<dmat,dmat> W2 = WFE::Single(100, 100, a_t, 0, 0.5, 1e-9, 1e-9);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));	
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
		}
	}
	
}

SCENARIO("Single Wide matrix are built correctly", "[single]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(200, 100, a_t, 0, 0.5, 1e-9, 1e-9, true, 0);
			std::pair<dmat, dmat> W2 = WFE::Single(200, 100, a_t, 0, 0.5, 1e-9, 1e-9);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));	
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
		}
	}
}

SCENARIO("Single Narrow matrix are built correctly", "[single]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(100, 200, a_t, 0, 0.5, 1e-9, 1e-9, true, 0);
			std::pair<dmat, dmat> W2 = WFE::Single(100, 200, a_t, 0, 0.5, 1e-9, 1e-9);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
		}
	}
}

SCENARIO("Equilibrium linear solution matrix is built correctly", "[equilibrium]") {
	WF::Matrix W1 = WF::Equilibrium(100, 0, 0.5, 1e-9, 1e-9, 0);
	std::pair<dmat, dmat> W2 = WFE::Equilibrium(100, 0, 0.5, 1e-9, 1e-9);
	REQUIRE(approx_eq(W1.Q.dense(), W2.first));
}

SCENARIO("Soft Selective Sweep model matrix is built correctly", "[soft]") {
    dvec s(2); s << 0, 0.1;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1 - 1e-3, 1e-3, 0, 1;

	WF::Matrix W1 = WF::NonAbsorbingToFixationOnly(100, s, h, mu, mu, sw, 0);
	std::pair<dmat, dmat> W2 = WFE::NonAbsorbingToFixationOnly(100, s, h, mu, mu, sw);

	// cout << W1.Q << endl;
	// cout << "--" << endl;
	// cout << W1.Q.dense() << endl;
	// cout << "--" << endl;
	// cout << W2.first << endl;

	REQUIRE(approx_eq(W1.Q.dense(), W2.first));
}

SCENARIO("Soft Selective Sweep Both Absorbing model matrix is built correctly", "[soft]") {
    dvec s(2); s << 0, 0.1;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1 - 1e-3, 1e-3, 0, 1;

	WF::Matrix W1 = WF::NonAbsorbingToBothAbsorbing(100, s, h, mu, mu, sw, 0);
	std::pair<dmat, dmat> W2 = WFE::NonAbsorbingToBothAbsorbing(100, s, h, mu, mu, sw);

	// cout << W1.Q << endl;
	// cout << "--" << endl;
	// cout << W1.Q.dense() << endl;
	// cout << "--" << endl;
	// cout << W2.first << endl;

	REQUIRE(approx_eq(W1.Q.dense(), W2.first));

	// cout << W1.R << endl;
	// cout << "--" << endl;
	// cout << W2.second << endl;
	// REQUIRE(approx_eq(W1.R), W2.second);
}

SCENARIO("Switching Up matrix are built correctly", "[switch]") {
	lvec N(2); N << 100, 200;
	dvec s(2); s << 0, 0.1;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1, 1, 1, 1;

	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {

		WF::absorption_type a_t = WF::absorption_type(i);

		GIVEN(absorption_type_desc(a_t)) {

			WF::Matrix W1            = WF::Switching(N, a_t, s, h, mu, mu, sw, 0);
			std::pair<dmat, dmat> W2 = WFE::SwitchingTwoByTwo(N, a_t, s, h, mu, mu, sw);

			REQUIRE(approx_eq(W1.Q.dense(), W2.first));
			if(a_t != WF::NON_ABSORBING) REQUIRE(approx_eq(W1.R, W2.second));
		}
	}
}

SCENARIO("Switching Down matrix are built correctly", "[switch]") {
	lvec N(2); N << 200, 100;
	dvec s(2); s << 0, 0.1;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1, 1, 1, 1;

	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Switching(N, a_t, s, h, mu, mu, sw, 0);
			std::pair<dmat, dmat> W2 = WFE::SwitchingTwoByTwo(N, a_t, s, h, mu, mu, sw);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}	
		}
	}
}

SCENARIO("Switching runs with large magnitude selection", "[switching-selection]") 
{
	lvec N(2); N << 200, 100;
    llong size = (2 * N(0) + 1) + (2 * N(1) + 1);
    dvec sel(4); sel << -0.99, -0.9, 5.7, 10.99;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1, 1, 1, 1;

    for (int i = 0; i < 4; i++) {
        GIVEN("Selection coefficient is " + to_string(sel(i))) {

            dvec s(2); s << sel(i), sel(i);
            WF::Matrix W = WF::Switching(N, WF::NON_ABSORBING, s, h, mu, mu, sw);

            WHEN("Building a Switching model matrix") {
                W.Q.get_diag_copy();
                REQUIRE(W.Q.get_diag_copy().size() == size);
            }
        }
    }
}

SCENARIO("Switching non-standard models with large magnitude selection", "[switching-selection]") 
{
    llong N = 100;
    dvec sel(4); sel << -0.99, -0.9, 5.7, 10.99;
	dvec h(2); h << 0.5, 0.5;
	dvec mu(2); mu << 1e-9, 1e-9;
	dmat sw(2,2); sw << 1, 1, 1, 1;

    for (int i = 0; i < 4; i++) {
        GIVEN("Selection coefficient is " + to_string(sel(i))) {

            dvec s(2); s << sel(i), sel(i);
            WF::Matrix W1 = WF::NonAbsorbingToFixationOnly(100, s, h, mu, mu, sw, 0);
            WHEN("Building a N->F model matrix") {
                W1.Q.get_diag_copy();
                REQUIRE(W1.Q.get_diag_copy().size() == 4*N+1);
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        GIVEN("Selection coefficient is " + to_string(sel(i))) {

            dvec s(2); s << sel(i), sel(i);
            WF::Matrix W2 = WF::NonAbsorbingToBothAbsorbing(100, s, h, mu, mu, sw, 0);
            WHEN("Building a N->B model matrix") {
                W2.Q.get_diag_copy();
                REQUIRE(W2.Q.get_diag_copy().size() == 4*N);
            }
        }
    }
}

TEST_CASE("start_indeces are computed", "[extra]") {
    lvec x(3); x << 199, 199, 399;
    lvec ex(3); ex << 0, 199, 398;
    REQUIRE(ex == start_indeces(x));
}

TEST_CASE("start_indeces with arithmetic", "[extra]") {
    lvec x(3); x << 100, 100, 200;
    lvec ex(3); ex << 0, 199, 398;
    REQUIRE(ex == start_indeces(2*x-lvec::Ones(x.size())));
}

TEST_CASE("boolean subset", "[extra]") {
    dvec x(6); x << 0.9, 0.8, 0.7, 0.6, 0.5, 0.4;
    dvec ex(3); ex << 0.9, 0.8, 0.7;
    REQUIRE((x.array()>0.5).count() == 4);
}

TEST_CASE("step range works", "[extra]") {
    lvec ex0(3); ex0 << 0, 2, 4;
    lvec ex1(3); ex1 << 1, 3, 5;
    REQUIRE(range_step(0,6,2) == ex0);
    REQUIRE(range_step(1,6,2) == ex1);
}

int main( int argc, char* argv[] ) {
	#ifdef OMP
		cout << "Using OMP" << endl;
		omp_set_num_threads(2);
	#endif

	mkl_set_num_threads(2);

	int result = Catch::Session().run( argc, argv );


	return result;
}
