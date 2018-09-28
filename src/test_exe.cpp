#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "WrightFisherElementary.hpp"
#include "WrightFisher.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
namespace WFE = WrightFisherElementary;

using namespace std;

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
    llong k = 100;
    SparseMatrix A(k, k);
    dmat B(k, k);
    for(llong i = 0; i < k; i++) {
        dvec range = dvec::LinSpaced(k, 1, 100);
        A.append_row(range, 0, k-1);
        B.row(i) = range;
    }
    A.compress_csr();

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
			WF::Matrix W1 = WF::Switching(N, a_t, s, h, mu, mu, sw, 0);
			std::pair<dmat, dmat> W2 = WFE::SwitchingTwoByTwo(N, a_t, s, h, mu, mu, sw);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
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

int main( int argc, char* argv[] ) {
	#ifdef OMP
		cout << "Using OMP" << endl;
		omp_set_num_threads(2);
	#endif

	mkl_set_num_threads(2);

	int result = Catch::Session().run( argc, argv );


	return result;
}
