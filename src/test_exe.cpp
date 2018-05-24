#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "WrightFisherElementary.hpp"
#include "WrightFisher.hpp"
#include "util.hpp"

namespace WF = WrightFisher;
namespace WFE = WrightFisherElementary;

using namespace std;

TEST_CASE("Psi Calculations", "[psi]") {
	REQUIRE(WF::psi_diploid(0, 100) == WFE::psi_diploid(0, 100));
}

TEST_CASE("Binomial row", "[binom-row]") {
	dvec a = WFE::binom_row(100, 0.5);
	dvec b = WF::binom_row(100, 0.5, 0).Q;
	REQUIRE(approx_eq(a, b));
}

// SINGLE 

SCENARIO("Single Square matrix are built correctly", "[single-square]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(100, 100, a_t, 0, 0.5, 1e-9, 1e-9, 0);
			std::pair<dmat,dmat> W2 = WFE::Single(100, 100, a_t, 0, 0.5, 1e-9, 1e-9);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));	
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
		}
	}
	
}

SCENARIO("Single Wide matrix are built correctly", "[single-wide]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(200, 100, a_t, 0, 0.5, 1e-9, 1e-9, 0);
			std::pair<dmat, dmat> W2 = WFE::Single(200, 100, a_t, 0, 0.5, 1e-9, 1e-9);
			REQUIRE(approx_eq(W1.Q.dense(), W2.first));	
			if(a_t != WF::NON_ABSORBING) {
				REQUIRE(approx_eq(W1.R, W2.second));
			}
		}
	}
}

SCENARIO("Single Narrow matrix are built correctly", "[single-narrow]") {
	for(int i = WF::NON_ABSORBING; i <= WF::BOTH_ABSORBING; i ++) {
		WF::absorption_type a_t = WF::absorption_type(i);
		GIVEN(absorption_type_desc(a_t)) {
			WF::Matrix W1 = WF::Single(100, 200, a_t, 0, 0.5, 1e-9, 1e-9, 0);
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
	// REQUIRE(approx_eq(W1.Q.dense(), W2.first));
}

SCENARIO("Switching Up matrix are built correctly", "[switch-up]") {
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

SCENARIO("Switching Down matrix are built correctly", "[switch-down]") {
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