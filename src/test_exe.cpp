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

TEST_CASE("Non-Absorbing WF - square", "[non-absorbing]") {
	WF::Matrix W1 = WF::SingleAlt(100, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Non-Absorbing WF - wide", "[non-absorbing]") {
	WF::Matrix W1 = WF::SingleAlt(100, 20, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 20, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Non-Absorbing WF - narrow", "[non-absorbing]") {
	WF::Matrix W1 = WF::SingleAlt(20, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(20, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}