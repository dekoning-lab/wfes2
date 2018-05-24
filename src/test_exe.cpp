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

// Non-absorbing 
TEST_CASE("Non-Absorbing WF - square", "[non-absorbing]") {
	WF::Matrix W1 = WF::Single(100, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Non-Absorbing WF - wide", "[non-absorbing]") {
	WF::Matrix W1 = WF::Single(100, 20, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 20, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Non-Absorbing WF - narrow", "[non-absorbing]") {
	WF::Matrix W1 = WF::Single(20, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(20, 100, WF::NON_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

// Both absorbing
TEST_CASE("Both-Absorbing WF - square", "[both-absorbing]") {
	WF::Matrix W1 = WF::Single(100, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Both-Absorbing WF - wide", "[both-absorbing]") {
	WF::Matrix W1 = WF::Single(100, 20, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 20, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Both-Absorbing WF - narrow", "[both-absorbing]") {
	WF::Matrix W1 = WF::Single(20, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(20, 100, WF::BOTH_ABSORBING, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

// Fixation

TEST_CASE("Fixation WF - square", "[fixation]") {
	WF::Matrix W1 = WF::Single(100, 100, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 100, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Fixation WF - wide", "[fixation]") {
	WF::Matrix W1 = WF::Single(100, 20, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 20, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Fixation WF - narrow", "[fixation]") {
	WF::Matrix W1 = WF::Single(20, 100, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(20, 100, WF::FIXATION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

//Extinction

TEST_CASE("Extinction WF - square", "[extinction]") {
	WF::Matrix W1 = WF::Single(100, 100, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 100, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Extinction WF - wide", "[extinction]") {
	WF::Matrix W1 = WF::Single(100, 20, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(100, 20, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}

TEST_CASE("Extinction WF - narrow", "[extinction]") {
	WF::Matrix W1 = WF::Single(20, 100, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9, 0);
	dmat W2 = WFE::Single(20, 100, WF::EXTINCTION_ONLY, 0, 0.5, 1e-9, 1e-9).first;
	REQUIRE(approx_eq(W1.Q.dense(), W2));
}