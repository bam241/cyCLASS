#include <gtest/gtest.h>

#include "bu_solver_mlp.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"

using cybam::BU_solver_MLP;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class BU_solver_MLPTest : public ::testing::Test {
 protected:
  cyclus::TestContext tc;
  BU_solver_MLP* facility;

  virtual void SetUp() {
    facility = new BU_solver_MLP(tc.get());
  }

  virtual void TearDown() {
    delete facility;
  }
};

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(BU_solver_MLPTest, InitialState) {
  // Test things about the initial state of the facility here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(BU_solver_MLPTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
  // Test BU_solver_MLP specific aspects of the print method here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(BU_solver_MLPTest, Tick) {
  ASSERT_NO_THROW(facility->Tick());
  // Test BU_solver_MLP specific behaviors of the Tick function here
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(BU_solver_MLPTest, Tock) {
  EXPECT_NO_THROW(facility->Tock());
  // Test BU_solver_MLP specific behaviors of the Tock function here
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
// Do Not Touch! Below section required for connection with Cyclus
cyclus::Agent* BU_solver_MLPConstructor(cyclus::Context* ctx) {
  return new BU_solver_MLP(ctx);
}
// Required to get functionality in cyclus agent unit tests library
#ifndef CYCLUS_AGENT_TESTS_CONNECTED
int ConnectAgentTests();
static int cyclus_agent_tests_connected = ConnectAgentTests();
#define CYCLUS_AGENT_TESTS_CONNECTED cyclus_agent_tests_connected
#endif  // CYCLUS_AGENT_TESTS_CONNECTED
INSTANTIATE_TEST_CASE_P(BU_solver_MLP, FacilityTests,
                        ::testing::Values(&BU_solver_MLPConstructor));
INSTANTIATE_TEST_CASE_P(BU_solver_MLP, AgentTests,
                        ::testing::Values(&BU_solver_MLPConstructor));
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
