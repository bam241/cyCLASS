#include <gtest/gtest.h>

#include "bu_solver_mlp.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"

using cybam::MLPBUsolver;

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
class BU_solver_MLPTest : public ::testing::Test {
 protected:
  MLPBUsolver* facility;

  virtual void SetUp() {
    facility = new MLPBUsolver();
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

