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

TEST(BU_solver_MLPTest, UpdateInputComposition) {

    CompMap fissil;
    fissil.insert(std::pair<Nuc, double>(942380000,94238));
    fissil.insert(std::pair<Nuc, double>(942390000,94239));
    fissil.insert(std::pair<Nuc, double>(942400000,94240));
    fissil.insert(std::pair<Nuc, double>(942410000,94241));
    fissil.insert(std::pair<Nuc, double>(942420000,94242));
    fissil.insert(std::pair<Nuc, double>(952410000,95241));

    CompMap fertil;
    fertil.insert(std::pair<Nuc, double>(922350000,92235));
    fertil.insert(std::pair<Nuc, double>(922380000,92238));

    MLPBUsolver* MyBUSolver;

    double BU = 41.235;

   // MyBUSolver->UpdateInputComposition(cyclus::Composition::CreateFromAtom(fissil),
     //                                  cyclus::Composition::CreateFromAtom(fertil),
       //                                BU);


}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST_F(BU_solver_MLPTest, Print) {
  EXPECT_NO_THROW(std::string s = facility->str());
  // Test BU_solver_MLP specific aspects of the print method here
}

