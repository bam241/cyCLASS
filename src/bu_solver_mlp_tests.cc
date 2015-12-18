#include <gtest/gtest.h>

#include "bu_solver_mlp.h"

#include "agent_tests.h"
#include "context.h"
#include "facility_tests.h"

using cybam::MLPBUsolver;
namespace cybam {


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(BU_solver_MLPTest, InitialState) {
  // Test things about the initial state of the facility here
}

TEST(BU_solver_MLPTest, GetEnrichment) {

    CompMap fissil_comp;
    fissil_comp.insert(std::pair<Nuc, double>(942380000,0.28));
    fissil_comp.insert(std::pair<Nuc, double>(942390000,4.65));
    fissil_comp.insert(std::pair<Nuc, double>(942400000,2.19));
    fissil_comp.insert(std::pair<Nuc, double>(942410000,1.06));
    fissil_comp.insert(std::pair<Nuc, double>(942420000,0.73));
    fissil_comp.insert(std::pair<Nuc, double>(952410000,0.1));
    fissil_comp =  NormalizeComp(fissil_comp);

    CompMap fertil_comp;
    fertil_comp.insert(std::pair<Nuc, double>(922350000,0.25));
    fertil_comp.insert(std::pair<Nuc, double>(922380000,99.75));
    fertil_comp =  NormalizeComp(fertil_comp);

    MLPBUsolver* MyBUSolver = new MLPBUsolver();

    double BU = 41.235;

    double Enrch = MyBUSolver->GetEnrichment(Composition::CreateFromAtom(fissil_comp),
                                             Composition::CreateFromAtom(fertil_comp),
                                             BU );


    std::cout << Enrch << std::endl;
    ASSERT_LE(0, Enrch);
    ASSERT_GE(1, Enrch);
    EXPECT_GT(0.2, Enrch);
}

TEST(BU_solver_MLPTest, GetBu) {

    CompMap fissil_comp;
    fissil_comp.insert(std::pair<Nuc, double>(942380000,0.28));
    fissil_comp.insert(std::pair<Nuc, double>(942390000,4.65));
    fissil_comp.insert(std::pair<Nuc, double>(942400000,2.19));
    fissil_comp.insert(std::pair<Nuc, double>(942410000,1.06));
    fissil_comp.insert(std::pair<Nuc, double>(942420000,0.73));
    fissil_comp.insert(std::pair<Nuc, double>(952410000,0.1));
    fissil_comp =  NormalizeComp(fissil_comp);

    CompMap fertil_comp;
    fertil_comp.insert(std::pair<Nuc, double>(922350000,0.25));
    fertil_comp.insert(std::pair<Nuc, double>(922380000,99.75));
    fertil_comp =  NormalizeComp(fertil_comp);

    MLPBUsolver* MyBUSolver = new MLPBUsolver();

    float BU = 41.235;

    double Enrch = MyBUSolver->GetEnrichment(Composition::CreateFromAtom(fissil_comp),
                                             Composition::CreateFromAtom(fertil_comp),
                                             BU );

    Composition::Ptr fuel = Composition::CreateFromAtom( fertil_comp*( 1-Enrch ) + fissil_comp*Enrch );

    float BU_cal = MyBUSolver->GetBU(fuel);

    std::cout << BU << " " << BU_cal << std::endl;

    ASSERT_LT(std::abs(BU - BU_cal)/(BU_cal/2.+BU/2.), 1e-5);
}



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TEST(BU_solver_MLPTest, Print) {
    MLPBUsolver* MyBUSolver = new MLPBUsolver();

  EXPECT_NO_THROW(std::string s = MyBUSolver->str());
  // Test BU_solver_MLP specific aspects of the print method here
}
} // namespace cybam

