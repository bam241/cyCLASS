#include "CLASSAdaptator.h"

#include "EquivalenceModel.hxx"
#include "EvolutionData.hxx"
#include "EvolutionData.hxx"
#include "IrradiationModel.hxx"
#include "XSModel.hxx"

#include "IM_Matrix.hxx"
#include "IM_RK4.hxx"

#include "XSM_MLP.hxx"

#include "EQ_OneParameter.hxx"

#include <string>

using cyclus::Nuc;
using cyclus::Material;
using cyclus::Composition;

namespace cyclass {

EquivalenceModel* EQmodelfor(std::string name, std::string command) {
  std::stringstream command_t;
  command_t << command;

  if (name == "EQmodelfor") {
    std::string WeightPath;
    getline(command_t, WeightPath, ',');

    std::string InformationFile;
    getline(command_t, InformationFile, ',');

    std::string buff;
    getline(command_t, buff, ',');
    int NumOfBatch = atoi(buff.c_str());

    getline(command_t, buff, ',');
    double CriticalityThreshold = atof(buff.c_str());

    EQ_OneParameter* myEQM = new EQ_OneParameter(WeightPath, InformationFile);
    myEQM->SetModelParameter("kThreshold", CriticalityThreshold);
    myEQM->SetModelParameter("NumberOfBatch", NumOfBatch);
    return myEQM;
  }
}

XSModel* XSmodelfor(std::string name, std::string command) {
  std::stringstream command_t;
  command_t << command;

  if (name == "XSM_MLP") {
    std::string buff;
    getline(command_t, buff, ',');
    string TMVA_Weight_Dir = buff;
    getline(command_t, buff, ',');
    string Info = buff;
    getline(command_t, buff, ',');

    return new XSM_MLP(TMVA_Weight_Dir.c_str(), Info.c_str());

  } else {
    stringstream msg;
    msg << "Bad XSModel keyword " << name;
    throw cyclus::ValidationError(msg.str());
  }
}

IrradiationModel* IMmodelfor(std::string name, std::string command) {
  if (name == "RK4") {
    return new IM_RK4();
  } else if (name == "MATRIX") {
    return new IM_Matrix();
  } else {
    stringstream msg;
    msg << "Bad IRModel keyword " << name;
    throw cyclus::ValidationError(msg.str());
  }
}

CLASSAdaptator::CLASSAdaptator(std::string EQModel, std::string EQcommand) {
  
  myPhysicsModel = new PhysicsModels();

  myPhysicsModel->SetEQM(EQmodelfor(EQModel, EQcommand));

  //TMVAWeightFile = EQcommand;
  //IsotopicVector IV_fissil =
  //    myPhysicsModel->GetEQM()->GetStreamList("Fissile");
 // if (IV_fissil.GetZAIIsotopicQuantity(94, 241, 0) > 0)
 //   IV_fissil += ZAI(95, 241, 0) * 1;

 // fissil_list = Composition::CreateFromAtom(CLASS2CYCLUS(IV_fissil));

 // fertil_list = Composition::CreateFromAtom(CLASS2CYCLUS(
 //     myPhysicsModel->GetEQM()->GetStreamList("Fertile")));
  
}

CLASSAdaptator::CLASSAdaptator(std::string EQModel, std::string EQcommand,
                               std::string XSModel, std::string XScommand,
                               std::string IMModel, std::string IMcommand) {
  
  TMVAWeightFile = EQcommand;
  myPhysicsModel = new PhysicsModels();
  

  myPhysicsModel->SetEQM(EQmodelfor(EQModel, EQcommand));
  myPhysicsModel->SetXSM(XSmodelfor(XSModel, XScommand));
  myPhysicsModel->SetIM(IMmodelfor(IMModel, IMcommand));
  

  
 // IsotopicVector IV_fissil =
 //     myPhysicsModel->GetEQM()->GetStreamList("Fissile");
 // if (IV_fissil.GetZAIIsotopicQuantity(94, 241, 0) > 0)
 //   IV_fissil += ZAI(95, 241, 0) * 1;

 // 

 // fissil_list = Composition::CreateFromAtom(CLASS2CYCLUS(IV_fissil));
 // 

 // fertil_list = Composition::CreateFromAtom(CLASS2CYCLUS(
 //     myPhysicsModel->GetEQM()->GetStreamList("Fertile")));
  
}

//________________________________________________________________________
float CLASSAdaptator::GetEnrichment(cyclus::Composition::Ptr c_fissil,
                                    cyclus::Composition::Ptr c_fertil,
                                    double bu_target, double eps) const {
  

  double val = 0.20;

  IsotopicVector iv_fissil = CYCLUS2CLASS(c_fissil);
  IsotopicVector iv_fertil = CYCLUS2CLASS(c_fertil);
  
  IsotopicVector iv_fuel = iv_fissil * val + (1-val) * iv_fertil;  
  
  float bu = myPhysicsModel->GetEQM()->CalculateTargetParameter(iv_fuel, "BurnUpMax");
  float val_p = val;
  float bu_p = bu;

  val = 0.03;
  iv_fuel = iv_fissil * val + (1-val) * iv_fertil;  
  bu = myPhysicsModel->GetEQM()->CalculateTargetParameter(iv_fuel, "BurnUpMax");
  int i = 0;
  while( abs(bu-bu_target) > eps ){
      // BU  = A * frac + B;
      float d_bu = bu-bu_target;
      double A = (bu - bu_p) / (val - val_p);
      double B = bu - A * val;
      if( bu == bu_p) return val;
      //old = new
      bu_p = bu;
      val_p = val;

      //target bu = bu_target
      val = (bu_target - B) / A;
      iv_fuel = iv_fissil * val + (1-val) * iv_fertil;
      bu = myPhysicsModel->GetEQM()->CalculateTargetParameter(iv_fuel, "BurnUpMax");
      i++;
  }
  return val;  //
}

float CLASSAdaptator::GetBU(cyclus::Composition::Ptr fuel, double eps) const{
  IsotopicVector iv_fuel = CYCLUS2CLASS(fuel);
  float bu = myPhysicsModel->GetEQM()->CalculateTargetParameter(iv_fuel, "BurnUpMax");
  return bu;
}


cyclus::Composition::Ptr CLASSAdaptator::GetCompAfterIrradiation(
    cyclus::Composition::Ptr InitialCompo, double power, double mass,
    double burnup) {
  IsotopicVector InitialIV = CYCLUS2CLASS(InitialCompo);
  double ratio = 1 / InitialIV.GetTotalMass() * mass * 1e-3;
  InitialIV *= ratio;

  cSecond finaltime = burnup * mass * 1e-3 / (power * 1e-3) * 3600 * 24;
  // cout << "time " << finaltime << endl;
  // cout << "power " << power << endl;
  EvolutionData myEvolution =
      myPhysicsModel->GenerateEvolutionData(InitialIV, finaltime, power * 1e6);
  IsotopicVector AfterIrradiationIV =
      myEvolution.GetIsotopicVectorAt(finaltime).GetActinidesComposition();
  double CLASS_mass_ratio = AfterIrradiationIV.GetTotalMass() /
                            InitialIV.GetActinidesComposition().GetTotalMass();
  double missing_mass = -(AfterIrradiationIV.GetTotalMass() -
                          InitialIV.GetActinidesComposition().GetTotalMass());
  //:InitialIV.Print();
  // AfterIrradiationIV.Print();
  double Avogadro = 6.02214129e23;
  AfterIrradiationIV +=
      missing_mass * Avogadro / 136.9070 * ZAI(55, 137, 0) * 1e6;

  AfterIrradiationIV *= 1 / ratio;
  return Composition::CreateFromAtom(CLASS2CYCLUS(AfterIrradiationIV));
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
double AtomIn(CompMap Source) {
  double total = 0;

  CompMap::iterator it;

  for (it = Source.begin(); it != Source.end(); it++) total += it->second;

  return total;
}

//________________________________________________________________________
double AtomIn(cyclus::Composition::Ptr Source) {
  return AtomIn(Source->atom());
};

//________________________________________________________________________
cyclus::Composition::Ptr ExtractAccordinglist(cyclus::Composition::Ptr source,
                                              cyclus::Composition::Ptr list) {
  
  // create Output Composition
  CompMap separatedCompo;

  
  // Extract Nuc map from source compo & list...
  CompMap sourceComp = source->atom();
  
  CompMap ListComp = list->atom();
  

  CompMap::iterator it;
  

  // Fill output composition
  for (it = ListComp.begin(); it != ListComp.end(); it++) {
    CompMap::iterator it2 = sourceComp.find(it->first);

    if (it2 != sourceComp.end())
      separatedCompo.insert(std::pair<Nuc, double>(it->first, it2->second));
  }
  

  return Composition::CreateFromAtom(separatedCompo);
}

//________________________________________________________________________
CompMap NormalizeComp(CompMap source, double norm) {
  double Total = AtomIn(source);
  CompMap::iterator it;
  for (it = source.begin(); it != source.end(); it++)
    it->second *= norm / Total;

  return source;
}

//________________________________________________________________________
CompMap operator+(CompMap const& IVa, CompMap const& IVb) {
  CompMap IVtmp = IVa;
  CompMap IVbtmp = IVb;
  CompMap::iterator it;

  for (it = IVbtmp.begin(); it != IVbtmp.end(); it++) {
    std::pair<CompMap::iterator, bool> IResult;
    IResult = IVtmp.insert(std::pair<Nuc, double>(it->first, it->second));
    if (!IResult.second) IResult.first->second += it->second;
  }
  return IVtmp;
}

//________________________________________________________________________
CompMap operator*(CompMap const& IVA, double F) {
  CompMap IVtmp = IVA;
  CompMap::iterator it;

  for (it = IVtmp.begin(); it != IVtmp.end(); it++) it->second *= F;

  return IVtmp;
}

//________________________________________________________________________
void Print(cyclus::Composition::Ptr compo) { Print(compo->atom()); }

//________________________________________________________________________
void Print(CompMap compo) {
  CompMap::iterator it;
  cout << "Printing Compo" << std::endl;
  for (it = compo.begin(); it != compo.end(); it++) {
    std::cout << it->first << " " << it->second << std::endl;
  }
}

//________________________________________________________________________
CompMap operator-(CompMap const& IVa, CompMap const& IVb) {
  return IVa + (-1 * IVb);
};

//________________________________________________________________________
CompMap operator/(CompMap const& IVA, double F) { return IVA * (1 / F); };

//________________________________________________________________________
CompMap operator*(double F, CompMap const& IVA) { return IVA * F; };

//________________________________________________________________________
IsotopicVector CYCLUS2CLASS(CompMap c_compo) {
  IsotopicVector IV;
  map<int, double>::iterator it;

  for (it = c_compo.begin(); it != c_compo.end(); it++) {
    int Z = (int)it->first / 10000000;
    int A = (int)(it->first / 10000) % 1000;
    int I = (int)it->first % 10000;

    IV += it->second * ZAI(Z, A, I);
  }

  return IV;
}

//________________________________________________________________________
IsotopicVector CYCLUS2CLASS(cyclus::Composition::Ptr c_compo) {
  return CYCLUS2CLASS(c_compo->atom());
}

//________________________________________________________________________
CompMap CLASS2CYCLUS(IsotopicVector const& IV) {
  CompMap myCompMap;

  map<ZAI, double> IsotopicMap = IV.GetIsotopicQuantity();
  map<ZAI, double>::iterator it;

  for (it = IsotopicMap.begin(); it != IsotopicMap.end(); it++) {
    if (it->first.Z() > 0 && it->first.A() > 0 && it->first.I() >= 0) {
      Nuc MyNuc =
          it->first.Z() * 10000000 + it->first.A() * 10000 + it->first.I();
      myCompMap.insert(pair<Nuc, double>(MyNuc, it->second));
    }
  }
  return myCompMap;
}

//________________________________________________________________________

}  // namespace cyclass
