#include "CLASSAdaptator.h"

#include "EvolutionData.hxx"

#include "Irradiation/IM_Matrix.hxx"
#include "Irradiation/IM_RK4.hxx"

#include "XS/XSM_MLP.hxx"


#include "Equivalence/EQM_FBR_BakerRoss_MOX.hxx"  
#include "Equivalence/EQM_FBR_MLP_Keff_BOUND.hxx"
#include "Equivalence/EQM_PWR_LIN_MOX.hxx"
#include "Equivalence/EQM_PWR_MLP_MOX_Am.hxx"
#include "Equivalence/EQM_PWR_QUAD_MOX.hxx"
#include "Equivalence/EQM_FBR_MLP_Keff.hxx"
#include "Equivalence/EQM_MLP_Kinf.hxx"
#include "Equivalence/EQM_PWR_MLP_MOX.hxx"
#include "Equivalence/EQM_PWR_POL_UO2.hxx"

#include <string>

using cyclus::Nuc;
using cyclus::Material;
using cyclus::Composition;


//#define cyDBGL		std::cout << __FILE__ << " : " << __LINE__ << " [" << __FUNCTION__ << "]" << std::endl;
#define cyDBGL ;


namespace cyclass {



    EquivalenceModel* EQmodelfor(std::string name, std::stringstream& command ) {

        if( name == "FBR_BakerRoss_MOX" ){

            std::string buff;
            getline( command, buff, ',' );
            double w_235u = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_238pu = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_240pu = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_241pu = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_242pu = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_241am = atof(buff.c_str());
            getline( command, buff, ',' );
            double w_eqfissil = atof(buff.c_str());


            return new EQM_FBR_BakerRoss_MOX(w_235u, w_238pu, w_240pu, w_241pu, w_242pu, w_241am, w_eqfissil );

        }else if( name == "FBR_MLP_Keff" ){

            std::string TMVAWeightPath;
            getline( command, TMVAWeightPath, ',' );

            std::string buff;
            getline( command, buff, ',' );
            double keff_target = atof(buff.c_str());

            std::string InformationFile;
            getline( command, InformationFile, ',' );

            return new EQM_FBR_MLP_Keff(TMVAWeightPath, keff_target, InformationFile);


        }else if( name == "FBR_MLP_Keff_BOUND" ){

            std::string TMVAWeightPath;
            getline( command, TMVAWeightPath, ',' );

            std::string buff;
            getline( command, buff, ',' );
            int NumOfBatch = atoi(buff.c_str());

            getline( command, buff, ',' );
            double LowerKeffective = atof(buff.c_str());

            getline( command, buff, ',' );
            double UpperKeffective = atof(buff.c_str());

            std::string InformationFile;
            getline( command, InformationFile, ',' );



            return new EQM_FBR_MLP_Keff_BOUND(TMVAWeightPath,
                                              NumOfBatch,
                                              LowerKeffective,
                                              UpperKeffective,
                                              InformationFile);


        }else if( name == "MLP_Kinf" ){


            std::string WeightPathAlpha0;
            getline( command, WeightPathAlpha0, ',' );

            std::string WeightPathAlpha1;
            getline( command, WeightPathAlpha1, ',' );

            std::string WeightPathAlpha2;
            getline( command, WeightPathAlpha2, ',' );

            std::string InformationFile;
            getline( command, InformationFile, ',' );

            std::string buff;
            getline( command, buff, ',' );
            int NumOfBatch = atoi(buff.c_str());

            getline( command, buff, ',' );
            double CriticalityThreshold = atof(buff.c_str());


            return new EQM_MLP_Kinf(WeightPathAlpha0,
                                    WeightPathAlpha1,
                                    WeightPathAlpha2,
                                    InformationFile,
                                    NumOfBatch,
                                    CriticalityThreshold);

        }else if( name == "PWR_LIN_MOX" ){

            return new EQM_PWR_LIN_MOX(command.str());

        }else if( name == "PWR_MLP_MOX_Am" ){

            return new EQM_PWR_MLP_MOX_AM(command.str());
        }else if( name == "PWR_QUAD_MOX" ){

            return new EQM_PWR_QUAD_MOX(command.str());
        }else if( name == "PWR_MLP_MOX" ){

            return new EQM_PWR_MLP_MOX(command.str());
        }else if( name == "PWR_POL_UO2" ){

            return new EQM_PWR_POL_UO2(command.str());
        }else {
            throw cyclus::ValidationError("Bad EqModel keyword");
        }
    }

    XSModel* XSmodelfor(std::string name, std::stringstream& command){

        if( name == "XSM_MLP" ){
            return new XSM_MLP(command.str());
        }else {
            throw cyclus::ValidationError("Bad XSModel keyword");
        }
    }

    IrradiationModel* IMmodelfor(std::string name, std::stringstream& command){

        if( name == "RK4" ){
            return new IM_RK4();
        }else if( name == "MATRIX" ){
            return new IM_Matrix();
        }else{
            throw cyclus::ValidationError("Bad IMModel keyword");
        }
        
    }
    CLASSAdaptator::CLASSAdaptator(std::string EQModel, std::stringstream& EQcommand,
                                   std::string XSModel, std::stringstream& XScommand,
                                   std::string IMModel, std::stringstream& IMcommand){

        cyDBGL
        myPhysicsModel = new PhysicsModels();

        myPhysicsModel->SetEquivlalenceModel(   EQmodelfor(EQModel, EQcommand) );
        myPhysicsModel->SetXSModel(             XSmodelfor(XSModel, XScommand) );
        myPhysicsModel->SetIrradiationModel(    IMmodelfor(IMModel, IMcommand) );

        cyDBGL

    }

    //________________________________________________________________________
    float CLASSAdaptator::GetEnrichment(cyclus::Composition::Ptr c_fissil,
                                        cyclus::Composition::Ptr c_fertil,
                                        double BurnUp) const{
        cyDBGL

        double val = 0;

        IsotopicVector IV_fissil = CYCLUS2CLASS(c_fissil);
        IsotopicVector IV_fertil = CYCLUS2CLASS(c_fertil);

        val= myPhysicsModel->GetEquivalenceModel()->GetFissileMolarFraction(IV_fissil, IV_fertil, BurnUp);

        cyDBGL

        return val; //
    }


    //________________________________________________________________________
    float CLASSAdaptator::GetBU(cyclus::Composition::Ptr fuel, double eps) const{
        cyDBGL

        float BU_max = 80;
        float BU_min = 10;
        cyDBGL

        cyclus::Composition::Ptr fuel_fissil = ExtractAccordinglist( fuel, fissil_list);
        cyclus::Composition::Ptr fuel_fertil = ExtractAccordinglist( fuel, fertil_list);



        if( std::abs(AtomIn(fuel_fertil) + AtomIn(fuel_fissil) - AtomIn(fuel)) > 1e-10 ){
            std::cout << "You fuel has nuclei that this model could not manage.."<< std::endl;
            exit(1);
        }


        float rho_target = AtomIn(fuel_fissil)/AtomIn(fuel);
        float rho_min = GetEnrichment(fuel_fissil, fuel_fertil, BU_min);
        float rho_max = GetEnrichment(fuel_fissil, fuel_fertil, BU_max);
        float BU_estimation = BU_max;
        float rho_estimated = rho_max;

        do {
            //Update BU_estimation

            if (rho_estimated > rho_target){
                rho_max = rho_estimated;
                BU_max = BU_estimation;

            } else {
                rho_min = rho_estimated;
                BU_min = BU_estimation;
            }

            BU_estimation = (BU_max+BU_min)/2.;

            rho_estimated = GetEnrichment(fuel_fissil, fuel_fertil, BU_estimation);


        }while( std::abs(rho_target - rho_estimated)/(rho_target/2.+rho_estimated/2.) > eps );

        cyDBGL
        return BU_estimation;


    }

    cyclus::Composition::Ptr CLASSAdaptator::GetCompAfterIrradiation(cyclus::Composition::Ptr InitialCompo, double power, double mass, double burnup){

        IsotopicVector InitialIV = CYCLUS2CLASS(InitialCompo);

        cSecond finaltime = burnup /(power*1e-3) /mass *3600*24;

        EvolutionData myEvolution = myPhysicsModel->GenerateEvolutionData(InitialIV, finaltime, power);

        IsotopicVector AfterIrradiationIV = myEvolution.GetIsotopicVectorAt(finaltime);

        return Composition::CreateFromAtom(CLASS2CYCLUS(AfterIrradiationIV));

    }


    //________________________________________________________________________
    //________________________________________________________________________
    //________________________________________________________________________
    double AtomIn(CompMap Source){

        double total = 0;

        CompMap::iterator it;

        for(it = Source.begin(); it != Source.end(); it++)
            total += it->second;

        return total;
    }

    //________________________________________________________________________
    double AtomIn(cyclus::Composition::Ptr Source) { return AtomIn(Source->atom());};

    //________________________________________________________________________
    cyclus::Composition::Ptr ExtractAccordinglist(
                                                  cyclus::Composition::Ptr source,
                                                  cyclus::Composition::Ptr list){

        //create Output Composition
        CompMap separatedCompo;

        // Extract Nuc map from source compo & list...
        CompMap sourceComp = source->atom();
        CompMap ListComp = list->atom();

        CompMap::iterator it;

        // Fill output composition
        for (it = ListComp.begin(); it != ListComp.end(); it++) {
            CompMap::iterator it2 = sourceComp.find( it->first );

            if(it2 != sourceComp.end())
                separatedCompo.insert( std::pair<Nuc,double>(it->first, it2->second) );

        }

        return Composition::CreateFromAtom(separatedCompo);
    }

    //________________________________________________________________________
    CompMap NormalizeComp( CompMap source, double norm ){
        double Total = AtomIn(source);
        CompMap::iterator it;
        for(it = source.begin(); it != source.end(); it++)
            it->second *= norm/Total;

        return source;
    }

    //________________________________________________________________________
    CompMap operator+( CompMap const& IVa, CompMap const& IVb){

        CompMap IVtmp = IVa;
        CompMap IVbtmp = IVb;
        CompMap::iterator it;

        for(it = IVbtmp.begin(); it != IVbtmp.end(); it++){

            std::pair<CompMap::iterator, bool> IResult;
            IResult = IVtmp.insert(std::pair<Nuc, double> (it->first, it->second));
            if(!IResult.second)
                IResult.first->second += it->second;

        }
        return IVtmp;
    }


    //________________________________________________________________________
    CompMap operator*(CompMap const& IVA, double F){

        CompMap IVtmp = IVA;
        CompMap::iterator it;

        for(it = IVtmp.begin(); it != IVtmp.end(); it++)
            it->second *= F;

        return IVtmp;
    }

    //________________________________________________________________________
    void Print(cyclus::Composition::Ptr compo ){
        Print(compo->atom());
    }

    //________________________________________________________________________
    void Print(CompMap compo){
        CompMap::iterator it;

        for (it = compo.begin(); it != compo.end(); it++){
            std::cout << it->first << " " << it->second << std::endl;
        }
    }

    //________________________________________________________________________
    CompMap operator-(CompMap const& IVa, CompMap const& IVb) { return IVa + (-1*IVb); };

    //________________________________________________________________________
    CompMap operator/(CompMap const& IVA, double F) { return IVA * (1/F); };

    //________________________________________________________________________
    CompMap operator*(double F, CompMap const& IVA) { return IVA*F; };

    //________________________________________________________________________
    IsotopicVector CYCLUS2CLASS(CompMap c_compo){

        IsotopicVector IV;
        map<int, double>::iterator it;

        for( it = c_compo.begin(); it != c_compo.end(); it++){
            int Z = (int)it->first/10000 ;
            int A = (int)(it->first/10)%1000;
            int I = (int) it->first%10%10;

            IV += it->second * ZAI(Z, A, I);
        }

        return IV;
    }

    //________________________________________________________________________
    IsotopicVector CYCLUS2CLASS(cyclus::Composition::Ptr c_compo){
        return CYCLUS2CLASS(c_compo->atom());

    }

    //________________________________________________________________________
    CompMap CLASS2CYCLUS(IsotopicVector const& IV){

        CompMap myCompMap;

        map<ZAI, double> IsotopicMap = IV.GetIsotopicQuantity();
        map<ZAI, double>::iterator it;

        for(it = IsotopicMap.begin(); it != IsotopicMap.end(); it++){

            Nuc MyNuc = it->first.Z()*10000 + it->first.A()*10 + it->first.I();
            myCompMap.insert( pair<Nuc, double> (MyNuc, it->second) );
        }
        return myCompMap;
    }


    //________________________________________________________________________



} // namespace cyclass
