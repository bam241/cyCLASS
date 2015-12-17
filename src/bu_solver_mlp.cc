#include "bu_solver_mlp.h"

#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"


#include <sstream>
#include <map>
#include <string>

using cyclus::Nuc;
using cyclus::Material;
using cyclus::Composition;

//#define DBGL		std::cout << __FILE__ << " : " << __LINE__ << " [" << __FUNCTION__ << "]" << std::endl;
#define DBGL		;

namespace cybam {

    //________________________________________________________________________
    MLPBUsolver::MLPBUsolver(std::string inputfile){
        DBGL
        TMVAWeightFile = inputfile;
        reader = new TMVA::Reader( "Silent" );

        reader->AddVariable( "BU"   		,&BU );
        reader->AddVariable( "U5_enrichment",&U5_enrichment );
        reader->AddVariable( "Pu8"  		,&Pu8 );
        reader->AddVariable( "Pu9"  		,&Pu9 );
        reader->AddVariable( "Pu10" 		,&Pu10);
        reader->AddVariable( "Pu11" 		,&Pu11);
        reader->AddVariable( "Pu12" 		,&Pu12);
        reader->AddVariable( "Am1"  		,&Am1 );

        // --- Book the MVA methods
        TString methodName = "MLP method";
        reader->BookMVA( methodName, TMVAWeightFile );


        CompMap fissil;
        fissil.insert(std::pair<Nuc, double>(942380000,1));
        fissil.insert(std::pair<Nuc, double>(942390000,1));
        fissil.insert(std::pair<Nuc, double>(942400000,1));
        fissil.insert(std::pair<Nuc, double>(942410000,1));
        fissil.insert(std::pair<Nuc, double>(942420000,1));
        fissil.insert(std::pair<Nuc, double>(952410000,1));
        fissil_list = Composition::CreateFromAtom(fissil);

        CompMap fertil;
        fertil.insert(std::pair<Nuc, double>(922350000,1));
        fertil.insert(std::pair<Nuc, double>(922380000,1));
        fertil_list = Composition::CreateFromAtom(fertil);
        DBGL
    }

    //________________________________________________________________________
    MLPBUsolver::~MLPBUsolver(){
        delete reader;
    }

    //________________________________________________________________________
    void MLPBUsolver::UpdateInputComposition(cyclus::Composition::Ptr c_fissil,
                                             cyclus::Composition::Ptr c_fertil,
                                             double BurnUp) {
        DBGL;

        Pu8  = 0;
        Pu9  = 0;
        Pu10 = 0;
        Pu11 = 0;
        Pu12 = 0;
        Am1 = 0;
        U5_enrichment   = 0;
        BU   = BurnUp;

        float U8   = 0;
        float U4   = 0;


        //Get Pu composition
        CompMap fissil_map = c_fissil->atom();

        CompMap::iterator it;
        for (it = fissil_map.begin(); it != fissil_map.end(); it++){
            Nuc nuc = it->first;
            double Q = it->second;

            if (nuc == 942380000) {
                Pu8 = Q;
            } else if (nuc == 942390000) {
                Pu9 = Q;
            } else if (nuc == 942400000) {
                Pu10 = Q;
            } else if (nuc == 942410000) {
                Pu11 = Q;
            } else if (nuc == 942420000) {
                Pu12 = Q;
            } else if (nuc == 952410000) {
                Am1 = Q;
            } else  {
                std::cout << "Pb exception... BAD Pu stream nuclei " << nuc << std::endl;
            }
        }

        // Normalize Pu composition
        double Pu = Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am1;
        Pu8  *= 1/Pu;
        Pu9  *= 1/Pu;
        Pu10 *= 1/Pu;
        Pu11 *= 1/Pu;
        Pu12 *= 1/Pu;
        Am1  *= 1/Pu;


        //Get U composition
        CompMap fertil_map = c_fertil->atom();

        for (it = fertil_map.begin(); it != fertil_map.end(); it++){
            Nuc nuc = it->first;
            double Q = it->second;

            if (nuc == 922380000) {
                U8 = Q;
            } else if (nuc == 922350000) {
                U5_enrichment = Q;
            } else if (nuc == 922340000) {
                U4 = Q;
            } else  {
                std::cout << "Pb exception... BAD U stream" << std::endl;
            }
        }

        // Normalize U5 enrchiments
        U5_enrichment *= 1/(U8 + U5_enrichment + U4);

        DBGL
    }


    //________________________________________________________________________
    double MLPBUsolver::GetEnrichment(cyclus::Composition::Ptr c_fissil,
                                      cyclus::Composition::Ptr c_fertil,
                                      double BurnUp){
        DBGL

        UpdateInputComposition(c_fissil, c_fertil, BurnUp);

        Float_t val = (reader->EvaluateRegression(  "MLP method" ))[0];

        DBGL

        return (double)val; //
    }


    //________________________________________________________________________
    double MLPBUsolver::GetBU(cyclus::Composition::Ptr fuel, double eps)
    {
    DBGL

    double BU_max = 60;
    double BU_min = 20;
    DBGL


    cyclus::Composition::Ptr fuel_fissil = ExtractAccordinglist( fuel, fissil_list);
    cyclus::Composition::Ptr fuel_fertil = ExtractAccordinglist( fuel, fertil_list);


    
    if( std::abs(AtomIn(fuel_fertil) + AtomIn(fuel_fissil) - AtomIn(fuel)) > 1e-6 ){
        std::cout << "You fuel has nuclei that this model could not manage.."<< std::endl;
        exit(1);
    }
    DBGL


    double rho_target = AtomIn(fuel_fissil)/AtomIn(fuel);

    DBGL

    double rho_min = GetEnrichment(fuel_fissil, fuel_fertil, BU_min);
    double rho_max = GetEnrichment(fuel_fissil, fuel_fertil, BU_max);
    double BU_estimation = 0;
    double rho_estimated = 0;
    DBGL


    do {
        //Update BU_estimation
        BU_estimation = (BU_max+BU_min)/2;

        rho_estimated = GetEnrichment(fuel_fissil, fuel_fertil, BU_estimation);

        if( rho_estimated == rho_target ){

            return BU_estimation;

        } else if (rho_estimated > rho_target){
            rho_max = rho_estimated;
            BU_max = BU_estimation;

        } else {
            rho_min = rho_estimated;
            BU_min = BU_estimation;


        }

    }while( std::abs(rho_target - rho_estimated)/rho_target > eps );

    DBGL
    return BU_estimation;

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
    double AtomIn(cyclus::Composition::Ptr Source) { return AtomIn(Source->atom());};

    //________________________________________________________________________
    cyclus::Composition::Ptr ExtractAccordinglist( cyclus::Composition::Ptr source, cyclus::Composition::Ptr list){

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

CompMap operator-(CompMap const& IVa, CompMap const& IVb) { return IVa + (-1*IVb); };
CompMap operator/(CompMap const& IVA, double F) { return IVA * (1/F); };
CompMap operator*(double F, CompMap const& IVA) { return IVA*F; };



//________________________________________________________________________





} // namespace cybam
