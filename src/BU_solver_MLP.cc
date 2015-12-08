#include "bu_solver_mlp.h"

//#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <sstream>
#include <map>



using cyclus::Material;
using cyclus::Composition;

typedef std::map<cyclus::Nuc, double> CompMap;

namespace cybam {

    //________________________________________________________________________
    MLPBUsolver::MLPBUsolver(std::string inputfile){

        TMVAWeightFile = inputfile;
        reader = new TMVA::Reader( "Silent" );


    }

    //________________________________________________________________________
    MLPBUsolver::~MLPBUsolver(){
        delete reader;

    }

    //________________________________________________________________________
    TTree* MLPBUsolver::CreateTMVAInputTree(cyclus::Composition::Ptr c_fissil,
                                            cyclus::Composition::Ptr c_fertil,
                                            double BurnUp) {
        float Pu8  = 0;
        float Pu9  = 0;
        float Pu10 = 0;
        float Pu11 = 0;
        float Pu12 = 0;
        float Am11  = 0;
        float U5   = 0;
        float U8   = 0;
        float U4   = 0;
        float BU   = BurnUp;


        //Prepare the input TTree with correct Branch
        TTree*   InputTree = new TTree("EQTMP", "EQTMP");
        InputTree->Branch(	"Pu8"	,&Pu8	,"Pu8/F"	);
        InputTree->Branch(	"Pu9"	,&Pu9	,"Pu9/F"	);
        InputTree->Branch(	"Pu10"	,&Pu10	,"Pu10/F"	);
        InputTree->Branch(	"Pu11"	,&Pu11	,"Pu11/F"	);
        InputTree->Branch(	"Pu12"	,&Pu12	,"Pu12/F"	);
        InputTree->Branch(	"Am1"	,&Am11	,"Am1/F"	);
        InputTree->Branch(	"U5_enrichment"	,&U5	,"U5_enrichment/F"	);
        InputTree->Branch(	"BU"	,&BU	,"BU/F"	);



        //Get Pu composition
        CompMap fissil_map = c_fissil->atom();

        cyclus::CompMap::iterator it;
        for (it = fissil_map.begin(); it != fissil_map.end(); it++){
            cyclus::Nuc nuc = it->first;
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
                Am11 = Q;
            } else  {
                std::cout << "Pb exception... BAD Pu stream" << std::endl;
            }
        }

        // Normalize Pu composition
        double Pu = Pu8 + Pu9 + Pu10 + Pu11 + Pu12 + Am11;
        Pu8  *= 1/Pu;
        Pu9  *= 1/Pu;
        Pu10 *= 1/Pu;
        Pu11 *= 1/Pu;
        Pu12 *= 1/Pu;
        Am11  *= 1/Pu;


        //Get U composition
        CompMap fertil_map = c_fertil->atom();

        for (it = fertil_map.begin(); it != fertil_map.end(); it++){
            cyclus::Nuc nuc = it->first;
            double Q = it->second;

            if (nuc == 922380000) {
                U8 = Q;
            } else if (nuc == 922350000) {
                U5 = Q;
            } else if (nuc == 922340000) {
                U4 = Q;
            } else  {
                std::cout << "Pb exception... BAD U stream" << std::endl;
            }
        }

        // Normalize U5 enrchiments
        U5 *= 1/(U8 + U5 + U4);




        //Fill the input TTree
        InputTree->Fill();

        return InputTree;
    }


    //________________________________________________________________________
    double MLPBUsolver::GetEnrichment(cyclus::Composition::Ptr c_fissil,
                                      cyclus::Composition::Ptr c_fertil,
                                      double BurnUp){

        // Create a set of variables and declare them to the reader
        // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
        Float_t Pu8,Pu9,Pu10,Pu11,Pu12,Am1,BU,U5_enrichment;

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

        // Book method MLP
        TTree* theTree = CreateTMVAInputTree(c_fissil, c_fertil, BurnUp);

        theTree->SetBranchAddress( "BU"   			,&BU 	);
        theTree->SetBranchAddress( "U5_enrichment"  ,&U5_enrichment  )	;
        theTree->SetBranchAddress( "Pu8"  			,&Pu8  );
        theTree->SetBranchAddress( "Pu9"  			,&Pu9  );
        theTree->SetBranchAddress( "Pu10" 			,&Pu10 );
        theTree->SetBranchAddress( "Pu11" 			,&Pu11 );
        theTree->SetBranchAddress( "Pu12" 			,&Pu12 );
        theTree->SetBranchAddress( "Am1"  			,&Am1  );
        theTree->GetEntry(0);

        Float_t val = (reader->EvaluateRegression( methodName ))[0];

        delete theTree;

        return (double)val; //
    }


    //________________________________________________________________________
    double MLPBUsolver::GetBU(cyclus::Composition::Ptr fuel, double eps)
    {

    double BU_max = 60;
    double BU_min = 20;

    cyclus::Composition::Ptr fuel_fissil = ExtractAccordinglist( fuel, fissil_list);
    cyclus::Composition::Ptr fuel_fertil = ExtractAccordinglist( fuel, fertil_list);

    if(AtomIn(fuel_fertil) + AtomIn(fuel_fissil) != AtomIn(fuel)){
        std::cout << "You fuel has nuclei that this model could not manage.."<< std::endl;
        exit(1);
    }


    double rho_target = AtomIn(fuel_fissil)/AtomIn(fuel);

    double rho_min = GetEnrichment(fuel_fissil, fuel_fertil, BU_min);
    double rho_max = GetEnrichment(fuel_fissil, fuel_fertil, BU_max);
    double BU_estimation = 0;
    double rho_estimated = 0;


    do {
        //Update BU_estimation
        BU_estimation = (BU_max+BU_min)/2;

        rho_estimated = GetEnrichment(fuel_fissil, fuel_fertil, BU_estimation);

        if(rho_estimated == rho_target){

            return BU_estimation;

        } else if (rho_estimated > rho_target){

            rho_min = rho_estimated;
            BU_min = BU_estimation;

        } else {

            rho_max = rho_estimated;
            BU_max = BU_estimation;

        }

    }while( std::abs(rho_target - rho_estimated) > eps );

    return BU_estimation;

    }


    //________________________________________________________________________
    double AtomIn(cyclus::Composition::Ptr Source){

        double total = 0;
        CompMap Source_map = Source->atom();

        CompMap::iterator it;

        for(it = Source_map.begin(); it != Source_map.end(); it++)
            total += it->second;

        return total;
    }

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
                separatedCompo.insert( std::pair<cyclus::Nuc,double>(it->first, it2->second) );
            
        }
        
        return Composition::CreateFromAtom(separatedCompo);
    }
    
    
    
    
}
