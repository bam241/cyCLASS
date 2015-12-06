#include "bu_solver_mlp.h"

//#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include <sstream>




using cyclus::Material;
using cyclus::Composition;
using pyne::simple_xs;

typedef std::map<int, double> CompMap;

namespace cycbam {


    //________________________________________________________________________
    TTree* MLPBUsolver::CreateTMVAInputTree(cyclus::Composition::Ptr c_fissil,
                                            cyclus::Composition::Ptr c_fertil,
                                            double BurnUp) {
        float Pu8  = 0;
        float Pu9  = 0;
        float Pu10 = 0;
        float Pu11 = 0;
        float Pu12 = 0;
        float Am1  = 0;
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
        InputTree->Branch(	"Am1"	,&Am1	,"Am1/F"	);
        InputTree->Branch(	"U5_enrichment"	,&U5	,"U5_enrichment/F"	);
        InputTree->Branch(	"BU"	,&BU	,"BU/F"	);



        //Get Pu composition
        CompMap fissil_map = f_fissil.atom();

        cyclus::CompMap::iterator it;
        for (it = f_fissil.begin(); it != f_fissil.end(); it++){
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
                Pu11 = Q;
            } else if (nuc == 952410000) {
                Am11 = Q;
            } else  {
                cout << "Pb exception... BAD Pu stream" << endl;
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
        CompMap fertil_map = f_fissil.atom();

        for (it = f_fissil.begin(); it != f_fissil.end(); it++){
            cyclus::Nuc nuc = it->first;
            double Q = it->second;

            if (nuc == 922380000) {
                U8 = Q;
            } else if (nuc == 922350000) {
                U5 = Q;
            } else if (nuc == 922340000) {
                U4 = Q;
            } else  {
                cout << "Pb exception... BAD U stream" << endl;
            }
        }

        // Normalize U5 enrchiments
        U5 *= 1/(U8 + U5 + U4);




        //Fill the input TTree
        InputTree->Fill();

        return InputTree;
    }


    //________________________________________________________________________
    double EQM_PWR_MLP_MOX::ExecuteTMVA(TTree* theTree)
    {
    // --- Create the Reader object
    TMVA::Reader *reader = new TMVA::Reader( "Silent" );
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

    // Book method MLP
    TString methodName = "MLP method";
    reader->BookMVA( methodName, TMVAWeightFile );
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
    
    delete reader;
    delete theTree;
    
    return (double)val; //retourne teneur
    }
    
}
