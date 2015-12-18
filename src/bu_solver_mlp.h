#ifndef CYBAM_SRC_MLP_BU_H_
#define CYBAM_SRC_MLP_BU_H_

#include "cyclus.h"
#include "TMVA/Reader.h"


class TTree;

using cyclus::Nuc;
typedef std::map<Nuc, double> CompMap;


/// For the moment, this is just an adaptation of the CLASS-MLP Model to CYCLUS
/// @code
/// [1] B. Leniau, B. Mouginot, N. Thiolliere et al. "A neural network approach
/// for burn-up calculation and its application to the dynamic fuel cycle code
/// CLASS" Annals of Nuclear Energy, Volume 81, July 2015
/// @endcode




namespace cybam {

    class MLPBUsolver {
    public:
        MLPBUsolver(std::string inputfile = "/Users/mouginot/scratch/App/CLASS/Database/PWR/MOX/EQModel/EQM_MLP_PWR_MOX_3batch.xml");
        ~MLPBUsolver();



        // Lauch the MLP (using TMVA) to predict the requeirt fissil enrichment (according to the compisition of fissil and fertil stream and the targeted burnup)
        float GetEnrichment(cyclus::Composition::Ptr c_fissil,
                           cyclus::Composition::Ptr c_fertil,
                           double BurnUp);

        // liner dicchotomy to determine the Burn-up reachable by a fuel depending of its composition...
        float GetBU(cyclus::Composition::Ptr fuel, double eps = 1e-6 );

        std::string str() { return TMVAWeightFile;};

    private:

        // Fissil and fertil list should be provided with the tmva input xml file
        cyclus::Composition::Ptr fissil_list; // list of nuclei composing the fissil stream
        cyclus::Composition::Ptr fertil_list; // list of nuclei composing the fertil stream

        float Pu8;
        float Pu9;
        float Pu10;
        float Pu11;
        float Pu12;
        float Am1;
        float BU;
        float U5_enrichment;

        TMVA::Reader *reader;
        std::string TMVAWeightFile;

        // Create a TMWA input with the compistion of both stream and the targeted BU..
        void UpdateInputComposition(cyclus::Composition::Ptr c_fissil,
                                    cyclus::Composition::Ptr c_fertil,
                                    double BurnUp);


    };

    double AtomIn(CompMap Source);
    double AtomIn(cyclus::Composition::Ptr Source);
    cyclus::Composition::Ptr ExtractAccordinglist( cyclus::Composition::Ptr source, cyclus::Composition::Ptr list);
    CompMap NormalizeComp( CompMap source, double norm = 1);
    void Print(cyclus::Composition::Ptr compo);
    void Print(CompMap compo);



    CompMap operator+(CompMap const& IVa, CompMap const& IVb);
    CompMap operator*(CompMap const& IVA, double F);
    CompMap operator*(double F, CompMap const& IVA);

    CompMap operator-(CompMap const& IVa, CompMap const& IVb);
    CompMap operator/(CompMap const& IVA, double F);



} // namespace cybam

#endif  // cybam_SRC_FUEL_FAB_H_
