#ifndef CYBAM_SRC_FUEL_FAB_H_
#define CYBAM_SRC_FUEL_FAB_H_

#include "TMVA/Reader.h"
#include "cyclus.h"

class TTree;
typedef std::map<cyclus::Nuc, double> CompMap;


/// For the moment, this is just an adaptation of the CLASS-MLP Model to CYCLUS
/// @code
/// [1] B. Leniau, B. Mouginota, N. Thiolliere et al. "A neural network approach
/// for burn-up calculation and its application to the dynamic fuel cycle code
/// CLASS" Annals of Nuclear Energy, Volume 81, July 2015
/// @endcode
namespace cybam {

    class MLPBUsolver {
    public:
        MLPBUsolver(std::string inputfile = "/Users/mouginot/scratch/App/CLASS/Database/PWR/MOX/EQModel/EQM_MLP_PWR_MOX_3batch.xml");
        ~MLPBUsolver();


        // Create a TMWA input with the compistion of both stream and the targeted BU..
        TTree* CreateTMVAInputTree(cyclus::Composition::Ptr c_fissil,
                                   cyclus::Composition::Ptr c_fertil,
                                   double BurnUp);

        // Lauch the MLP (using TMVA) to predict the requeirt fissil enrichment (according to the compisition of fissil and fertil stream and the targeted burnup)
        double GetEnrichment(cyclus::Composition::Ptr c_fissil,
                           cyclus::Composition::Ptr c_fertil,
                           double BurnUp);

        // liner dicchotomy to determine the Burn-up reachable by a fuel depending of its composition...
        double GetBU(cyclus::Composition::Ptr fuel, double eps = 1e-6 );

        std::string str() { return TMVAWeightFile;};

    private:

        // Fissil and fertil list should be provided with the tmva input xml file
        cyclus::Composition::Ptr fissil_list; // list of nuclei composing the fissil stream
        cyclus::Composition::Ptr fertil_list; // list of nuclei composing the fertil stream

        TMVA::Reader *reader;
        std::string TMVAWeightFile;


    };

    double AtomIn(CompMap Source);
    double AtomIn(cyclus::Composition::Ptr Source) { return AtomIn(Source->atom());};
    cyclus::Composition::Ptr ExtractAccordinglist( cyclus::Composition::Ptr source, cyclus::Composition::Ptr list);
    CompMap NormalizeComp( CompMap source, double norm = 1);

} // namespace cybam

#endif  // cybam_SRC_FUEL_FAB_H_
