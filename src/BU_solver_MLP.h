#ifndef CYCBAM_SRC_FUEL_FAB_H_
#define CYCBAM_SRC_FUEL_FAB_H_

#include <string>
#include "cyclus.h"

namespace cycbam {

        class MLPBUsolver : {
#pragma cyclus note { \
}
    public:
        MLPBUsolver();
        virtual ~MLPBUsolver(){};

        virtual std::string version() { return CYCBAM_VERSION; }


            TTree* MLPBUsolver::CreateTMVAInputTree(cyclus::Composition::Ptr c_fissil,
                                                    cyclus::Composition::Ptr c_fertil,
                                                    double BurnUp);
            double EQM_PWR_MLP_MOX::ExecuteTMVA(TTree* theTree)


#pragma cyclus


    private:
#pragma cyclus var { \
"doc": "file containing the TMVA weight fot the MLP" \
"uilabel": "TMVA weight file name", \
"default": "", \
"uitype": "", \
}
            std::string TMVAWeightFile;

} // namespace cycbam


#endif  // cycbam_SRC_FUEL_FAB_H_
