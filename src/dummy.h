#ifndef CYBAM_SRC_FUEL_FAB_H_
#define CYBAM_SRC_FUEL_FAB_H_



#include "bu_solver_mlp.h"

namespace cybam {

     class dummy  {
    public:
        dummy();
        virtual ~dummy(){};



    private:

        MLPBUsolver* MyBUSolver;

    };
    
} // namespace cybam


#endif  // CYBAM_SRC_FUEL_FAB_H_
