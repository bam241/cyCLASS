#include "dummy.h"

#include <iostream>


#define SHOW(X)                                                     \
std::cout << std::setprecision(17) << __FILE__ << ":" << __LINE__ \
<< ": " #X " = " << X << "\n"

namespace cybam {

    dummy::dummy() {

        std::cout << " toto" << std::endl;
    }

}  // namespace cybam
