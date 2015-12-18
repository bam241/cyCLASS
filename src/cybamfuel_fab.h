#ifndef CYBAM_SRC_FUEL_FAB_H_
#define CYBAM_SRC_FUEL_FAB_H_



#include "bu_solver_mlp.h"

#include <string>
#include "cyclus.h"

namespace cybam {

    /// cybamFuelFab takes in 2 streams of material and mixes them in ratios in order to
    /// supply material that matches some neutronics properties of reqeusted
    /// material.  
    /// @endcode
    class cybamFuelFab : public cyclus::Facility {
#pragma cyclus note { \
"niche": "fabrication", \
"doc": \
"cybamFuelFab takes in 2 streams of material and mixes them in ratios in order to" \
"", \
}
    public:
        cybamFuelFab(cyclus::Context* ctx);
        virtual ~cybamFuelFab(){};


#pragma cyclus

        virtual void Tick(){};
        virtual void Tock(){};
        virtual void EnterNotify();

        virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(
                                                                                  cyclus::CommodMap<cyclus::Material>::type& commod_requests);

        virtual void GetMatlTrades(
                                   const std::vector<cyclus::Trade<cyclus::Material> >& trades,
                                   std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                   cyclus::Material::Ptr> >& responses);

        virtual void AcceptMatlTrades(const std::vector<std::pair<
                                      cyclus::Trade<cyclus::Material>, cyclus::Material::Ptr> >& responses);

        virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
        GetMatlRequests();

    private:

        MLPBUsolver* MyBUSolver;


#pragma cyclus var { \
"doc": "Ordered list of commodities on which to requesting filler stream material.", \
"uilabel": "Filler Stream Commodities", \
"uitype": ["oneormore", "incommodity"], \
}
        std::vector<std::string> fill_commods;
#pragma cyclus var { \
"default": [], \
"uilabel": "Filler Stream Preferences", \
"doc": "Filler stream commodity request preferences for each of the given filler commodities (same order)." \
" If unspecified, default is to use 1.0 for all preferences.", \
}
        std::vector<double> fill_commod_prefs;
#pragma cyclus var { \
"doc": "Name of recipe to be used in filler material stream requests.", \
"uilabel": "Filler Stream Recipe", \
"uitype": "recipe", \
}
        std::string fill_recipe;
#pragma cyclus var { \
"doc": "Size of filler material stream inventory.", \
"uilabel": "Filler Stream Inventory Capacity", \
"units": "kg", \
}
        double fill_size;
#pragma cyclus var {"capacity": "fill_size"}
        cyclus::toolkit::ResBuf<cyclus::Material> fill;

#pragma cyclus var { \
"doc": "Ordered list of commodities on which to requesting fissile stream material.", \
"uilabel": "Fissile Stream Commodities", \
"uitype": ["oneormore", "incommodity"], \
}
        std::vector<std::string> fiss_commods;
#pragma cyclus var { \
"default": [], \
"uilabel": "Fissile Stream Preferences", \
"doc": "Fissile stream commodity request preferences for each of the given fissile commodities (same order)." \
" If unspecified, default is to use 1.0 for all preferences.", \
}
        std::vector<double> fiss_commod_prefs;
#pragma cyclus var { \
"doc": "Name for recipe to be used in fissile stream requests." \
" Empty string results in use of an empty dummy recipe.", \
"uitype": "recipe", \
"uilabel": "Fissile Stream Recipe", \
"default": "", \
}
        std::string fiss_recipe;
#pragma cyclus var { \
"doc": "Size of fissile material stream inventory.", \
"uilabel": "Fissile Stream Inventory Capacity", \
"units": "kg", \
}
        double fiss_size;
#pragma cyclus var {"capacity": "fiss_size"}
        cyclus::toolkit::ResBuf<cyclus::Material> fiss;

#pragma cyclus var { \
"doc": "Commodity on which to request material for top-up stream." \
" This MUST be set if 'topup_size > 0'.", \
"uilabel": "Top-up Stream Commodity", \
"default": "", \
"uitype": "incommodity", \
}

#pragma cyclus var { \
"doc": "Commodity on which to offer/supply mixed fuel material.", \
"uilabel": "Output Commodity", \
"uitype": "outcommodity", \
}
        std::string outcommod;

#pragma cyclus var { \
"doc": "Maximum number of kg of fuel material that can be supplied per time step.", \
"uilabel": "Maximum Throughput", \
"units": "kg", \
}
        double throughput;


        // intra-time-step state - no need to be a state var
        // map<request, inventory name>
        std::map<cyclus::Request<cyclus::Material>*, std::string> req_inventories_;
    };

    double AtomToMassFrac(double atomfrac, cyclus::Composition::Ptr c1, cyclus::Composition::Ptr c2);
    
} // namespace cybam


#endif  // CYBAM_SRC_FUEL_FAB_H_
