#ifndef CYCLASS_SRC_REACTOR_H_
#define CYCLASS_SRC_REACTOR_H_

#include "CLASSAdaptator.h"
#include "cyclus.h"

namespace cyclass {
using cyclus::toolkit::MatVec;
using cyclus::toolkit::ResBuf;

class Reactor : public cyclus::Facility,
                public cyclus::toolkit::CommodityProducer {
  #pragma cyclus note { "niche" : "bamreactor", "doc" : "", }

 public:
  Reactor(cyclus::Context* ctx);
  virtual ~Reactor(){};

  virtual void Tick();
  virtual void Tock();
  virtual void EnterNotify();

  virtual void AcceptMatlTrades(
      const std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                                  cyclus::Material::Ptr>>& responses);

  virtual std::set<cyclus::RequestPortfolio<cyclus::Material>::Ptr>
  GetMatlRequests();

  virtual std::set<cyclus::BidPortfolio<cyclus::Material>::Ptr> GetMatlBids(
      cyclus::CommodMap<cyclus::Material>::type& commod_requests);

  virtual void GetMatlTrades(
      const std::vector<cyclus::Trade<cyclus::Material>>& trades,
      std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                            cyclus::Material::Ptr>>& responses);

  #pragma cyclus decl

 private:
  std::string fuel_incommod(cyclus::Material::Ptr m);
  std::string fuel_outcommod(cyclus::Material::Ptr m);
  double fuel_pref(cyclus::Material::Ptr m);
  bool FullCore();
  bool Discharged();
  bool Refueling();
  bool InCycle();

  CLASSAdaptator* MyCLASSAdaptator;

  bool retired() {
    return exit_time() != -1 && context()->time() >= exit_time();
  }

  /// Store fuel info index for the given resource received on incommod.
  void index_res(cyclus::Resource::Ptr m, std::string incommod);

  /// Discharge a batch from the core if there is room in the spent fuel
  /// inventory.  Returns true if a batch was successfully discharged.
  bool Discharge(int i);

  /// Top up core inventory as much as possible.
  void Load(int i);

  /// Transmute the specified number of assemblies in the core to their
  /// fully burnt state as defined by their outrecipe.
  void Transmute(int n_assem);

  /// Records a reactor event to the output db with the given name and note val.
  void Record(std::string name, std::string val);

  /// Complement of PopSpent - must be called with all materials passed that
  /// were not traded away to other agents.
  void PushSpent(std::map<std::string, MatVec> leftover);

  /// Returns all spent assemblies indexed by outcommod - removing them from
  /// the spent fuel buffer.
  std::map<std::string, MatVec> PopSpent();

  /// Returns all spent assemblies indexed by outcommod without removing them
  /// from the spent fuel buffer.
  std::map<std::string, MatVec> PeekSpent();

/////// fuel specifications /////////
  #pragma cyclus var {                                                     \
    "uitype" :                                                             \
    ["oneormore", "incommodity"], "uilabel" : "Fresh Fuel Commodity List", \
     "doc" : "Ordered list of input commodities on which to requesting fuel.", }
  std::vector<std::string> fuel_incommods;

  #pragma cyclus var { \
    "default": [], \
    "uilabel": "Fresh Fuel Preference List", \
    "doc": "The preference for each type of fresh fuel requested corresponding"\
    " to each input commodity (same order).  If no preferences are " \
    "specified, 1.0 is used for all fuel " \
    "requests (default).", \
  }
  std::vector<double> fuel_prefs;

  #pragma cyclus var { \
    "uitype": ["oneormore", "outcommodity"], \
    "uilabel": "Spent Fuel Commodity List", \
    "doc": "Output commodities on which to offer spent fuel originally " \
    "received as each particular input commodity (same order)." \
  }
  std::vector<std::string> fuel_outcommods;

  #pragma cyclus var {                                                      \
    "doc" : "Mass (kg) of a Batch assembly.", "uilabel" : "Batch Mass", \
    "units" : "kg", }
  double batch_size;

  #pragma cyclus var {                     \
    "uilabel" : "Number of batches in Core", \
    "units" : "batches",                                    \
    "doc" : "Number of batches that constitute a full core.", }
  int n_batch_core;

  #pragma cyclus var {                                         \
    "default" : 0, "uilabel" : "Minimum Fresh Fuel Inventory", \
    "units" : "kg",                                    \
    "doc" : "Mass of fresh fuel for each batch to keep on-hand if possible.", }
  double n_batch_fresh;

  #pragma cyclus var { \
    "default": 1000000000, \
    "uilabel": "Maximum Spent Fuel Inventory", \
    "units": "kg", \
    "doc": "Mass of spent fuel that can be stored on-site before" \
    " reactor operation stalls.", \
  }
  double m_batch_spent;

///////// Model params ///////////

  #pragma cyclus var { \
    "doc" : "Cross Secction Model name.", "uilabel" : "Name of the XS Model", }
    std::string xs_model;

  #pragma cyclus var {                      \
    "doc" : "Cross Section Parameter line", \
    "uilabel" : "Parameter line for the cross section model", }
  std::string xs_command;

  #pragma cyclus var {                                     \
    "default" : "mox", "doc" : "fuel type Parameter line", \
    "uilabel" : "Parameter to define the fuel type", }
  std::string fuel_type;

  #pragma cyclus var {                \
    "doc" : "Irradition Model name.", \
    "uilabel" : "Name of the irradiation Model", }
  std::string ir_model;
  #pragma cyclus var {                   \
    "doc" : "Irradition Parameter line", \
    "uilabel" : "Parameter line for the irradiation model", }
  std::string ir_command;

  #pragma cyclus var {                                   \
    "default" : "eq", "doc" : "Equivalence Model name.", \
    "uilabel" : "Name of the equivalence Model", }
  std::string eq_model;

  #pragma cyclus var {                    \
    "doc" : "equivalence Parameter line", \
    "uilabel" : "Parameter line for the equivalence model", }
  std::string eq_command;

///////// cycle params ///////////
  #pragma cyclus var {                                           \
    "doc" : "Discharge burnup.", "uilabel" : "Discharge burnup", \
    "units" : "GWd/t", }
  double burnup;

  #pragma cyclus var { \
    "doc": "The duration of a full operational cycle (excluding refueling " \
    "time) in time steps.", \
    "uilabel": "Cycle Length", \
    "units": "time steps", \
  }
  int cycle_time;

  #pragma cyclus var { \
    "doc": "The duration of a full refueling period - the minimum time between"\
    " the end of a cycle and the start of the next cycle.", \
    "uilabel": "Refueling Outage Duration", \
    "units": "time steps", \
  }
  int refuel_time;

  int cycle_step;
  int refueling_step;

//////////// power params ////////////
  #pragma cyclus var { \
    "doc": "Amount of thermal power the facility produces when operating " \
    "normally.", \
    "uilabel": "Thermal Reactor Power", \
    "units": "MWe", \
  }
  double power;

  #pragma cyclus var { \
    "default": 0, \
    "doc": "Amount of electrical power the facility produces when operating " \
    "normally.", \
    "uilabel": "Nominal Reactor Power", \
    "units": "MWe", \
  }
  double power_cap;

  #pragma cyclus var { \
    "default": "power", \
    "uilabel": "Power Commodity Name", \
    "doc": "The name of the 'power' commodity used in conjunction with a " \
    "deployment curve.", \
  }
  std::string power_name;

/////////// preference changes ///////////
  #pragma cyclus var { \
    "default": [], \
    "uilabel": "Time to Change Fresh Fuel Preference", \
    "doc": "A time step on which to change the request preference for a " \
    "particular fresh fuel type.", \
  }
  std::vector<int> pref_change_times;

  #pragma cyclus var { \
    "default": [], \
    "doc": "The input commodity for a particular fuel preference change.  " \
    "Same order as and direct correspondence to the specified " \
    "preference change times.", \
    "uilabel": "Commodity for Changed Fresh Fuel Preference", \
    "uitype": ["oneormore", "incommodity"], \
  }
  std::vector<std::string> pref_change_commods;

  #pragma cyclus var { \
    "default": [], \
    "uilabel": "Changed Fresh Fuel Preference",                        \
    "doc": "The new/changed request preference for a particular fresh fuel." \
    " Same order as and direct correspondence to the specified " \
    "preference change times.", \
  }
  std::vector<double> pref_change_values;

  // Resource inventories - these must be defined AFTER/BELOW the member vars
  // referenced (e.g. n_batch_fresh, assem_size, etc.).
  map<std::string, ResBuf<cyclus::Material>> fresh;
  map<std::string, ResBuf<cyclus::Material>> core;
  #pragma cyclus var {"capacity": "m_batch_spent"}
  ResBuf<cyclus::Material> spent;

// should be hidden in ui (internal only). True if fuel has already been
// discharged this cycle.
  vector<bool> discharged;

// This variable should be hidden/unavailable in ui.  Maps resource object
// id's to the index for the incommod through which they were received.
  #pragma cyclus var {                                           \
    "default" : {}, "doc" : "This should NEVER be set manually", \
               "internal" : True }
  std::map<int, int> res_indexes;

  // intra-time-step state - no need to be a state var
  // map<request, inventory name>
  std::map<cyclus::Request<cyclus::Material>*, std::string> req_inventories_;

  // populated lazily and no need to persist.
  std::set<std::string> uniq_outcommods_;
  };

  }  // namespace cyclass

#endif  // CYCLASS_SRC_REACTOR_H_
