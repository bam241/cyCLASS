#include "reactor.h"

using cyclus::Material;
using cyclus::Composition;
using cyclus::toolkit::ResBuf;
using cyclus::toolkit::MatVec;
using cyclus::KeyError;
using cyclus::ValueError;
using cyclus::Request;
using cyclus::Nuc;

#include <ostream>
typedef std::map<Nuc, double> CompMap;

//#define cyDBGL		std::cout << __FILE__ << " : " << __LINE__ << "
//["
//<< __FUNCTION__ << "]" << std::endl;
#define;

namespace cyclass {

Reactor::Reactor(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      n_assem_batch(0),
      assem_size(0),
      n_assem_core(0),
      burnup(0),
      n_assem_spent(0),
      n_assem_fresh(0),
      cycle_time(0),
      refuel_time(0),
      cycle_step(0),
      power_cap(0),
      power_name("power"),
      xs_model("xs"),
      xs_command("xs"),
      ir_model("ir"),
      ir_command("ir"),
      eq_model("eq"),
      eq_command("eq"),
      discharged(false) {
  cyclus::Warn<cyclus::EXPERIMENTAL_WARNING>(
      "the Reactor archetype "
      "is experimental");

  MyCLASSAdaptator = 0;
}

#pragma cyclus def clone cyclass::Reactor

#pragma cyclus def schema cyclass::Reactor

#pragma cyclus def annotations cyclass::Reactor

#pragma cyclus def infiletodb cyclass::Reactor

#pragma cyclus def snapshot cyclass::Reactor

#pragma cyclus def snapshotinv cyclass::Reactor

#pragma cyclus def initinv cyclass::Reactor

//________________________________________________________________________
void Reactor::InitFrom(Reactor* m) {
#pragma cyclus impl initfromcopy cyclass::Reactor
  cyclus::toolkit::CommodityProducer::Copy(m);
}

//________________________________________________________________________
void Reactor::InitFrom(cyclus::QueryableBackend* b) {
#pragma cyclus impl initfromdb cyclass::Reactor

  namespace tk = cyclus::toolkit;
  tk::CommodityProducer::Add(tk::Commodity(power_name),
                             tk::CommodInfo(power_cap, power_cap));
}

//________________________________________________________________________
void Reactor::EnterNotify() {
  cyclus::Facility::EnterNotify();

  MyCLASSAdaptator = new CLASSAdaptator(eq_model, eq_command, xs_model,
                                        xs_command, ir_model, ir_command);
  // initialisation internal variable
  for (int i = 0; i < n_batches; i++) {
    std::string batch_name = "batch_" + std::to_string(i);
    fresh[batch_name].capacity(n_batch_fresh);
    core[batch_name].capacity(batch_size);
    spent[batch_name].capacity(n_batch_spent);
    cycle_step.push_back((n_batches - i) * cycle_time);
  }

  // If the user ommitted fuel_prefs, we set it to zeros for each fuel
  // type.  Without this segfaults could occur - yuck.
  if (fuel_prefs.size() == 0) {
    for (int i = 0; i < fuel_outcommods.size(); i++) {
      fuel_prefs.push_back(cyclus::kDefaultPref);
    }
  }
}
//________________________________________________________________________
bool Reactor::CheckDecommissionCondition() {
  return core[batch_name].quantity() == 0 && spent[batch_name].quantity() == 0;
}

bool Reactor::InCycle(i){
  std::string batch_name = "batch_" + std::to_string(i);

  return cycle_step[i] >= 0 && cycle_step[i] < cycle_time* n_batch_core &&
         core[batch_name].quantity() == batch_size;
}

//________________________________________________________________________
void Reactor::Tick() {
  // The following code must go in the Tick so they fire on the time step
  // following the cycle_step update - allowing for the all reactor events
  // to
  // occur and be recorded on the "beginning" of a time step.  Another
  // reason
  // they
  // can't go at the beginnin of the Tock is so that resource exchange has a
  // chance to occur after the discharge on this same time step.
  using cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>;

  if (retired()) {
    Record("RETIRED", "");

    // record the last time series entry if the reactor was operating at the
    // time of retirement.
    if (exit_time() == context()->time()) {
      if (cycle_step > 0 && cycle_step <= cycle_time &&
          core[batch_name].quantity() == n_assem_core) {
        RecordTimeSeries<cyclus::toolkit::POWER>(this, power_cap);
      } else {
        RecordTimeSeries<cyclus::toolkit::POWER>(this, 0);
      }
    }

    if (context()->time() == exit_time()) {  // only need to transmute once
      for (int i = 0; i < n_batch_core; i++) {
        Transmute(i);
      }
    }
    for (int i = 0; i < n_batch_core; i++) {
      std::string batch_name = "batch_" + std::to_string(i);
      while (core[batch_name].quantity() > 0) {
        if (!Discharge()) {
          break;
        }
      }
    }
    // in case a cycle lands exactly on our last time step, we will need to
    // burn a batch from fresh inventory on this time step.  When retired,
    // this batch also needs to be discharged to spent fuel inventory.
    for (int i = 0; i < n_batch_core; i++) {
      std::string batch_name = "batch_" + std::to_string(i);
      while (fresh[batch_name].quantity() > 0 &&
             spent[batch_name].space() >= batch_size) {
        spent[batch_name].Push(fresh[batch_name].Pop());
      }
      return;
    }
  }

  for (int i = 0; i < n_batch_core; i++) {
    if (cycle_step[i] == cycle_time * n_batch_core) {
      Transmute(i);
      Record("CYCLE_END", "");
    }

    if (cycle_step[i] >= cycle_time * n_batch_core && !discharged[i]) {
      discharged[i] = Discharge(i);
    }
    if (cycle_step[i] >= cycle_time * n_batch_core) {
      Load(i);
    }
  }

  int t = context()->time();

  // update preferences
  for (int i = 0; i < pref_change_times.size(); i++) {
    int change_t = pref_change_times[i];
    if (t != change_t) {
      continue;
    }

    std::string incommod = pref_change_commods[i];
    for (int j = 0; j < fuel_incommods.size(); j++) {
      if (fuel_incommods[j] == incommod) {
        fuel_prefs[j] = pref_change_values[i];
        break;
      }
    }
  }
}

//________________________________________________________________________
std::set<cyclus::RequestPortfolio<Material>::Ptr> Reactor::GetMatlRequests() {
  using cyclus::RequestPortfolio;

  std::set<RequestPortfolio<Material>::Ptr> ports;
  Material::Ptr m;

  // second min expression reduces assembles to amount needed until
  // retirement if it is near.
  int n_assem_order =
      n_assem_core - core[batch_name].quantity() + n_assem_fresh - fresh[batch_name].quantity();

  if (exit_time() != -1) {
    // the +1 accounts for the fact that the reactor is alive and gets to
    // operate during its exit_time time step.
    int t_left = exit_time() - context()->time() + 1;
    int t_left_cycle = cycle_time + refuel_time - cycle_step;
    double n_cycles_left = static_cast<double>(t_left - t_left_cycle) /
                           static_cast<double>(cycle_time + refuel_time);
    n_cycles_left = ceil(n_cycles_left);
    int n_need = std::max(0.0, n_cycles_left * n_assem_batch - n_assem_fresh +
                                   n_assem_core - core[batch_name].quantity());
    n_assem_order = std::min(n_assem_order, n_need);
  }

  if (n_assem_order == 0) {
    return ports;
  } else if (retired()) {
    return ports;
  }

  for (int i = 0; i < n_assem_order; i++) {
    RequestPortfolio<Material>::Ptr port(new RequestPortfolio<Material>());
    std::vector<Request<Material>*> mreqs;

    for (int j = 0; j < fuel_incommods.size(); j++) {
      std::string commod = fuel_incommods[j];
      double pref = fuel_prefs[j];

      // Defining the Stream composition

      CompMap fissil_comp;
      if (fuel_type == "mox") {
        fissil_comp.insert(std::pair<Nuc, double>(942390000, 100));
      } else if (fuel_type == "uox") {
        fissil_comp.insert(std::pair<Nuc, double>(922350000, 100));
      }
      fissil_comp = NormalizeComp(fissil_comp);
      CompMap fertil_comp;
      fertil_comp.insert(std::pair<Nuc, double>(922380000, 100));
      fertil_comp = cyclass::NormalizeComp(fertil_comp);

      Composition::Ptr fissil_stream = Composition::CreateFromAtom(fissil_comp);
      Composition::Ptr fertil_stream = Composition::CreateFromAtom(fertil_comp);

      double Enrch =
          MyCLASSAdaptator->GetEnrichment(fissil_stream, fertil_stream, burnup);
      // std::cout << "Enrich required: " << Enrch << std::endl;
      Composition::Ptr fuel = Composition::CreateFromAtom(
          fertil_comp * (1 - Enrch) + fissil_comp * Enrch);

      m = Material::CreateUntracked(assem_size, fuel);

      Request<Material>* r = port->AddRequest(m, this, commod, pref, true);
      mreqs.push_back(r);
    }
    port->AddMutualReqs(mreqs);
    ports.insert(port);
  }

  return ports;
}

//________________________________________________________________________
void Reactor::GetMatlTrades(
    const std::vector<cyclus::Trade<Material> >& trades,
    std::vector<std::pair<cyclus::Trade<Material>, Material::Ptr> >&
        responses) {
  using cyclus::Trade;

  std::map<std::string, MatVec> mats = PopSpent();
  for (int i = 0; i < trades.size(); i++) {
    std::string commod = trades[i].request->commodity();
    Material::Ptr m = mats[commod].back();
    mats[commod].pop_back();
    responses.push_back(std::make_pair(trades[i], m));
    res_indexes.erase(m->obj_id());
  }
  PushSpent(mats);  // return leftovers back to spent buffer
}

//________________________________________________________________________
void Reactor::AcceptMatlTrades(
    const std::vector<std::pair<cyclus::Trade<Material>, Material::Ptr> >&
        responses) {
  std::vector<std::pair<cyclus::Trade<cyclus::Material>,
                        cyclus::Material::Ptr> >::const_iterator trade;

  std::stringstream ss;
  int nload = std::min((int)responses.size(), n_assem_core - core[batch_name].quantity());
  if (nload > 0) {
    ss << nload << " assemblies";
    Record("LOAD", ss.str());
  }

  for (trade = responses.begin(); trade != responses.end(); ++trade) {
    std::string commod = trade->first.request->commodity();
    Material::Ptr m = trade->second;
    index_res(m, commod);

    if (core[batch_name].quantity() < n_assem_core) {
      core[batch_name].Push(m);
    } else {
      fresh[batch_name].Push(m);
    }
  }
}

//________________________________________________________________________
std::set<cyclus::BidPortfolio<Material>::Ptr> Reactor::GetMatlBids(
    cyclus::CommodMap<Material>::type& commod_requests) {
  using cyclus::BidPortfolio;

  std::set<BidPortfolio<Material>::Ptr> ports;

  std::map<std::string, MatVec> all_mats;

  if (uniq_outcommods_.empty()) {
    for (int i = 0; i < fuel_outcommods.size(); i++) {
      uniq_outcommods_.insert(fuel_outcommods[i]);
    }
  }

  std::set<std::string>::iterator it;
  for (it = uniq_outcommods_.begin(); it != uniq_outcommods_.end(); ++it) {
    std::string commod = *it;
    std::vector<Request<Material>*>& reqs = commod_requests[commod];
    if (reqs.size() == 0) {
      continue;
    } else {
      all_mats = PeekSpent();
    }

    MatVec mats = all_mats[commod];
    if (mats.size() == 0) {
      continue;
    }

    BidPortfolio<Material>::Ptr port(new BidPortfolio<Material>());

    for (int j = 0; j < reqs.size(); j++) {
      Request<Material>* req = reqs[j];
      double tot_bid = 0;
      for (int k = 0; k < mats.size(); k++) {
        Material::Ptr m = mats[k];
        tot_bid += m->quantity();
        port->AddBid(req, m, this, true);
        if (tot_bid >= req->target()->quantity()) {
          break;
        }
      }
    }

    double tot_qty = 0;
    for (int j = 0; j < mats.size(); j++) {
      tot_qty += mats[j]->quantity();
    }
    cyclus::CapacityConstraint<Material> cc(tot_qty);
    port->AddConstraint(cc);
    ports.insert(port);
  }

  return ports;
}

//________________________________________________________________________
void Reactor::Tock() {
  using cyclus::toolkit::RecordTimeSeries<cyclus::toolkit::POWER>;

  if (retired()) {
    return;
  }
  for (int i = 0; i < n_batches; i++) {
    std::string batch_name = "batch_" + std::to_string(i);

    if (cycle_step >= cycle_time *n_batch_core + refuel_time &&
        core[batch_name].quantity() == batch_size) {
      discharged = false;
      cycle_step = 0;
    }

    if (cycle_step == 0 && core[batch_name].quantity() == n_batch_core) {
      Record("CYCLE_START", "");
    }

    if (InCycle(i)) {
      RecordTimeSeries<cyclus::toolkit::POWER>(this, power_cap);
    } else {
      RecordTimeSeries<cyclus::toolkit::POWER>(this, 0);
    }

    // "if" prevents starting cycle after initial deployment until core is full
    // even though cycle_step is its initial zero.
    if (cycle_step[i] > 0 || core[batch_name].quantity() == batch_size) {
      cycle_step[i]++;
    }
  }
}


//________________________________________________________________________
void Reactor::Transmute(int n_batch) {
  std::string batch_name = "batch_" + std::to_string(i);
  MatVec old = core[batch_name].Pop(std::min(core[batch_name].quantity()));
  core[batch_name].Push(old);

  if (!MyCLASSAdaptator) {
    MyCLASSAdaptator = new CLASSAdaptator(eq_model, eq_command, xs_model,
                                          xs_command, ir_model, ir_command);
  }


  std::stringstream ss;
  ss << old.quantity() << " kg from batch " << n_batch;
  Record("TRANSMUTE", ss.str());

  for (int i = 0; i < old.size(); i++) {V
    double mass = old[i]->quantity();
    double bu = power*cycle_step[n_batch]/batch_size;
    cyclus::Composition::Ptr compo = old[i]->comp();
    old[i]->Transmute(
        MyCLASSAdaptator->GetCompAfterIrradiation(compo, power, mass, burnup));
  }
}

//________________________________________________________________________
std::map<std::string, MatVec> Reactor::PeekSpent() {
  std::map<std::string, MatVec> mapped;
  for (int j = 0; j < n_batch_core; j++){
    std::string batch_name = "batch_" + std::to_string(j);
    MatVec mats = spent[batch_name].Pop(spent[batch_name].quantity());
    spent[batch_name].Push(mats);
    for (int i = 0; i < mats.size(); i++) {
      std::string commod = fuel_outcommod(mats[i]);
      mapped[commod].push_back(mats[i]);
    }
  }
  return mapped;
}

//________________________________________________________________________
bool Reactor::Discharge() {
  int npop = std::min(n_assem_batch, core[batch_name].quantity());
  if (n_assem_spent - spent[batch_name].quantity() < npop) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }

  std::stringstream ss;
  ss << npop << " assemblies";
  Record("DISCHARGE", ss.str());

  spent[batch_name].Push(core[batch_name].Pop(npop));
  return true;
}

//________________________________________________________________________
void Reactor::Load() {
  int n = std::min(n_assem_core - core[batch_name].quantity(), fresh[batch_name].quantity());
  if (n == 0) {
    return;
  }

  std::stringstream ss;
  ss << n << " assemblies";
  Record("LOAD", ss.str());
  core[batch_name].Push(fresh[batch_name].Pop(n));
}

//________________________________________________________________________
std::string Reactor::fuel_incommod(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_incommods.size()) {
    throw KeyError("cyclass::Reactor - no incommod for material object");
  }
  return fuel_incommods[i];
}

//________________________________________________________________________
std::string Reactor::fuel_outcommod(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_outcommods.size()) {
    throw KeyError("cyclass::Reactor - no outcommod for material object");
  }
  return fuel_outcommods[i];
}

//________________________________________________________________________
double Reactor::fuel_pref(Material::Ptr m) {
  int i = res_indexes[m->obj_id()];
  if (i >= fuel_prefs.size()) {
    return 0;
  }
  return fuel_prefs[i];
}

//________________________________________________________________________
void Reactor::index_res(cyclus::Resource::Ptr m, std::string incommod) {
  for (int i = 0; i < fuel_incommods.size(); i++) {
    if (fuel_incommods[i] == incommod) {
      res_indexes[m->obj_id()] = i;
      return;
    }
  }
  throw ValueError("cyclass::Reactor - received unsupported incommod material");
}

//________________________________________________________________________
std::map<std::string, MatVec> Reactor::PopSpent() {
    MatVec mats = spent[batch_name].Pop(spent[batch_name].quantity());
    std::map<std::string, MatVec> mapped;
    for (int i = 0; i < mats.size(); i++) {
      std::string commod = fuel_outcommod(mats[i]);
      mapped[commod].push_back(mats[i]);
    }

  // needed so we trade away oldest assemblies first
  std::map<std::string, MatVec>::iterator it;
  for (it = mapped.begin(); it != mapped.end(); ++it) {
    std::reverse(it->second.begin(), it->second.end());
  }
  return mapped;
}

//________________________________________________________________________
void Reactor::PushSpent(std::map<std::string, MatVec> leftover) {
  std::map<std::string, MatVec>::iterator it;
  for (it = leftover.begin(); it != leftover.end(); ++it) {
    // undo reverse in PopSpent to make sure oldest assemblies come out first
    std::reverse(it->second.begin(), it->second.end());
    spent[batch_name].Push(it->second);
  }
}

//________________________________________________________________________
void Reactor::Record(std::string name, std::string val) {
  context()
      ->NewDatum("ReactorEvents")
      ->AddVal("AgentId", id())
      ->AddVal("Time", context()->time())
      ->AddVal("Event", name)
      ->AddVal("Value", val)
      ->Record();
}

extern "C" cyclus::Agent* ConstructReactor(cyclus::Context* ctx) {
  return new Reactor(ctx);
}

}  // namespace cyclass
