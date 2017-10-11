#include "reactor.h"
#include <ostream>

using cyclus::Material;
using cyclus::Composition;
using cyclus::toolkit::ResBuf;
using cyclus::toolkit::MatVec;
using cyclus::KeyError;
using cyclus::ValueError;
using cyclus::Request;
using cyclus::Nuc;

typedef std::map<Nuc, double> CompMap;

namespace cyclass {

Reactor::Reactor(cyclus::Context* ctx)
    : cyclus::Facility(ctx),
      batch_size(0),
      n_batch_core(0),
      burnup(0),
      m_batch_spent(0),
      n_batch_fresh(0),
      cycle_time(0),
      refuel_time(0),
      cycle_step(0),
      refueling_step(-1),
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

cyclus::Inventories Reactor::SnapshotInv() {
  cyclus::Inventories invs;

  // these inventory names are intentionally convoluted so as to not clash
  // with the user-specified stream commods that are used as the separations
  // streams inventory names.
  invs["spent-inv-name"] = spent.PopNRes(spent.count());
  spent.Push(invs["spent-inv-name"]);
  std::map<std::string, ResBuf<Material> >::iterator it;
  for (it = fresh.begin(); it != fresh.end(); ++it) {
    invs[it->first] = it->second.PopNRes(it->second.count());
    it->second.Push(invs[it->first]);
  }
  for (it = core.begin(); it != core.end(); ++it) {
    invs[it->first] = it->second.PopNRes(it->second.count());
    it->second.Push(invs[it->first]);
  }
  return invs;
}
void Reactor::InitInv(cyclus::Inventories& inv) {
  spent.Push(inv["spent-inv-name"]);

  cyclus::Inventories::iterator it;
  for (it = inv.begin(); it != inv.end(); ++it) {
    std::map<std::string, ResBuf<Material> >::iterator it2;
    it2 = fresh.find(it->first);
    if (it2 != fresh.end()){
      fresh[it->first].Push(it->second);
    }
    it2 = core.find(it->first);
    if (it2 != core.end()){
      core[it->first].Push(it->second);
    }
  }
}

//________________________________________________________________________
void Reactor::EnterNotify() {
  cyclus::Facility::EnterNotify();

  MyCLASSAdaptator = new CLASSAdaptator(eq_model, eq_command, xs_model,
                                        xs_command, ir_model, ir_command);
  // initialisation internal variable
  for (int i = 0; i < n_batch_core; i++) {
    std::string batch_name = "batch_" + std::to_string(i);
    std::string fresh_name = "fresh_" + std::to_string(i);
    fresh[fresh_name].capacity(n_batch_fresh*batch_size);
    core[batch_name].capacity(batch_size);
    discharged.push_back(true);
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
bool Reactor::FullCore(){
  bool full_core = true;
  for (int i = 0; i < n_batch_core; i++) {
    std::string batch_name = "batch_" + std::to_string(i);
    if (core[batch_name].quantity() < batch_size){
      full_core = false;
    }
  }

  return full_core;
}

//________________________________________________________________________
bool Reactor::Discharged(){
  for (int i = 0; i < n_batch_core; i++){
    if (!discharged[i]){
      return false;
    }
  }
  return true;
}

//________________________________________________________________________
bool Reactor::Refueling(){
  if (refueling_step == -1){
    return false;
  } else {
    return true;
  }
}
//________________________________________________________________________
bool Reactor::InCycle(){
  if (Refueling() || !Discharged()) {
    return false;
  }
  return cycle_step >= 0 && FullCore();
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
  using cyclus::toolkit::RecordTimeSeries;
  using cyclus::toolkit::POWER;

  if (retired()) {
    Record("RETIRED", "");

    // record the last time series enter if the reactor was operating at the
    // time of retirement.
    if (exit_time() == context()->time()) {
      if (FullCore() && cycle_step > 0 && cycle_step <= cycle_time * n_batch_core ) {
        RecordTimeSeries<POWER>(this, power_cap);
      } else {
        RecordTimeSeries<POWER>(this, 0);
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
        if (!Discharge(i)) {
          break;
        }
      }
    }
    // in case a cycle lands exactly on our last time step, we will need to
    // burn a batch from fresh inventory on this time step.  When retired,
    // this batch also needs to be discharged to spent fuel inventory.
    for (int i = 0; i < n_batch_core; i++) {
    std::string fresh_name = "fresh_" + std::to_string(i);
      while (fresh[fresh_name].quantity() > 0 &&
             spent.space() >= m_batch_spent) {
        spent.Push(fresh[fresh_name].Pop());
      }
      return;
    }
  }

  if(FullCore()){
    if ( cycle_step !=0 && cycle_step % cycle_time == 0 && !Refueling()){
        int batch = cycle_step / cycle_time -1;
        Transmute(batch);
        discharged[batch] = false;
        std::stringstream ss;
        ss << "batch " << batch;
        Record("CYCLE_END", ss.str());
    }
  }


  if (cycle_step !=0 && cycle_step % cycle_time == 0 && !Discharged() ){
      int batch = cycle_step / cycle_time -1;
      discharged[batch] = Discharge(batch);
    }
  if (cycle_step % cycle_time == 0 && Discharged() && !Refueling()){
      int batch = cycle_step / cycle_time -1;
      while( batch < 0 ){
        batch += n_batch_core;
      }
      Load(batch);
      if (FullCore()){
        refueling_step = 0;
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
  if (retired()) {
    return ports;
  }

  Material::Ptr m;

  // second min expression reduces assembles to amount needed until
  // retirement if it is near.
  for (int u = 0; u < n_batch_core; u++) {
    // Deal first with the Core
    std::string batch_name = "batch_" + std::to_string(u);
    std::string fresh_name = "fresh_" + std::to_string(u);
    double mass_in_core_n = core[batch_name].quantity();
    double mass_to_order = batch_size - mass_in_core_n;

    if (mass_to_order != 0) {
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

        Composition::Ptr fissil_stream =
            Composition::CreateFromAtom(fissil_comp);
        Composition::Ptr fertil_stream =
            Composition::CreateFromAtom(fertil_comp);
        double required_burnup = burnup;
        if (context()->time() - enter_time() < cycle_time && !InCycle() &&
            u != 0) {
          required_burnup = burnup / n_batch_core * u;
        }
        double enrich = MyCLASSAdaptator->GetEnrichment(
            fissil_stream, fertil_stream, required_burnup);
        CompMap fuel_comp = fertil_comp * (1 - enrich) + fissil_comp * enrich;
        Composition::Ptr fuel = Composition::CreateFromAtom(fuel_comp);

        m = Material::CreateUntracked(mass_to_order, fuel);

        Request<Material>* r = port->AddRequest(m, this, commod, pref, true);
        req_inventories_[r] = batch_name;
        mreqs.push_back(r);
      }
      port->AddMutualReqs(mreqs);
      ports.insert(port);
    }

    double mass_fresh_n = fresh[fresh_name].quantity();

    mass_to_order = n_batch_fresh * batch_size - mass_fresh_n;

    // reduce the amount required if the reactor is about to be decommisionned
    if (exit_time() != -1) {
      // the +1 accounts for the fact that the reactor is alive and gets to
      // operate during its exit_time time step.
      int t_left = exit_time() - context()->time() + 1;
      int irradiation_time = cycle_step;
      if (context()->time() - enter_time() > n_batch_core * cycle_time) {
        irradiation_time += u * cycle_time;
        if (irradiation_time > n_batch_core * cycle_time) {
          irradiation_time -= n_batch_core * cycle_time;
        }
      }
      int t_left_cycle =
          cycle_time * n_batch_core + refuel_time - irradiation_time;

      double n_cycles_left =
          static_cast<double>(t_left - t_left_cycle) /
          static_cast<double>(cycle_time * n_batch_core + refuel_time);
      n_cycles_left = ceil(n_cycles_left);
      double mass_need =
          std::max(0.0, n_cycles_left * batch_size - mass_fresh_n);
      mass_to_order = std::min(mass_to_order, mass_need);
    }

    if (mass_to_order != 0) {
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

        Composition::Ptr fissil_stream =
            Composition::CreateFromAtom(fissil_comp);
        Composition::Ptr fertil_stream =
            Composition::CreateFromAtom(fertil_comp);

        double enrich = MyCLASSAdaptator->GetEnrichment(fissil_stream,
                                                        fertil_stream, burnup);
        CompMap fuel_comp = fertil_comp * (1 - enrich) + fissil_comp * enrich;
        Composition::Ptr fuel = Composition::CreateFromAtom(fuel_comp);

        m = Material::CreateUntracked(mass_to_order, fuel);

        Request<Material>* r = port->AddRequest(m, this, commod, pref, true);
        req_inventories_[r] = fresh_name;
        mreqs.push_back(r);
      }
      port->AddMutualReqs(mreqs);
      ports.insert(port);
    }
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
  double mass_in_core = 0;
  for (int i = 0; i < n_batch_core; i++) {
    std::string batch_name = "batch_" + std::to_string(i);
    mass_in_core += core[batch_name].quantity();
  }
  double m_responses = 0;
  for (int i = 0; i < (int)responses.size(); i++) {
    m_responses += responses[i].second->quantity();
  }
  double mload =
      std::min(m_responses, n_batch_core * batch_size - mass_in_core);
  if (mload > 0) {
    ss << mload << " kg ";
    Record("LOAD", ss.str());
  }
  for (trade = responses.begin(); trade != responses.end(); ++trade) {
    std::string batch_name = req_inventories_[trade->first.request];
    std::string commod = trade->first.request->commodity();

    bool fresh_r = false;
    if (batch_name.at(0) == 'f') {
      fresh_r = true;
    }

    Material::Ptr m = trade->second;
    index_res(m, commod);
    if (!fresh_r) {
      if (m->quantity() <= batch_size - core[batch_name].quantity()) {
        core[batch_name].Push(m);
        if (core[batch_name].quantity() == batch_size) {
          refueling_step = 0;
        }
      } else {
        fresh_r = true;
      }
    }
    if (fresh_r) {
      fresh[batch_name].Push(m);
    }
  }
  req_inventories_.clear();
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
  using cyclus::toolkit::RecordTimeSeries;
  using cyclus::toolkit::POWER;

  if (retired()) {
    return;
  }

  if (FullCore()) {
    if (refueling_step >= refuel_time) {
      refueling_step = -1;
    }

    if (!Refueling() && cycle_step % cycle_time == 0) {
      int batch = cycle_step / cycle_time - 1;
      std::stringstream ss;
      ss << "new cycle for batch " << batch;
      Record("CYCLE_START", ss.str());

      // if last batch reset cycle_step
      if ( batch == n_batch_core - 1) {
        cycle_step = 0;
      }
    }

    if (InCycle()) {
      RecordTimeSeries<POWER>(this, power_cap);
      cycle_step++;
    } else {
      RecordTimeSeries<POWER>(this, 0);
    }

    if (Refueling()) {
      refueling_step++;
    }
  }
}

//________________________________________________________________________
void Reactor::Transmute(int n_batch) {
  std::string batch_name = "batch_" + std::to_string(n_batch);
  MatVec old = core[batch_name].PopN(core[batch_name].count());
  core[batch_name].Push(old);

  if (!MyCLASSAdaptator) {
    MyCLASSAdaptator = new CLASSAdaptator(eq_model, eq_command, xs_model,
                                          xs_command, ir_model, ir_command);
  }

  double old_mass = 0;
  for (int i = 0; i < old.size(); i++) {
    old_mass += old[i]->quantity();
  }
  std::stringstream ss;
  ss << old_mass << " kg from batch " << n_batch;
  Record("TRANSMUTE", ss.str());

  // Get time inside the reactor
  // Get last batch number:
  int irradiation_time = cycle_step;
  if (context()->time() - enter_time() > n_batch_core * cycle_time) {
    irradiation_time += n_batch * cycle_time;
    if (irradiation_time > n_batch_core * cycle_time) {
      irradiation_time -= n_batch_core * cycle_time;
    }
  }

  for (int i = 0; i < old.size(); i++) {
    double mass = old[i]->quantity();
    double bu = power * irradiation_time / batch_size;
    cyclus::Composition::Ptr compo = old[i]->comp();
    old[i]->Transmute(
        MyCLASSAdaptator->GetCompAfterIrradiation(compo, power, mass, burnup));
  }
}

//________________________________________________________________________
std::map<std::string, MatVec> Reactor::PeekSpent() {
  std::map<std::string, MatVec> mapped;
    MatVec mats = spent.PopN(spent.count());
    spent.Push(mats);
    for (int i = 0; i < mats.size(); i++) {
      std::string commod = fuel_outcommod(mats[i]);
      mapped[commod].push_back(mats[i]);
    }
  return mapped;
}


//________________________________________________________________________
bool Reactor::Discharge(int i) {
  std::string batch_name = "batch_" + std::to_string(i);
  double npop = core[batch_name].quantity();
  if (m_batch_spent - spent.quantity() < npop) {
    Record("DISCHARGE", "failed");
    return false;  // not enough room in spent buffer
  }

  std::stringstream ss;
  ss << npop << " batches";
  Record("DISCHARGE", ss.str());

  spent.Push(core[batch_name].Pop(npop));
  return true;
}


//________________________________________________________________________
void Reactor::Load(int i) {
  std::string batch_name = "batch_" + std::to_string(i);
  std::string fresh_name = "fresh_" + std::to_string(i);
  double n = std::min(batch_size - core[batch_name].quantity(),
                      fresh[fresh_name].quantity());
  if (n != 0) {
    std::stringstream ss;
    ss << n << " batches";
    Record("LOAD", ss.str());
    core[batch_name].Push(fresh[fresh_name].Pop(n));
  }
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
std::map<std::string, MatVec > Reactor::PopSpent() {
  std::map<std::string, MatVec> mapped;

  MatVec mats = spent.PopN(spent.count());

  for (int i = 0; i < mats.size(); i++) {
    std::string commod = fuel_outcommod(mats[i]);
    mapped[commod].push_back(mats[i]);
  }

  // needed so we trade away oldest batches first
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
    // undo reverse in PopSpent to make sure oldest batches come out first
    std::reverse(it->second.begin(), it->second.end());
    spent.Push(it->second);
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
