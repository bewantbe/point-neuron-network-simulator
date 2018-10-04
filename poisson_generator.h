#ifndef HEADER_POISSON_GENERATOR
#define HEADER_POISSON_GENERATOR

#include "common_header.h"

/* For generating and store poisson events.
 * No idea why this simple task got so complex.
 */

// Used to generate poisson events
class TyPoissonSource
{
  double rate, strength, t;
  std::exponential_distribution<> exp_dis;
public:
  double CurrentEventTime() const
  { return t; }

  double Rate() const
  { return rate; }

  double Strength() const
  { return strength; }

  double NextEventTime()
  { return t += exp_dis(rand_eng); }

  void Set(double _rate, double _strength, double t0)
  {
    rate = _rate;
    strength = _strength;
    t = t0;

	if (rate == 0) rate = std::numeric_limits<double>::min();
    exp_dis = std::exponential_distribution<>(rate);
  }

  double Init(double _rate, double _strength, double t0)
  {
    assert(rate >= 0);
    Set(_rate, _strength, t0);
    if (_rate == 0) {
      t = Inf;
      return Inf;  // avoid a call to rand_eng
    } else {
      return NextEventTime();
    }
  }

  TyPoissonSource()
  { Set(0.0, 0.0, Inf); }

  TyPoissonSource(double _rate, double _strength, double t0)
    : rate(_rate), strength(_strength), t(t0), exp_dis(_rate)
  {
    if (rate == 0) {
      t = Inf;
    } else {
      NextEventTime();
    }
  }
};

struct EventTimeStrength
{
  double time, strength;
  EventTimeStrength() = default;
  EventTimeStrength(double _t, double _s)
    :time(_t), strength(_s)
  {}

  bool operator < (const EventTimeStrength &b) const
  { return time < b.time; }
  bool operator > (const EventTimeStrength &b) const
  { return time > b.time; }
  bool operator == (const EventTimeStrength &b) const
  { return time == b.time && strength == b.strength; }
};

typedef std::vector<EventTimeStrength> TyTimeStrengthVec;

// Used to generate Poisson events for one neuron.
class TyPoissonTimeSeq: public TyTimeStrengthVec
{
public:
  typedef size_t TyIdx;
  TyIdx id_seq = 0;      // point to current event
  const size_t EVENT_VEC_SIZE_LIMIT = 256;
  const size_t EVENT_GEN_CHUNK = 12;

  TyPoissonTimeSeq() = default;
  TyPoissonTimeSeq(const TyPoissonTimeSeq &pts)
  : TyTimeStrengthVec(pts),
    id_seq(pts.id_seq)
  {}

  void Shrink()
  {
    assert(id_seq <= size());
    erase(begin(), begin()+id_seq);
    id_seq = 0;
  }

  const EventTimeStrength & Front() const
  {
    return operator[](id_seq);
  }

  void InitEventVec(TyTimeStrengthVec &event_vec, double rate, double strength, double t0) const
  {
    assert(rate >= 0);
    event_vec.clear();
    std::exponential_distribution<> exp_dis(rate);
    event_vec.emplace_back(t0 + exp_dis(rand_eng), strength);
  }

  // Requirements:
  //   rate >= 0, rate < Inf
  //   t_until <= Inf
  //   event_vec.size() > 0
  void AddEventsUntilTime(TyTimeStrengthVec &event_vec, double rate, double strength, double t_until) const
  {
    assert(rate >= 0);
    std::exponential_distribution<> exp_dis(rate);
    while (event_vec.back().time < t_until) {
      event_vec.emplace_back(event_vec.back().time + exp_dis(rand_eng), strength);
    }
  }

  // Used for E type Poisson input only.
  void Init(double rate, double strength, double t0)
  {
    InitEventVec(*this, rate, strength, t0);
    id_seq = 0;
  }

  void PrintEvents()
  {
    for (size_t i = 0; i < size(); i++) {
      printf("event[%lu] = %.17g\n", i, (*this)[i].time);
    }
  }

  // For E type poisson input.
  // Use of the two versions of PopAndFill must NOT mixed.
  // If you do want to mix, call Init when change version.
  void PopAndFillOld(double rate, double strength, bool auto_shrink = false)
  {
    assert(id_seq < size());
    id_seq++;
    if (id_seq == size()) {
      assert(size() > 0);
      if (std::isfinite(back().time)) {
        if (auto_shrink) {
          Init(rate, strength, back().time);
        }
        AddEventsUntilTime(*this, rate, strength, back().time + EVENT_GEN_CHUNK / rate);
      } else {
        id_seq--;
      }
    }
  }

  // Poisson event generators for E and I type.
  TyPoissonSource poisson_src1, poisson_src2;

  // Fill poisson events upto t_until and no less than n_least_fill new
  // events when possible.
  void FillEvents(TyTimeStrengthVec &event_vec, double t_until, int n_least_fill)
  {
    // Assume the last status are stored in corresponding PoissonSources.
    double tp[2] = {poisson_src1.CurrentEventTime(),
                    poisson_src2.CurrentEventTime()};
    int id_smaller = 1*(tp[1] <= tp[0]);
    if (!std::isfinite(tp[id_smaller])) { // Both sources give Inf
      if (event_vec.size() == 0 || std::isfinite(event_vec.back().time)) {
        // Add tailing Inf when the queue do not have one.
        event_vec.emplace_back(Inf, 0);
      }
      return;
    }
    TyPoissonSource * const psrc[2] = {&poisson_src1, &poisson_src2};
    double sp[2] = {poisson_src1.Strength(), poisson_src2.Strength()};
    while (tp[id_smaller] < t_until || --n_least_fill >= 0) {
      event_vec.emplace_back(tp[id_smaller], sp[id_smaller]);
      tp[id_smaller] = psrc[id_smaller]->NextEventTime();
      id_smaller = 1*(tp[1] <= tp[0]);
    }
  }

  void FillEvents(double t_until, int n_least_fill)
  {
    FillEvents(*this, t_until, n_least_fill);
  }

  // Use for E and I type Poisson input.
  void Init(double rate1, double strength1, double rate2, double strength2, double t1, double t2)
  {
    /* // Old behaviour
    Init(rate1, strength1, t0);
    poisson_src1.Set(rate1, strength1, 0);
    */

    poisson_src1.Init(rate1,  strength1, t1);
    poisson_src2.Init(rate2, -strength2, t2);  // Negative for inhibitory.
    clear();
    id_seq = 0;
    FillEvents(-Inf, 1);  // Only one event is filled.
  }

  // Use for E and I type Poisson input.
  void PopAndFill(double t_until, bool auto_shrink = false)
  {
    id_seq++;
    if (id_seq == size()) {
      if (std::isfinite(back().time)) {
        //double t_end = back().time + EVENT_GEN_CHUNK / poisson_src1.Rate();
        //FillEvents(t_end, 0);
        if (auto_shrink && size() > EVENT_VEC_SIZE_LIMIT) {
          Shrink();
        }
        FillEvents(t_until, EVENT_GEN_CHUNK);
      } else {
        id_seq--;
      }
    }
  }

  void PopAndFill()
  {
    PopAndFill(-Inf);
    //PopAndFillOld(poisson_src1.Rate(), poisson_src1.Strength());  // Old behaviour
  }
};

// Used to generate Poisson events for each neurons
class TyPoissonTimeVec: public std::vector<TyPoissonTimeSeq>
{
  std::vector<TyPoissonTimeSeq::TyIdx> id_seq_vec;  // point to current event
public:
  void Init(const TyArrVals &arr_pr,  const TyArrVals &arr_ps,
            const TyArrVals &arr_pri, const TyArrVals &arr_psi, double t0)
  {
    // each TyPoissonTimeSeq should have at least one event
    resize(arr_pr.size());
    for (size_t j = 0; j < arr_pr.size(); j++) {
      operator[](j).Init(arr_pr[j], arr_ps[j], arr_pri[j], arr_psi[j], t0, t0);
    }
    id_seq_vec.resize(arr_pr.size());
    std::fill(id_seq_vec.begin(), id_seq_vec.end(), 0);
  }

  TyPoissonTimeVec() = default;  // enable all other default constructors

  void RemoveEvents()
  {
    assert(size() == id_seq_vec.size());
    for (auto &i : *this) {
      i.clear();
      i.id_seq = 0;
    }
    std::fill(id_seq_vec.begin(), id_seq_vec.end(), 0);
  }

  // require: call Init first.
  void FillEvents(double t_until, int n_least_fill)
  {
    for (size_t j = 0; j < size(); j++) {
      operator[](j).FillEvents(t_until, n_least_fill);
    }
  }
  
  void ClearAndContinueFill(double t_until, int n_least_fill)
  {
//    std::vector<TyPoissonSource> poisson_src1_vec(size());
//    std::vector<TyPoissonSource> poisson_src2_vec(size());
//    for (size_t j = 0; j < size(); j++) {
//      poisson_src1_vec[j] = operator[](j).poisson_src1;
//      poisson_src2_vec[j] = operator[](j).poisson_src2;
//    }
    RemoveEvents();
    FillEvents(t_until, n_least_fill);
  }

  size_t EventCount() const
  {
    size_t cnt = 0;
    for (size_t j = 0; j < size(); j++) {
      cnt += operator[](j).size();
    }
    return cnt;
  }

  void RestoreIdx()
  {
    for (size_t j = 0; j < size(); j++) {
      operator[](j).id_seq = id_seq_vec[j];
    }
  }

  void RestoreIdx(const std::vector<int> &ids)
  {
    for (const int &id : ids) {
      operator[](id).id_seq = id_seq_vec[id];
    }
  }

  void SaveIdxAndClean()
  {
    int j = 0;
    for (iterator it = begin(); it != end(); ++it, ++j) {
      if (it->size() - it->id_seq < it->id_seq / 7) {  // The factor here is non-critical
        it->Shrink();
        dbg_printf("  ( %d-th TyPoissonTimeSeq size shrinked )\n", j);
      }
      id_seq_vec[j] = it->id_seq;
    }
  }

  void SaveIdxAndClean(const std::vector<int> &ids)
  {
    for (const int &id : ids) {
      auto it = begin() + id;
      if (it->size() - it->id_seq < it->id_seq / 7) {  // The factor here is non-critical
        it->Shrink();
        dbg_printf("  ( %d-th TyPoissonTimeSeq size shrinked )\n", id);
      }
      id_seq_vec[id] = it->id_seq;
    }
  }
};

size_t FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps);
void SavePoissonInput(std::ofstream &fout, TyPoissonTimeVec &poisson_time_vec, double t_step_end);

#endif
