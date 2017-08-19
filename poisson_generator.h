#ifndef HEADER_POISSON_GENERATOR
#define HEADER_POISSON_GENERATOR

#include "common_header.h"

// Used to generate poisson events
class TyPoissonSource {
  double rate, strength, t;
  std::exponential_distribution<> exp_dis;
public:
  double CurrentEventTime() const
  { return t; }

  double Strength() const
  { return strength; }

  double NextEventTime()
  { return t += exp_dis(rand_eng); }

  double Init(double _rate, double _strength, double t0)
  {
    assert(rate >= 0);
    rate = _rate;
    strength = _strength;
    t = t0;
    exp_dis = std::exponential_distribution<>(rate);
    return NextEventTime();
  }

  TyPoissonSource()
  { Init(0.0, 0.0, 0.0); }

  TyPoissonSource(double _rate, double _strength, double t0)
    : rate(_rate), strength(_strength), t(t0), exp_dis(_rate)
  { NextEventTime(); }
};

struct EventTimeStrength {
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

// Used to generate Poisson events for one neuron
class TyPoissonTimeSeq: public TyTimeStrengthVec
{
public:
  typedef size_t TyIdx;
  TyIdx id_seq = 0;      // point to current event
  size_t EVENT_VEC_SIZE_LIMIT = 256;

  TyPoissonTimeSeq() = default;
  TyPoissonTimeSeq(const TyPoissonTimeSeq &pts)
  : TyTimeStrengthVec(pts),
    id_seq(pts.id_seq)
  {}

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
      tp[id_smaller] = (psrc[id_smaller])->NextEventTime();
      id_smaller = 1*(tp[1] <= tp[0]);
    }
  }

  // Use for E and I type Poisson input.
  void Init(double rate1, double strength1, double rate2, double strength2, double t0)
  {
    poisson_src1.Init(rate1, strength1, t0);
    poisson_src2.Init(rate2, strength2, t0);
    clear();
    id_seq = 0;
    FillEvents(*this, t0, 1);  // t_until is not important here
  }

  // Use for E type Poisson input.
  void Init(double rate, double strength, double t0)
  {
    InitEventVec(*this, rate, strength, t0);
    id_seq = 0;
  }

  void PopAndFill(double t_until, bool auto_shrink = false)
  {
    id_seq++;
    if (id_seq == size()) {
      if (std::isfinite(back().time)) {
        if (auto_shrink && size() > EVENT_VEC_SIZE_LIMIT) {
          Shrink();
        }
        FillEvents(*this, t_until, 12);
      } else {
        id_seq--;
      }
    }
  }

  const EventTimeStrength & Front() const
  {
    return operator[](id_seq);
  }

  // For E type poisson input.
  // Use of the two versions of PopAndFill must NOT mixed.
  // If you do want to mix, call Init when change version.
  void PopAndFill(double rate, double strength, bool auto_shrink = false)
  {
    assert(id_seq < size());
    id_seq++;
    if (id_seq == size()) {
      assert(size() > 0);
      if (std::isfinite(back().time)) {
        if (auto_shrink) {
          Init(rate, strength, back().time);
        }
        AddEventsUntilTime(*this, rate, strength, back().time + 12.0 / rate);
      } else {
        id_seq--;
      }
    }
  }

  void Shrink()
  {
    assert(id_seq < size());
    erase(begin(), begin()+id_seq);
    id_seq = 0;
  }
};

// Used to generate Poisson events for each neurons
class TyPoissonTimeVec: public std::vector<TyPoissonTimeSeq>
{
  std::vector<TyPoissonTimeSeq::TyIdx> id_seq_vec;  // point to current event
public:
  void Init(const TyArrVals &rate_vec, const TyArrVals &arr_ps, double t0)
  {
    // each TyPoissonTimeSeq should have at least one event
    resize(rate_vec.size());
    for (size_t j = 0; j < rate_vec.size(); j++) {
      operator[](j).Init(rate_vec[j], arr_ps[j], t0);
    }
    id_seq_vec.resize(rate_vec.size());
  }

  TyPoissonTimeVec(const TyArrVals &rate_vec, const TyArrVals &arr_ps, double t0 = 0)
  {
    Init(rate_vec, arr_ps, t0);
  }

  TyPoissonTimeVec() = default;  // enable all other default constructors

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

void FillPoissonEventsFromFile(TyPoissonTimeVec &poisson_time_vec, const char *path, const TyArrVals &arr_ps);
void SavePoissonInput(std::ofstream &fout, TyPoissonTimeVec &poisson_time_vec, double t_step_end, const TyArrVals &arr_pr, const TyArrVals &arr_ps);

#endif
