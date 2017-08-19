#ifndef HEADER_POISSON_GENERATOR
#define HEADER_POISSON_GENERATOR

#include "common_header.h"

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

  // Used for generating excititory-inhibitory mixed events.
  TyTimeStrengthVec buf_event_vec1, buf_event_vec2;

  void InitEventVec(TyTimeStrengthVec &event_vec, double rate, double strength, double t0) const
  {
    assert(rate >= 0);
    event_vec.clear();
    std::exponential_distribution<> exp_dis(rate);
    event_vec.emplace_back(t0 + exp_dis(rand_eng), strength);
  }

  // Used before refill an event_vec.
  void KeepOnlyLatestEvent(TyTimeStrengthVec &event_vec) const
  {
    event_vec[0] = event_vec.back();
    event_vec.resize(1);
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

  void RefillEvents(TyTimeStrengthVec &event_vec, double rate, double strength, double t_until) const
  {
  }


  void Init(double rate, double strength, double t0)
  {
    InitEventVec(*this, rate, strength, t0);
    id_seq = 0;
  }

  void Init(double rate1, double strength1, double rate2, double strength2, double t0)
  {
    InitEventVec(buf_event_vec1, rate1, strength1, t0);
    InitEventVec(buf_event_vec2, rate2, strength2, t0);
    this->resize(2);
    std::merge(buf_event_vec1.begin(), buf_event_vec1.end(),
               buf_event_vec2.begin(), buf_event_vec2.end(),
               this->begin());
    id_seq = 0;
  }

  const EventTimeStrength & Front() const
  {
    return operator[](id_seq);
  }

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

  // Assume buf_event_vec1 and buf_event_vec2 each contains exactly one element.
  // Assume this->size()>1, the two tailing events should be belong to buf_event_vec1 and buf_event_vec2 each.
  void PopAndFill(double rate1, double strength1, double rate2, double strength2, double t_until, bool auto_shrink = false)
  {
    assert(size() >= 2);
    assert(id_seq + 2 <= size());
    if (id_seq + 2 == size()) {  // only two elements left
      assert(buf_event_vec1.size() == 1);
      assert(buf_event_vec2.size() == 1);
      if (auto_shrink && size() > EVENT_VEC_SIZE_LIMIT) {
        Shrink();
      }
      // Generate new poisson events until t_until.
      AddEventsUntilTime(buf_event_vec1, rate1, strength1, t_until);
      AddEventsUntilTime(buf_event_vec2, rate2, strength2, t_until);
      // merge to main event vec
      this->resize(id_seq + buf_event_vec1.size() + buf_event_vec2.size());
      std::merge(buf_event_vec1.begin(), buf_event_vec1.end(),
                 buf_event_vec2.begin(), buf_event_vec2.end(),
                 this->begin()+id_seq);
      /*if (buf_event_vec1.size() == 1 && buf_event_vec1*/
      // TODO: when no new event added.
      // clear buf
      KeepOnlyLatestEvent(buf_event_vec1);
      KeepOnlyLatestEvent(buf_event_vec2);
    }
    id_seq++;
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
