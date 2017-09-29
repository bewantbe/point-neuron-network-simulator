#include <iostream>
#include <random>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdio>
#include "../tictocs.cpp"
using std::cout;
using std::cerr;
using std::endl;

#define dbg_printf(...) ((void)0);
//#define dbg_printf printf

static const double Inf = std::numeric_limits<double>::infinity();

std::mt19937 rand_eng(1);

double g_rand()
{
  static std::uniform_real_distribution<double> udis(0, 1);
  return udis(rand_eng);
}

class TyPoissonSource {
  double rate, strength, t;
  std::exponential_distribution<> exp_dis;
public:
  double Init(double _rate, double _strength, double t0)
  {
    rate = _rate;
    strength = _strength;
    t = t0;
    exp_dis = *(new std::exponential_distribution<>(rate));
    return t = t0 + exp_dis(rand_eng);
  }

  TyPoissonSource()
  {
    Init(0.0, 0.0, 0.0);
  }

  TyPoissonSource(double _rate, double _strength, double t0)
    : rate(_rate), strength(_strength), exp_dis(_rate)
  {
    t = t0 + exp_dis(rand_eng);
  }

  double NextEventTime()
  {
    return t += exp_dis(rand_eng);
    //return t += rand_eng();
    //return t += g_rand();
    //return t += exp(g_rand());
  }

  inline double CurrentEventTime() const
  { return t; }

  inline double Strength() const
  { return strength; }
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
  void AddEventsUntilTime(TyTimeStrengthVec &event_vec, double rate, double strength, double t_until)
  {
    assert(rate >= 0);
    std::exponential_distribution<> exp_dis(rate);
    while (event_vec.back().time < t_until) {
      //event_vec.emplace_back(event_vec.back().time + exp_dis(rand_eng), strength);
      event_vec.emplace_back(poisson_src1.NextEventTime(), strength);
    }
  }

  TyPoissonSource poisson_src1, poisson_src2;

  // Fill poisson events upto t_until and no less than n_least_fill new events when possible
  void FillEvents(TyTimeStrengthVec &event_vec, double t_until, int n_least_fill)
  {
    // It should gurantee that event_vec was filled upto some t_until_old.
    // And that last status is still stored in corresponding PoissonSource.
    double tp[2] = {poisson_src1.CurrentEventTime(),
                    poisson_src2.CurrentEventTime()};
    int id_smaller = 1*(tp[1] <= tp[0]);
    dbg_printf("tp[%d] = %f\n", id_smaller, tp[id_smaller]);
    if (!std::isfinite(tp[id_smaller])) {
      // Both source give Inf
      if (event_vec.size() == 0 || std::isfinite(event_vec.back().time)) {
        // Add tailing Inf in these cases
        event_vec.emplace_back(Inf, 0);
      }
      return;
    }
    TyPoissonSource * const psrc[2] = {&poisson_src1, &poisson_src2};
    double sp[2] = {poisson_src1.Strength(), poisson_src2.Strength()};
    //if (!std::isfinite(tp[1-id_smaller])) {
      //// One source give Inf, the other give finite value
      //// Possible simpler algorithm
      //TyPoissonSource &poi_src = *psrc[id_smaller];
      //double strength = sp[id_smaller];
      //double t = poi_src.CurrentEventTime();
      //while (t < t_until || --n_least_fill >= 0) {
        //event_vec.emplace_back(t, strength);
        //t = poi_src.NextEventTime();
      //}
      //return;
    //}
    while (tp[id_smaller] < t_until || --n_least_fill >= 0) {
      dbg_printf("tp[%d] = %f\n", id_smaller, tp[id_smaller]);
      event_vec.emplace_back(tp[id_smaller], sp[id_smaller]);
      tp[id_smaller] = (psrc[id_smaller])->NextEventTime();
      id_smaller = 1*(tp[1] <= tp[0]);
    }
  }

  void FillEvents2(TyTimeStrengthVec &event_vec, double rate1, double strength1, double rate2, double strength2, double t_until, int n_least_fill)
  {
    // It should gurantee that event_vec was filled upto some t_until_old.
    // And that last status is still stored in corresponding PoissonSource.
    double t1 = poisson_src1.CurrentEventTime();
    double t2 = poisson_src2.CurrentEventTime();
    double t_smaller;
    int id_smaller;
    if (t1 < t2) {
      t_smaller = t1;
      id_smaller = 1;
    } else {
      t_smaller = t2;
      id_smaller = 2;
    }
    if (!std::isfinite(t_smaller)) {
      if (event_vec.size() == 0 || std::isfinite(event_vec.back().time)) {
        // Add tailing Inf in these cases
        event_vec.emplace_back(Inf, 0);
      }
      return;
    }
    while (t_smaller < t_until || --n_least_fill >= 0) {
      if (id_smaller == 1) {
        event_vec.emplace_back(t_smaller, strength1);
        t1 = poisson_src1.NextEventTime();
      } else {
        event_vec.emplace_back(t_smaller, strength2);
        t2 = poisson_src2.NextEventTime();
      }
      if (t1 < t2) {
        t_smaller = t1;
        id_smaller = 1;
      } else {
        t_smaller = t2;
        id_smaller = 2;
      }
    }
  }

  void FillEvents3(TyTimeStrengthVec &event_vec, double rate1, double strength1, double rate2, double strength2, double t_until, int n_least_fill)
  {
    // It should gurantee that event_vec was filled upto some t_until_old.
    // And that last status is still stored in corresponding PoissonSource.
    double t1 = poisson_src1.CurrentEventTime();
    double t2 = poisson_src2.CurrentEventTime();
    double t_smaller;
    if (t1 < t2) {
      t_smaller = t1;
    } else {
      t_smaller = t2;
    }
    if (!std::isfinite(t_smaller)) {
      if (event_vec.size() == 0 || std::isfinite(event_vec.back().time)) {
        // Add tailing Inf in these cases
        event_vec.emplace_back(Inf, 0);
      }
      return;
    }
    while (true) {
      if (t1 < t2) {
        if (t1 >= t_until && --n_least_fill < 0) break;
        event_vec.emplace_back(t1, strength1);
        t1 = poisson_src1.NextEventTime();
      } else {
        if (t2 >= t_until && --n_least_fill < 0) break;
        event_vec.emplace_back(t2, strength2);
        t2 = poisson_src2.NextEventTime();
      }
    }
  }

  void Init(double rate, double strength, double t0)
  {
    InitEventVec(*this, rate, strength, t0);
    id_seq = 0;
    poisson_src1.Init(rate, strength, t0);
  }

  void Init(double rate1, double strength1, double rate2, double strength2, double t0)
  {
    double t1 = poisson_src1.Init(rate1, strength1, t0);
    double t2 = poisson_src2.Init(rate2, strength2, t0);
    clear();
    if (t1 < t2) {
      emplace_back(t1, strength1);
      if (std::isfinite(t2)) {
        emplace_back(t2, strength2);
      }
    } else {
      emplace_back(t2, strength2);
      if (std::isfinite(t1)) {
        emplace_back(t1, strength1);
      }
    }
    poisson_src1.NextEventTime();
    poisson_src2.NextEventTime();
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
    dbg_printf("PopAndFill_1: id_seq=%lu, t=%f\n", id_seq, Front().time);
    id_seq++;
    if (id_seq == size()) {
      assert(size() > 0);
      if (std::isfinite(back().time)) {
        if (auto_shrink) {
          //double t0 = back().time;
          resize(0);
          //std::exponential_distribution<> exp_dis(rate);
          //emplace_back(t0 + exp_dis(rand_eng), strength);
          id_seq = 0;
          //Init(rate, strength, back().time);
        }
        AddEventsUntilTime(*this, rate, strength, back().time + 12.0 / rate);
      } else {
        id_seq--;
      }
    }
  }

  // Assume buf_event_vec1 and buf_event_vec2 each contains exactly one element.
  // Assume this->size()>1, the two tailing events should be belong to buf_event_vec1 and buf_event_vec2 each.
  void PopAndFill(double t_until, bool auto_shrink = false)
  {
    dbg_printf("PopAndFill_2: id_seq=%lu, t=%f\n", id_seq, Front().time);
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

  void Shrink()
  {
    assert(id_seq <= size());
    erase(begin(), begin()+id_seq);
    id_seq = 0;  // after Shrink(), one should check size()>0
  }
};

int main()
{
  size_t nit = (size_t)4e6;
  double s = 0;

  double r1 = 1.5;
  double s1 = 1.0;
  tic();
  s = 0;
  TyPoissonSource ps(r1, s1, 0.0);
  for (size_t i = 0; i < nit; i++) {
    s += ps.NextEventTime();
  }
  toc("poisson_source");
  cout << "________s = " << s << "\n";

  tic();
  s = 0;
  TyPoissonTimeSeq poisson_seq;
  poisson_seq.Init(r1, s1, 0.0);
  for (size_t i = 0; i < nit; i++) {
    s += poisson_seq.Front().time;
    poisson_seq.PopAndFill(r1, s1, true);
  }
  toc("E");
  cout << "________s = " << s << "\n";

  double r2 = 0.0;
  double s2 = 0.0;
  tic();
  s = 0;
  poisson_seq.Init(r1, s1, r2, s2, 0.0);
  for (size_t i = 0; i < nit; i++) {
    s += poisson_seq.Front().time;
    poisson_seq.PopAndFill(0.0, true);
  }
  toc("EI");
  cout << "________s = " << s << "\n";

  return 0;
}
