#ifndef EPIDEMIC_HPP
#define EPIDEMIC_HPP

#include <vector>
template<typename T>
struct State {
  T sus_;
  T inf_;
  T rim_;
};


class Epidemic {
 private:
  std::vector<State<int>> round_memory;
  std::vector<State<double>> raw_memory;
  const double beta;
  const double gamma;
  const int tot;
  const int days;

 public:
  Epidemic(int s, int i, int r, double b, double y, int d);

  State<int> Evolve(const State<double>&  start);

  void Update();

  int S_get(int i) const;

  int I_get(int i) const;

  int R_get(int i) const;
};
#endif
