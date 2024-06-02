#ifndef EPIDEMIC_HPP
#define EPIDEMIC_HPP

#include <vector>

struct State {
  int sus_;
  int inf_;
  int rim_;
};

class Epidemic {
 private:
  std::vector<State> memory;
  const double beta;
  const double gamma;
  const int tot;
  const int days;

 public:
  Epidemic(int s, int i, int r, double b, double y, int d);

  State Evolve(State const &start);

  void Update();

  int S_get(int i) const;

  int I_get(int i) const;

  int R_get(int i) const;
};
#endif
