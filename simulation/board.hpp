#ifndef BOARD_HPP
#define BOARD_HPP

#include <vector>

enum class Diag : char { s, i, r };
struct Diag_count {
  int sus;
  int inf;
  int rim;
  int sum() { return sus + rim + inf; }
};

struct Person {
  Diag diag;
  bool is_vaccinated;
  bool is_masked;
  int contact_time;
};
using line = std::vector<Person>;
class Board {
 private:
  const int dimension;
  const int cell_number;
  std::vector<Person> grid;
  const double beta_;
  const double gamma_;
  int count_time;
  Diag_count count;
  const int days;
  const double effectiveness;
  const int mask_percentage;
  const int vaccine_percentage;
  double t_max;
  double c;
  double a;

 public:
  Board(int s, int i, int r, double b, double y, int days, double e, int m,
        int v);
  void evolve();
  void mask_distribution();
  void mask_null();
  void vaccine_distribution();
  void vaccine_null();
  void counting(const Diag& diag);
  bool beta_random(int n, int t, bool b) const;
  bool gamma_random() const;
  void draw();
};

#endif
