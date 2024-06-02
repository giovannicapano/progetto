#include "epidemic.hpp"

#include <cassert>
#include <cmath>
#include <random>
inline std::default_random_engine generator;
inline std::uniform_real_distribution<double> distr(0.0, 1.0);
Epidemic::Epidemic(int s, int i, int r, double b, double y, int d)
    : memory{State{s, i, r}}, beta{b}, gamma{y}, tot{s + i + r}, days{d} {}

State Epidemic::Evolve(State const &start) {
  if (start.sus_ < 0 || start.inf_ < 0 || start.rim_ < 0) {
    throw std::runtime_error{
        "il numero delle varie popolazioni non può essere negativo"};
  }
  if (days < 0) {
    throw std::runtime_error{"Il numero di infetti deve essere non negativo"};
  }
  if (beta < 0 || beta > 1) {
    throw std::runtime_error{
        "Beta deve essere non negativo e minore o uguale di 1"};
  }
  if (gamma < 0 || gamma > 1) {
    throw std::runtime_error{
        "Gamma deve essere non negativo e minore o uguale a 1"};
  }
  State end;
  double susd = start.sus_ - beta * start.sus_ * start.inf_ / tot;
  double infd =
      start.inf_ + beta * start.sus_ * start.inf_ / tot - gamma * start.inf_;
  double rimd = tot - start.sus_ - start.inf_ + gamma * start.inf_;
  end.sus_ = std::round(susd);
  end.inf_ = std::round(infd);
  end.rim_ = std::round(rimd);
  assert((end.sus_ >= 0));
  assert((end.inf_ >= 0));
  assert((end.rim_ >= 0));
  if (end.inf_ == start.inf_ && end.inf_ != 0) {
    double rand = distr(generator);
    end.inf_ = (rand < gamma) ? end.inf_ - 1 : end.inf_;
  }
  return end;
}
void Epidemic::Update() {
  for (int i = 0; i < days; ++i) {
    memory.push_back(Evolve(memory[i]));
  }
}

int Epidemic::S_get(int i) const { return memory[i].sus_; }
int Epidemic::I_get(int i) const { return memory[i].inf_; }
int Epidemic::R_get(int i) const { return memory[i].rim_; }
