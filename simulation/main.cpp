#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

#include "board.hpp"
#include "epidemic.hpp"
bool is_number(std::string& str) {
  for (char ch : str) {
    if (!(ch >= 48 && ch <= 57)) {
      if (ch != 46 && ch != 44 && ch != 43 && ch != 45) {
        return false;
      }
    }
  }
  return true;
};
bool is_integer(std::string& str) {
  for (char ch : str) {
    if (!(ch >= 48 && ch <= 57)) {
      if (ch != 43 && ch != 45) {
        return false;
      }
    }
  }
  return true;
};

int main() {
  std::cout << "Per popolazioni maggiori di 100000, la simulazione grafica "
               "potrebbe risultare lenta."
            << '\n';

  std::string n_ss;
  std::cout << "Inserisci il numero di suscettibili:";
  std::cin >> n_ss;
  while (is_integer(n_ss) == false) {
    std::cout << "L'input deve essere un numero intero: ";
    std::cin >> n_ss;
  }
  int n_s = std::stoi(n_ss);
  std::string n_ii;
  std::cout << "Inserisci il numero di infetti:";
  std::cin >> n_ii;
  while (is_integer(n_ii) == false) {
    std::cout << "L'input deve essere un numero intero: ";
    std::cin >> n_ii;
  }
  int n_i = std::stoi(n_ii);
  std::string n_rr;
  std::cout << "Inserisci il numero di rimossi:";
  std::cin >> n_rr;
  while (is_integer(n_rr) == false) {
    std::cout << "L'input deve essere un numero intero: ";
    std::cin >> n_rr;
  }
  int n_r = std::stoi(n_rr);
  std::string n_dd;
  std::cout << "Inserisci il numero di giorni di evoluzione della epidemia:";
  std::cin >> n_dd;
  while (is_integer(n_dd) == false) {
    std::cout << "L'input deve essere un numero intero: ";
    std::cin >> n_dd;
  }
  int n_d = std::stoi(n_dd);
  std::string betaa_;
  std::cout
      << "Inserisci il parametro beta (valore compreso tra 0 e 1 inclusi):";
  std::cin >> betaa_;
  while (is_number(betaa_) == false) {
    std::cout << "L'input deve essere un numero: ";
    std::cin >> betaa_;
  }
  double beta_ = std::stod(betaa_);
  std::string gammaa_;
  std::cout
      << "Inserisci il parametro gamma (valore compreso tra 0 e 1 inclusi):";

  std::cin >> gammaa_;
  while (is_number(gammaa_) == false) {
    std::cout << "L'input deve essere un numero: ";
    std::cin >> gammaa_;
  }
  double gamma_ = std::stod(gammaa_);
  std::string n_mm;
  std::cout << "Inserisci la percentuale di popolazione con la "
               "mascherina(numero fra 0 e 100):";

  std::cin >> n_mm;
  while (is_number(n_mm) == false) {
    std::cout << "L'input deve essere un numero: ";
    std::cin >> n_mm;
  }
  int n_m = std::stoi(n_mm);
  std::string n_vv;
  std::cout << "Inserisci la percentuale di popolazione con vaccinati(numero "
               "fra 0 e 100):";

  std::cin >> n_vv;
  while (is_number(n_vv) == false) {
    std::cout << "L'input deve essere un numero: ";
    std::cin >> n_vv;
  }
  int n_v = std::stoi(n_vv);
  std::string n_ee;
  std::cout << "Inserisci la l'efficacia del vaccino(numero fra 0 e 1):";

  std::cin >> n_ee;
  while (is_number(n_ee) == false) {
    std::cout << "L'input deve essere un numero: ";
    std::cin >> n_ee;
  }
  double n_e = std::stod(n_ee);

  try {
    Epidemic p(n_s, n_i, n_r, beta_, gamma_, n_d);
    p.Update();
    std::cout << std::setw(6) << "Giorno" << '\t' << std::setw(11)
              << "Suscettibili" << std::setw(11) << "Infetti" << '\t'
              << std::setw(11) << "Rimossi" << '\n';
    for (int i = 0; i <= n_d; ++i) {
      p.S_get(i);
      p.I_get(i);
      p.R_get(i);

      std::cout << std::setw(6) << i << '\t' << std::setw(11) << p.S_get(i)
                << std::setw(11) << p.I_get(i) << '\t' << std::setw(11)
                << p.R_get(i) << '\n';
    };
    std::cout << "\n\n";
    std::cout
        << "Nella finestra grafica che si apre i suscettibili sono verdi, "
           "gli infetti rossi e i rimossi blu."
        << '\n';
    Board b(n_s, n_i, n_r, beta_, gamma_, n_d, n_e, n_m, n_v);
    b.draw();
  } catch (std::runtime_error const& e) {
    std::cerr << "Caught exception: " << e.what() << '\n';
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "Caught unknown exception\n";
    return EXIT_FAILURE;
  }
}
