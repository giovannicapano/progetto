
#include "board.hpp"

#include <SFML/Graphics.hpp>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>
#include <thread>

#include "font.hpp"

inline std::default_random_engine generator;
inline std::uniform_real_distribution<double> distr(0.0, 1.0);
bool vaccination_clicked = false;
bool masking_clicked = false;

Board::Board(int s, int i, int r, double b, double y, int days, double e, int m,
             int v)
    : dimension{static_cast<int>(sqrt(s + i + r))},
      cell_number{static_cast<int>(pow(dimension, 2))},
      grid(cell_number),
      beta_{b},
      gamma_{y},
      count_time{0},
      count{s, i, r},
      days{days},
      effectiveness{e},
      mask_percentage{m},
      vaccine_percentage{v} {
  if (count.sus < 0) {
    throw std::runtime_error{
        "Il numero di suscettibili deve essere non negativo"};
  };
  if (count.inf < 0) {
    throw std::runtime_error{"Il numero di infetti deve essere non negativo"};
  }
  if (count.rim < 0) {
    throw std::runtime_error{"Il numero di rimossi deve essere non negativo"};
  }
  if (effectiveness < 0 || effectiveness > 1) {
    throw std::runtime_error{
        "L'efficacia del vaccino deve essere non negativa"};
  }
  if (vaccine_percentage < 0 || vaccine_percentage > 100) {
    throw std::runtime_error{
        "La percentuale di individui vaccinati deve essere non negativa e non "
        "maggiore di 100"};
  }
  if (mask_percentage < 0 || mask_percentage > 100) {
    throw std::runtime_error{
        "La percentuale di individui con la mascherina deve essere non "
        "negativa e non maggiore di 100"};
  }
  int S_M = static_cast<int>(y / b);
  t_max = 1 / (beta_) *
          ((i + s + y / b * log((i * s + (S_M - s) ^ 2) / (i * y / b))));
  c = 4 * pow(t_max, 2);
  a = pow(M_E, 1 / 4) / sqrt(t_max);

  int l = 0;
  while (l < count.inf) {
    int ran = distr(generator) *cell_number;
    if (grid[ran].diag != Diag::i) {
      grid[ran].diag = Diag::i;
      l++;
    }
  }


  int f = 0;
  while (f < count.rim) {
    int ran = distr(generator) * cell_number;
    if (grid[ran].diag != Diag::i || grid[ran].diag != Diag::r) {
      grid[ran].diag = Diag::r;
      f++;
    }
  }
}

bool Board::beta_random(int n, int t, bool b) const {
  double prob_raw = (beta_ / gamma_ - 1) * beta_ / (gamma_) +
                    (1 - ((1 - pow(1 - (beta_ * t), n)))) * a *
                        sqrt(count_time) * pow(M_E, -pow(count_time, 2) / c);
  double prob_vacc =
      (2 * beta_ / gamma_ - 1) * beta_ / (1 + gamma_) * (1 - effectiveness) +
      (1 - pow(1 - (beta_ * t * (1 - effectiveness)), n)) * a *
          sqrt(count_time) * (1 - effectiveness) *
          pow(M_E, -pow(count_time, 2) / c);
  if (b == false) {
    bool x(distr(generator) < prob_raw);
    return x;
  } else {
    bool y(distr(generator) < prob_vacc);
    return y;
  }
  assert((prob_vacc < prob_raw));
}

bool Board::gamma_random() const {
  bool y(distr(generator) < gamma_);
  return y;
}
void Board::mask_distribution() {
  int l = 0;
  int perc = mask_percentage * 0.01 * (count.sum());
  while (l < perc) {
    int ran = distr(generator) * cell_number;

    if (grid[ran].is_masked == false) {
      grid[ran].is_masked = true;
      l++;
    }
  }
}

void Board::mask_null() {
  for (auto& cell : grid) {
    {
      if (cell.is_masked == true) {
        cell.is_masked = false;
      }
    }
  }
}
void Board::vaccine_distribution() {
  int l = 0;
  int perc = (vaccine_percentage * 0.01) * (count.sum());
  while (l < perc) {
    int ran = distr(generator) * cell_number;

    if (grid[ran].is_vaccinated == false) {
      grid[ran].is_vaccinated = true;
      l++;
    }
  }
}

void Board::vaccine_null() {
  for (auto& cell : grid) {
    if (cell.is_vaccinated == true) {
      cell.is_vaccinated = false;
    }
  }
  
}

void Board::counting(const Diag& diag) {
  switch (diag) {
    case Diag::s:
      count.sus += 1;
      break;
    case Diag::i:
      count.inf += 1;
      break;
    case Diag::r:
      count.rim += 1;
      break;
  }
};
void Board::evolve() {
  count = {0, 0, 0};
  int count_i = 0;
  std::vector<Person> virtual_grid = grid;
  int ncell = -1;
  for (auto& cell : virtual_grid) {
    ncell += 1;
    switch (cell.diag) {
      case Diag::s:
        if (cell.is_masked == true) {
          cell.diag = Diag::s;
        } else {
          for (int i : {-1, 0, 1}) {
            for (int j : {-1, 0, 1}) {
              if (ncell + dimension * j + i >= 0 &&
                  ncell + dimension * j + i < cell_number) {
                if (!(i == 0 && j == 0)) {
                  if (virtual_grid[ncell + dimension * j + i].diag ==
                          Diag::i &&
                      virtual_grid[ncell + dimension * j + i].is_masked ==
                          false) {
                    count_i += 1;
                  }
                }
              }
            }
          }

          if (count_i >= 1) {
            cell.contact_time += 1;
            if (beta_random(count_i, cell.contact_time, cell.is_vaccinated) ==
                true) {
              grid[ncell].diag = Diag::i;

            } else {
              grid[ncell].diag = Diag::s;
            }
          } else {
            grid[ncell].diag = Diag::s;
          }
        }
        break;
      case Diag::i:
        if (gamma_random() == true) {
          grid[ncell].diag = Diag::r;
        } else {
          grid[ncell].diag = Diag::i;
        }
        break;
      case Diag::r:
        grid[ncell].diag = Diag::r;
        break;
    }
    count_i = 0;
  }

  std::for_each(grid.begin(), grid.end(), [&](Person a) { counting(a.diag); });

  count_time += 1;
}

void Board::draw() {
  const float win_size = 1000.;
  float bit_size = win_size / (dimension);
  
  sf::RenderWindow window(sf::VideoMode(win_size * 2, win_size + win_size / 5),
                          "epidemia");
  sf::RectangleShape sus_bit(sf::Vector2f(bit_size, bit_size));
  sf::RectangleShape inf_bit(sf::Vector2f(bit_size, bit_size));
  sf::RectangleShape rec_bit(sf::Vector2f(bit_size, bit_size));

  sus_bit.setFillColor(sf::Color::Green);
  inf_bit.setFillColor(sf::Color::Red);
  rec_bit.setFillColor(sf::Color::Blue);
  
  sf::Vertex x_axis[2];
  x_axis[0].position = sf::Vector2f(0, win_size);
  x_axis[1].position = sf::Vector2f(win_size * 2, win_size);
  
  sf::Vertex y_axis[2];
  y_axis[0].position = sf::Vector2f(win_size, 0);
  y_axis[1].position = sf::Vector2f(win_size, win_size);
  y_axis[1].color = sf::Color::White;
  
  sf::CircleShape arrowy(10, 3);
  arrowy.setPosition(sf::Vector2f(win_size - 10, 4));
  arrowy.setFillColor(sf::Color::White);
  arrowy.setOutlineColor(sf::Color::Black);
  
  sf::CircleShape arrowx(10, 3);
  arrowx.setRotation(90.f);
  arrowx.setPosition(2 * win_size - 4, win_size - 10);
  arrowx.setFillColor(sf::Color::White);
  arrowx.setOutlineColor(sf::Color::Black);
  
  sf::Font font;
  font.loadFromMemory(chintzy_ttf, chintzy_ttf_len);
  
  sf::CircleShape marker[3];
  marker[0].setFillColor(sf::Color::Green);
  marker[1].setFillColor(sf::Color::Red);
  marker[2].setFillColor(sf::Color::Blue);
  
  sf::Text leg[3];
  leg[0].setString("SUSCEPTIBLE:");
  leg[1].setString("INFECTED:");
  leg[2].setString("REMOVED:");
  
  sf::VertexArray inf_graph(sf::LinesStrip);
  sf::VertexArray sus_graph(sf::LinesStrip);
  sf::VertexArray rim_graph(sf::LinesStrip);
  sf::Text num[3];
  sf::Text graph_label[2];
  sf::RectangleShape buttons[2];
  sf::Text switch_[2];
  sf::Text state[2];
  
  switch_[0].setPosition(sf::Vector2f(0.2 * win_size + 17.5, win_size + 60));
  switch_[1].setPosition(sf::Vector2f(0.5 * win_size + 17.5, win_size + 60));
  
  state[0].setString("MASKING");
  state[1].setString("VACCINATION");
  state[0].setPosition(sf::Vector2f(0.2 * win_size - 10, win_size + 120));
  state[1].setPosition(sf::Vector2f(0.5 * win_size - 20, win_size + 120));
  
  buttons[0].setPosition(sf::Vector2f(0.2 * win_size, win_size + 50));
  buttons[1].setPosition(sf::Vector2f(0.5 * win_size, win_size + 50));
  
  graph_label[0].setFont(font);
  graph_label[0].setString("TIME");
  graph_label[0].setPosition(sf::Vector2f(2 * win_size - win_size / 10 - 10,
                                          win_size + win_size / 120));
                                          graph_label[1].setFont(font);
  graph_label[1].setString("RELATIVE NUMER OF CASES");
  graph_label[1].setPosition(sf::Vector2f(win_size + win_size / 20, 20));
  graph_label[1].setRotation(90.f);
  
  for (int i = 0; i < 2; i++) {
    switch_[i].setFillColor(sf::Color::Black);
    switch_[i].setFont(font);
    switch_[i].setString("OFF");
    switch_[i].setCharacterSize(win_size / 40);
    state[i].setFillColor(sf::Color::White);
    state[i].setFont(font);
    state[i].setCharacterSize(win_size / 40);
    buttons[i].setFillColor(sf::Color::White);
    buttons[i].setSize(sf::Vector2f(100, 50));
    buttons[i].setOutlineColor(sf::Color(187, 179, 179));
    buttons[i].setOutlineThickness(10);
    graph_label[i].setFillColor(sf::Color::White);
    graph_label[i].setCharacterSize(win_size / 40);
    x_axis[i].color = sf::Color::White;
    y_axis[i].color = sf::Color::White;
  }
  for (int i = 0; i < 3; i++) {
    marker[i].setPosition(
        sf::Vector2f(win_size + win_size * 0.4 - 30, win_size + 75 * i + 10));
    marker[i].setRadius(10);
    leg[i].setPosition(
        sf::Vector2f(win_size + win_size * 0.4, win_size + 75 * i));
    leg[i].setFont(font);
    leg[i].setCharacterSize(win_size / 40);
    leg[i].setFillColor(sf::Color::White);
  }
  while (window.isOpen()) {
    sf::Event event;
    while (window.pollEvent(event)) {
      switch (event.type) {
        case sf::Event::Closed:
          window.close();
          break;
        case sf::Event::MouseButtonPressed: {
          sf::Vector2i mousePos = sf::Mouse::getPosition(window);
          sf::Vector2f mousePosF(mousePos.x, mousePos.y);
          if (vaccination_clicked == false) {
            if (buttons[1].getGlobalBounds().contains(mousePosF)) {
              buttons[1].setFillColor(sf::Color::Yellow);
              buttons[1].setOutlineColor(sf::Color(252, 190, 59));
              switch_[1].setString("ON");
              vaccine_distribution();
              vaccination_clicked = true;
            }
          } else {
            if (buttons[1].getGlobalBounds().contains(mousePosF)) {
              buttons[1].setFillColor(sf::Color::White);
              buttons[1].setOutlineColor(sf::Color(187, 179, 179));
              switch_[1].setString("OFF");
              vaccine_null();
              vaccination_clicked = false;
            }
          }
          if (masking_clicked == false) {
            if (buttons[0].getGlobalBounds().contains(mousePosF)) {
              buttons[0].setFillColor(sf::Color::Yellow);
              buttons[0].setOutlineColor(sf::Color(252, 190, 59));
              switch_[0].setString("ON");
              mask_distribution();
              masking_clicked = true;
            }
          } else {
            if (buttons[0].getGlobalBounds().contains(mousePosF)) {
              buttons[0].setFillColor(sf::Color::White);
              buttons[0].setOutlineColor(sf::Color(187, 179, 179));
              switch_[0].setString("OFF");
              mask_null();
              masking_clicked = false;
            }
          }
        }
        default:
          break;
      }
    }

    window.clear();
    window.setFramerateLimit(60);

    for (int i = 0; i < cell_number; ++i) {
      int x = i % dimension;

      int y = (i - x) / dimension;

      switch (grid[i].diag) {
        case Diag::s:
          sus_bit.setPosition(bit_size * x, bit_size * y);
          window.draw(sus_bit);
          break;
        case Diag::i:
          inf_bit.setPosition(bit_size * x, bit_size * y);
          window.draw(inf_bit);
          break;
        case Diag::r:
          rec_bit.setPosition(bit_size * x, bit_size * y);
          window.draw(rec_bit);
          break;
      }
    }

    std::string numi = std::to_string(count.inf);
    std::string nums = std::to_string(count.sus);
    std::string numr = std::to_string(count.rim);
    num[0].setString(nums);
    num[1].setString(numi);
    num[2].setString(numr);
    num[0].setFillColor(sf::Color::Green);
    num[1].setFillColor(sf::Color::Red);
    num[2].setFillColor(sf::Color::Blue);
    for (int i = 0; i != 3; i++) {
    }
    window.draw(x_axis, 2, sf::Lines);
    window.draw(y_axis, 2, sf::Lines);

    inf_graph.resize(count_time + 1);
    inf_graph[count_time].position = sf::Vector2f(
        win_size + win_size * count_time / days,
        win_size - win_size * count.inf / (count.sus + count.inf + count.rim));
    inf_graph[count_time].color = sf::Color::Red;
    sus_graph.resize(count_time + 1);
    sus_graph[count_time].position = sf::Vector2f(
        win_size + win_size * count_time / days,
        win_size - win_size * count.sus / (count.sus + count.inf + count.rim));
    sus_graph[count_time].color = sf::Color::Green;
    rim_graph.resize(count_time + 1);
    rim_graph[count_time].position = sf::Vector2f(
        win_size + win_size * count_time / days,
        win_size - win_size * count.rim / (count.sus + count.inf + count.rim));
    rim_graph[count_time].color = sf::Color::Blue;
    window.draw(inf_graph);
    window.draw(sus_graph);
    window.draw(rim_graph);
    window.draw(arrowx);
    window.draw(arrowy);
    for (int i = 0; i != 2; i++) {
      window.draw(buttons[i]);
      window.draw(switch_[i]);
      window.draw(state[i]);
      window.draw(graph_label[i]);
    }
    for (int i = 0; i != 3; i++) {
      window.draw(marker[i]);
      window.draw(num[i]);
      window.draw(leg[i]);
      num[i].setPosition(sf::Vector2f(
          win_size + win_size * 0.4 + leg[i].getGlobalBounds().width + 5,
          win_size + 75 * i));
      num[i].setFont(font);
      num[i].setCharacterSize(win_size / 40);
    }
    

    window.display();
    if (count_time == days + 1) {
      window.clear();
      sf::Text text_end;
      text_end.setFont(font);
      text_end.setFillColor(sf::Color::White);
      text_end.setString("SIMULATION ENDED:");
      text_end.setPosition(0.35 * win_size, 0.45 * win_size);
      text_end.setCharacterSize(0.1 * win_size);
      num[0].setPosition(1.5 * win_size, 0.25 * win_size);
      num[1].setPosition(1.5 * win_size, 0.5 * win_size);
      num[2].setPosition(1.5 * win_size, 0.75 * win_size);
      leg[0].setPosition(1.5 * win_size, 0.2 * win_size);
      leg[1].setPosition(1.5 * win_size, 0.45 * win_size);
      leg[2].setPosition(1.5 * win_size, 0.7 * win_size);
      for (int i = 0; i != 3; i++) {
        num[i].setCharacterSize(0.1 * win_size);
        leg[i].setCharacterSize(0.06 * win_size);
        window.draw(leg[i]);
        window.draw(num[i]);
      }
      window.draw(text_end);

      window.display();
      std::this_thread::sleep_for(std::chrono::seconds(7));
      window.close();
    }

    std::this_thread::sleep_for(std::chrono::seconds(1));
    evolve();
  }
}
