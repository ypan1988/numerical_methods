#include <iostream>
#include <vector>

using namespace std;

#include "bs.h"
#include "fdm.h"

int main() {
  const vector<double> Svec = {0.01,  2.00,  4.00,  6.00,  7.00,  8.00,  9.00,
                               10.00, 11.00, 12.00, 13.00, 14.00, 15.00, 16.00};

  const double E = 10.0;
  const double r = 0.05;
  const double time = 0.5;
  const double sigma = 0.20;

  transformed_bsm de(E, r, time, sigma, 'P');

  const vector<double> alpha = {0.50, 1.00, 5.00};

  const int Nplus = 700;
  const int Nminus = -7000;
  const double dx = 0.01;

  double dt = 0.0;
  int M = 0;
  vector<double> dvec1, dvec2, dvec3;

  cout << "alpha = 0.50" << endl;
  dt = alpha[0] * dx * dx;
  M = 0.5 * sigma * sigma * time / dt;
  implicit_fd1(dvec1, dx, dt, M, Nminus, Nplus, de);

  cout << "alpha = 1.00" << endl;
  dt = alpha[1] * dx * dx;
  M = 0.5 * sigma * sigma * time / dt;
  implicit_fd1(dvec2, dx, dt, M, Nminus, Nplus, de);

  cout << "alpha = 5.00" << endl;
  dt = alpha[2] * dx * dx;
  M = 0.5 * sigma * sigma * time / dt;
  implicit_fd1(dvec3, dx, dt, M, Nminus, Nplus, de);

  for (auto S : Svec) {
    double s25 = transform_u_back_to_option_value(dvec1, dx, Nminus, Nplus, S,
                                                  E, r, time, sigma);
    double s50 = transform_u_back_to_option_value(dvec2, dx, Nminus, Nplus, S,
                                                  E, r, time, sigma);
    double s52 = transform_u_back_to_option_value(dvec3, dx, Nminus, Nplus, S,
                                                  E, r, time, sigma);
    double exact = blackscholes(S, E, r, time, sigma, 'P');

    printf("%8.2f %10.4f %10.4f %10.4f %10.4f \n", S, s25, s50, s52, exact);
  }

  return 0;
}
