#ifndef _BS_H_
#define _BS_H_

#include <cassert>
#include <cmath>
#include <cstddef>

double normalCDF(double x) { return std::erfc(-x / std::sqrt(2)) / 2; }

double blackscholes(double S, double E, double r, double time, double sigma,
                    char ch = 'C') {
  double d1 = (log(S / E) + (r + 0.5 * sigma * sigma) * time) /
              (sigma * std::sqrt(time));
  double d2 = d1 - sigma * std::sqrt(time);

  double ans = 0.0;
  if (ch == 'C')
    ans = S * normalCDF(d1) - E * std::exp(-r * time) * normalCDF(d2);
  else if (ch == 'P')
    ans = E * std::exp(-r * time) * normalCDF(-d2) - S * normalCDF(-d1);

  return ans;
}

class transformed_bsm {
 public:
  transformed_bsm(double E, double r, double time, double sigma, char c)
      : E_(E),
        r_(r),
        time_(time),
        sigma_(sigma),
        k_(r / (0.5 * sigma * sigma)),
        ch_(c) {}

  double pay_off(double x) const {
    double ans = 0.0;
    if (ch_ == 'C') ans = exp(0.5 * (k_ + 1) * x) - exp(0.5 * (k_ - 1) * x);
    if (ch_ == 'P') ans = exp(0.5 * (k_ - 1) * x) - exp(0.5 * (k_ + 1) * x);

    return ans > 0 ? ans : 0;
  }

  double u_m_inf(double x, double tau) const { return u(x, tau); }

  double u_p_inf(double x, double tau) const { return u(x, tau); }

  // private:
  double u(double x, double tau) const {
    double tmp = sqrt(2 * tau);
    double d1 = x / tmp + 0.5 * tmp * (k_ + 1);
    double d2 = x / tmp + 0.5 * tmp * (k_ - 1);

    double mult1 = exp(0.5 * (k_ + 1) * x + 0.25 * (k_ + 1) * (k_ + 1) * tau);
    double mult2 = exp(0.5 * (k_ - 1) * x + 0.25 * (k_ - 1) * (k_ - 1) * tau);

    double ans = 0.0;
    if (ch_ == 'C') ans = mult1 * normalCDF(d1) - mult2 * normalCDF(d2);
    if (ch_ == 'P') ans = mult2 * normalCDF(-d2) - mult1 * normalCDF(-d1);

    return ans;
  }

 private:
  double E_;
  double r_;
  double time_;
  double sigma_;

  double k_;
  char ch_;
};

double transform_u_back_to_option_value(std::vector<double> &u, double dx,
                                        int Nminus, int Nplus, double S,
                                        double E, double r, double time,
                                        double sigma) {
  assert(log(S / E) <= Nplus * dx);

  double x = log(S / E);
  double k = r / (0.5 * sigma * sigma);
  int idx = round((x - Nminus * dx) / dx);

  return std::pow(E, 0.5 * (1 + k)) * std::pow(S, 0.5 * (1 - k)) *
         std::exp(-0.125 * (k + 1) * (k + 1) * sigma * sigma * time) * u[idx];
}

#endif
