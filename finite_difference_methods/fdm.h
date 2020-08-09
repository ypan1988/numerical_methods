#ifndef _FDM_H_
#define _FDM_H_

#include <vector>

struct gindex {
  gindex(int hi, int lo) : hi_(hi), lo_(lo) {}
  int operator()(int i) {
    assert(i <= hi_);
    return i - lo_;
  }
  int hi_, lo_;
};

template <typename DiffusionEquation>
void explicit_fd(std::vector<double> &values, double dx, double dt, int M,
                 int Nminus, int Nplus, DiffusionEquation &de) {
  gindex idx(Nplus, Nminus);

  values.resize(Nplus - Nminus + 1, 0);
  std::vector<double> oldu(Nplus - Nminus + 1, 0);
  std::vector<double> newu(Nplus - Nminus + 1, 0);

  double a = dt / (dx * dx);

  for (int n = Nminus; n <= Nplus; ++n) {
    oldu[idx(n)] = de.pay_off(n * dx);
  }

  for (int m = 1; m <= M; ++m) {
    double tau = m * dt;

    newu[idx(Nminus)] = de.u_m_inf(Nminus * dx, tau);
    newu[idx(Nplus)] = de.u_p_inf(Nplus * dx, tau);

    for (int n = Nminus + 1; n < Nplus; ++n) {
      newu[idx(n)] = oldu[idx(n)] + a * (oldu[idx(n - 1)] - 2 * oldu[idx(n)] +
                                         oldu[idx(n + 1)]);
    }

    for (int n = Nminus; n <= Nplus; ++n) {
      oldu[idx(n)] = newu[idx(n)];
    }
  }

  for (int n = Nminus; n <= Nplus; ++n) {
    values[idx(n)] = oldu[idx(n)];
  }
}

int lu_find_y(std::vector<double> &y, double a, int Nminus, int Nplus) {
  gindex idx(Nplus, Nminus);

  y.resize(Nplus - Nminus + 1, 0);

  double asq = a * a;
  y[idx(Nminus + 1)] = 1 + 2 * a;

  for (int n = Nminus + 2; n < Nplus; ++n) {
    y[idx(n)] = 1 + 2 * a - asq / y[idx(n - 1)];
    if (y[idx(n)] == 0.0) return -1;
  }

  return 0;
}

/* must call lu_find_y before using this */
void lu_solver(std::vector<double> &u, std::vector<double> &b,
               std::vector<double> &y, double a, int Nminus, int Nplus) {
  gindex idx(Nplus, Nminus);

  u.resize(Nplus - Nminus + 1, 0);

  std::vector<double> q(Nplus - Nminus + 1, 0);
  q[idx(Nminus + 1)] = b[idx(Nminus + 1)];

  for (int n = Nminus + 2; n < Nplus; ++n) {
    q[idx(n)] = b[idx(n)] + a * q[idx(n - 1)] / y[idx(n - 1)];
  }

  u[idx(Nplus - 1)] = q[idx(Nplus - 1)] / y[idx(Nplus - 1)];

  for (int n = Nplus - 2; n > Nminus; --n) {
    u[idx(n)] = (q[idx(n)] + a * u[idx(n + 1)]) / y[idx(n)];
  }
}

template <typename DiffusionEquation>
void implicit_fd1(std::vector<double> &values, double dx, double dt, int M,
                  int Nminus, int Nplus, DiffusionEquation &de) {
  gindex idx(Nplus, Nminus);

  values.resize(Nplus - Nminus + 1, 0);

  double a = dt / (dx * dx);

  for (int n = Nminus; n <= Nplus; ++n) {
    values[idx(n)] = de.pay_off(n * dx);
  }

  std::vector<double> y;
  lu_find_y(y, a, Nminus, Nplus);

  std::vector<double> b(Nplus - Nminus + 1, 0);
  for (int m = 1; m <= M; ++m) {
    double tau = m * dt;

    for (int n = Nminus + 1; n < Nplus; ++n) {
      b[idx(n)] = values[idx(n)];
    }

    values[idx(Nminus)] = de.u_m_inf(Nminus * dx, tau);
    values[idx(Nplus)] = de.u_p_inf(Nplus * dx, tau);
    b[idx(Nminus + 1)] += a * values[idx(Nminus)];
    b[idx(Nplus - 1)] += a * values[idx(Nplus)];

    lu_solver(values, b, y, a, Nminus, Nplus);
  }
}

int SOR_solver(std::vector<double> &u, std::vector<double> &b, int Nminus,
               int Nplus, double a, double omega, double eps) {
  gindex idx(Nplus, Nminus);

  double error, y;

  int loops = 0;
  do {
    error = 0.0;
    for (int n = Nminus + 1; n < Nplus; ++n) {
      y = (b[idx(n)] + a * (u[idx(n - 1)] + u[idx(n + 1)])) / (1 + 2 * a);
      y = u[idx(n)] + omega * (y - u[idx(n)]);
      error += (u[idx(n)] - y) * (u[idx(n)] - y);
      u[idx(n)] = y;
    }
    ++loops;
  } while (error > eps);

  return loops;
}

template <typename DiffusionEquation>
void implicit_fd2(std::vector<double> &values, double dx, double dt, int M,
                  int Nminus, int Nplus, DiffusionEquation &de) {
  gindex idx(Nplus, Nminus);

  values.resize(Nplus - Nminus + 1, 0);

  double a = dt / (dx * dx);
  double eps = 1.0e-8;
  double omega = 1.0;
  double domega = 0.05;
  int oldloops = 10000;

  for (int n = Nminus; n <= Nplus; ++n) {
    values[idx(n)] = de.pay_off(n * dx);
  }

  std::vector<double> b(Nplus - Nminus + 1, 0);
  for (int m = 1; m <= M; ++m) {
    double tau = m * dt;

    for (int n = Nminus + 1; n < Nplus; ++n) {
      b[idx(n)] = values[idx(n)];
    }

    values[idx(Nminus)] = de.u_m_inf(Nminus * dx, tau);
    values[idx(Nplus)] = de.u_p_inf(Nplus * dx, tau);

    int loops = SOR_solver(values, b, Nminus, Nplus, a, omega, eps);
    if (loops > oldloops) domega *= -1.0;
    omega += domega;
    oldloops = loops;
  }
}

template <typename DiffusionEquation>
void crank_fd1(std::vector<double> &val, double dx, double dt, int M,
               int Nminus, int Nplus, DiffusionEquation &de) {
  gindex idx(Nplus, Nminus);

  val.resize(Nplus - Nminus + 1, 0);

  double a = dt / (dx * dx);
  double a2 = a / 2.0;

  for (int n = Nminus; n <= Nplus; ++n) {
    val[idx(n)] = de.pay_off(n * dx);
  }

  std::vector<double> y;
  lu_find_y(y, a2, Nminus, Nplus);

  std::vector<double> b(Nplus - Nminus + 1, 0);
  for (int m = 1; m <= M; ++m) {
    double tau = m * dt;

    for (int n = Nminus + 1; n < Nplus; ++n) {
      b[idx(n)] =
          (1 - a) * val[idx(n)] + a2 * (val[idx(n + 1)] + val[idx(n - 1)]);
    }

    val[idx(Nminus)] = de.u_m_inf(Nminus * dx, tau);
    val[idx(Nplus)] = de.u_p_inf(Nplus * dx, tau);
    b[idx(Nminus + 1)] += a2 * val[idx(Nminus)];
    b[idx(Nplus - 1)] += a2 * val[idx(Nplus)];

    lu_solver(val, b, y, a2, Nminus, Nplus);
  }
}

template <typename DiffusionEquation>
void crank_fd2(std::vector<double> &val, double dx, double dt, int M,
               int Nminus, int Nplus, DiffusionEquation &de) {
  gindex idx(Nplus, Nminus);

  val.resize(Nplus - Nminus + 1, 0);

  double a = dt / (dx * dx);
  double a2 = a / 2.0;
  double eps = 1.0e-8;
  double omega = 1.0;
  double domega = 0.05;
  int oldloops = 10000;

  for (int n = Nminus; n <= Nplus; ++n) {
    val[idx(n)] = de.pay_off(n * dx);
  }

  std::vector<double> b(Nplus - Nminus + 1, 0);
  for (int m = 1; m <= M; ++m) {
    double tau = m * dt;

    for (int n = Nminus + 1; n < Nplus; ++n) {
      b[idx(n)] =
          (1 - a) * val[idx(n)] + a2 * (val[idx(n + 1)] + val[idx(n - 1)]);
    }

    val[idx(Nminus)] = de.u_m_inf(Nminus * dx, tau);
    val[idx(Nplus)] = de.u_p_inf(Nplus * dx, tau);

    int loops = SOR_solver(val, b, Nminus, Nplus, a, omega, eps);
    if (loops > oldloops) domega *= -1.0;
    omega += domega;
    oldloops = loops;
  }
}

#endif
