#include <typedefs.h>
#include <params.h>
#include <math.h>

template <typename T>
KOKKOS_INLINE_FUNCTION T max(T x1, T x2) {
  return (x1 > x2) ? x1 : x2;
}

struct collide_stream {

  const Params params;
  const Double tau = params.tau;
  const Double tau_inv = 1. / params.tau;
  const Double wo = 1. / 3.;
  const Double ws = 1. / 18.;
  const Double wd = 1. / 36.;
  const Double omtau_inv = 1.0 - tau_inv;

  const DistributionField fA, fB;

  collide_stream(DistributionField fB, DistributionField fA, Params params_) :
      fB(fB), fA(fA), params(params_) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {
    // load distributions
    Double f0 = fA(i, j, k, 0);
    Double f1 = fA(i, j, k, 1);
    Double f2 = fA(i, j, k, 2);
    Double f3 = fA(i, j, k, 3);
    Double f4 = fA(i, j, k, 4);
    Double f5 = fA(i, j, k, 5);
    Double f6 = fA(i, j, k, 6);
    Double f7 = fA(i, j, k, 7);
    Double f8 = fA(i, j, k, 8);
    Double f9 = fA(i, j, k, 9);
    Double f10 = fA(i, j, k, 10);
    Double f11 = fA(i, j, k, 11);
    Double f12 = fA(i, j, k, 12);
    Double f13 = fA(i, j, k, 13);
    Double f14 = fA(i, j, k, 14);
    Double f15 = fA(i, j, k, 15);
    Double f16 = fA(i, j, k, 16);
    Double f17 = fA(i, j, k, 17);
    Double f18 = fA(i, j, k, 18);

    // compute density
    Double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
    Double density_inv = 1. / density;

    // compute velocities
    Double u = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
    Double v = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
    Double w = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;

    // update distribubtions
    Double tw0r = wo * density * tau_inv; // w[0]*rho
    Double twsr = ws * density * tau_inv; // w[1-6]*rho
    Double twdr = wd * density * tau_inv; // w[7-18]*rho

    Double tu = 3.0 * u;
    Double tv = 3.0 * v;
    Double tw = 3.0 * w;

    Double a = 1.0 - 1.5 * (u * u + v * v + w * w);

    Double udot7 = tu + tv;
    Double udot9 = tu + tw;
    Double udot11 = tv + tw;
    Double udot13 = tu - tv;
    Double udot15 = tu - tw;
    Double udot17 = tv - tw;

    Double feq0 = tw0r * a;
    Double feq1 = twsr * (a + tu * (1.0 + 0.5 * tu));
    Double feq2 = twsr * (a - tu * (1.0 - 0.5 * tu));
    Double feq3 = twsr * (a + tv * (1.0 + 0.5 * tv));
    Double feq4 = twsr * (a - tv * (1.0 - 0.5 * tv));
    Double feq5 = twsr * (a + tw * (1.0 + 0.5 * tw));
    Double feq6 = twsr * (a - tw * (1.0 - 0.5 * tw));
    Double feq7 = twdr * (a + udot7 * (1.0 + 0.5 * udot7));
    Double feq8 = twdr * (a - udot7 * (1.0 - 0.5 * udot7));
    Double feq9 = twdr * (a + udot9 * (1.0 + 0.5 * udot9));
    Double feq10 = twdr * (a - udot9 * (1.0 - 0.5 * udot9));
    Double feq11 = twdr * (a + udot11 * (1.0 + 0.5 * udot11));
    Double feq12 = twdr * (a - udot11 * (1.0 - 0.5 * udot11));
    Double feq13 = twdr * (a + udot13 * (1.0 + 0.5 * udot13));
    Double feq14 = twdr * (a - udot13 * (1.0 - 0.5 * udot13));
    Double feq15 = twdr * (a + udot15 * (1.0 + 0.5 * udot15));
    Double feq16 = twdr * (a - udot15 * (1.0 - 0.5 * udot15));
    Double feq17 = twdr * (a + udot17 * (1.0 + 0.5 * udot17));
    Double feq18 = twdr * (a - udot17 * (1.0 - 0.5 * udot17));

    Double fB0 = omtau_inv * f0 + feq0;
    Double fB1 = omtau_inv * f1 + feq1;
    Double fB2 = omtau_inv * f2 + feq2;
    Double fB3 = omtau_inv * f3 + feq3;
    Double fB4 = omtau_inv * f4 + feq4;
    Double fB5 = omtau_inv * f5 + feq5;
    Double fB6 = omtau_inv * f6 + feq6;
    Double fB7 = omtau_inv * f7 + feq7;
    Double fB8 = omtau_inv * f8 + feq8;
    Double fB9 = omtau_inv * f9 + feq9;
    Double fB10 = omtau_inv * f10 + feq10;
    Double fB11 = omtau_inv * f11 + feq11;
    Double fB12 = omtau_inv * f12 + feq12;
    Double fB13 = omtau_inv * f13 + feq13;
    Double fB14 = omtau_inv * f14 + feq14;
    Double fB15 = omtau_inv * f15 + feq15;
    Double fB16 = omtau_inv * f16 + feq16;
    Double fB17 = omtau_inv * f17 + feq17;
    Double fB18 = omtau_inv * f18 + feq18;

    // stream updated distributions
    fB(i, j, k, 0) = fB0;
    fB(i, j + 1, k, 1) = fB1;
    fB(i, j - 1, k, 2) = fB2;
    fB(i + 1, j, k, 3) = fB3;
    fB(i - 1, j, k, 4) = fB4;
    fB(i, j, k + 1, 5) = fB5;
    fB(i, j, k - 1, 6) = fB6;
    fB(i + 1, j + 1, k, 7) = fB7;
    fB(i - 1, j - 1, k, 8) = fB8;
    fB(i, j + 1, k + 1, 9) = fB9;
    fB(i, j - 1, k - 1, 10) = fB10;
    fB(i + 1, j, k + 1, 11) = fB11;
    fB(i - 1, j, k - 1, 12) = fB12;
    fB(i - 1, j + 1, k, 13) = fB13;
    fB(i + 1, j - 1, k, 14) = fB14;
    fB(i, j + 1, k - 1, 15) = fB15;
    fB(i, j - 1, k + 1, 16) = fB16;
    fB(i + 1, j, k - 1, 17) = fB17;
    fB(i - 1, j, k + 1, 18) = fB18;
  }
};

struct bb_left {

  const int j = 1;
  const DistributionField fB;

  bb_left(DistributionField fB) :
      fB(fB) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int k) const {
    // 2, 8, 10, 14, 16
    fB(i, j, k, 1) = fB(i, j - 1, k, 2);
    fB(i, j, k, 7) = fB(i - 1, j - 1, k, 8);
    fB(i, j, k, 9) = fB(i, j - 1, k - 1, 10);
    fB(i, j, k, 13) = fB(i + 1, j - 1, k, 14);
    fB(i, j, k, 15) = fB(i, j - 1, k + 1, 16);

  }
};

struct bb_right {

  const int j;
  const DistributionField fB;

  bb_right(DistributionField fB) :
      fB(fB), j(fB.extent(1) - 2) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int k) const {

    // 1, 7, 9, 13, 15
    fB(i, j, k, 2) = fB(i, j + 1, k, 1);
    fB(i, j, k, 8) = fB(i + 1, j + 1, k, 7);
    fB(i, j, k, 10) = fB(i, j + 1, k + 1, 9);
    fB(i, j, k, 14) = fB(i - 1, j + 1, k, 13);
    fB(i, j, k, 16) = fB(i, j + 1, k - 1, 15);

  }
};

struct bb_bottom {

  const int i = 1;

  const DistributionField fB;

  bb_bottom(DistributionField fB) :
      fB(fB) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, const int k) const {

    // 4, 8 , 12, 13, 18
    fB(i, j, k, 3) = fB(i - 1, j, k, 4);
    fB(i, j, k, 7) = fB(i - 1, j - 1, k, 8);
    fB(i, j, k, 11) = fB(i - 1, j, k - 1, 12);
    fB(i, j, k, 14) = fB(i - 1, j + 1, k, 13);
    fB(i, j, k, 17) = fB(i - 1, j, k + 1, 18);

  }
};

struct bb_top {

  const Double wd = 1. / 36;
  const Double u;
  const Double w;
  const int i;
  const DistributionField fB;

  bb_top(DistributionField fB, Params params_) :
      fB(fB), u(params_.lid_u), w(params_.lid_w), i(fB.extent(0) - 2) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int j, const int k) const {

    Double rhs_7 = 6.0 * wd * u;
    Double rhs_11 = 6.0 * wd * w;
    Double rhs_14 = -6.0 * wd * u;
    Double rhs_17 = -6.0 * wd * w;

    // 3, 7, 11, 14, 17
    fB(i, j, k, 4) = fB(i + 1, j, k, 3);
    fB(i, j, k, 8) = fB(i + 1, j + 1, k, 7) - rhs_7;
    fB(i, j, k, 12) = fB(i + 1, j, k + 1, 11) - rhs_11;
    fB(i, j, k, 13) = fB(i + 1, j - 1, k, 14) - rhs_14;
    fB(i, j, k, 18) = fB(i + 1, j, k - 1, 17) - rhs_17;

  }
};

struct bb_front {

  const int k;
  const DistributionField fB;

  bb_front(DistributionField fB) :
      fB(fB), k(fB.extent(2) - 2) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    // 5, 9, 11, 16, 18
    fB(i, j, k, 6) = fB(i, j, k + 1, 5);
    fB(i, j, k, 10) = fB(i, j + 1, k + 1, 9);
    fB(i, j, k, 12) = fB(i + 1, j, k + 1, 11);
    fB(i, j, k, 15) = fB(i, j - 1, k + 1, 16);
    fB(i, j, k, 17) = fB(i - 1, j, k + 1, 18);

  }
};

struct bb_back {

  const int k = 1;
  const DistributionField fB;

  bb_back(DistributionField fB) :
      fB(fB) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j) const {

    // 6,10,12,15,17
    fB(i, j, k, 5) = fB(i, j, k - 1, 6);
    fB(i, j, k, 9) = fB(i, j - 1, k - 1, 10);
    fB(i, j, k, 11) = fB(i - 1, j, k - 1, 12);
    fB(i, j, k, 16) = fB(i, j + 1, k - 1, 15);
    fB(i, j, k, 18) = fB(i + 1, j, k - 1, 17);

  }
};

struct compute_macroscopic {

  const ScalarField u, v, w, rho;
  const DistributionField f;

  compute_macroscopic(DistributionField f, ScalarField u, ScalarField v, ScalarField w, ScalarField rho) :
      f(f), u(u), v(v), w(w), rho(rho) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    // load distributions
    Double f0 = f(i, j, k, 0);
    Double f1 = f(i, j, k, 1);
    Double f2 = f(i, j, k, 2);
    Double f3 = f(i, j, k, 3);
    Double f4 = f(i, j, k, 4);
    Double f5 = f(i, j, k, 5);
    Double f6 = f(i, j, k, 6);
    Double f7 = f(i, j, k, 7);
    Double f8 = f(i, j, k, 8);
    Double f9 = f(i, j, k, 9);
    Double f10 = f(i, j, k, 10);
    Double f11 = f(i, j, k, 11);
    Double f12 = f(i, j, k, 12);
    Double f13 = f(i, j, k, 13);
    Double f14 = f(i, j, k, 14);
    Double f15 = f(i, j, k, 15);
    Double f16 = f(i, j, k, 16);
    Double f17 = f(i, j, k, 17);
    Double f18 = f(i, j, k, 18);

    // compute density
    Double density = f0 + f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10 + f11 + f12 + f13 + f14 + f15 + f16 + f17 + f18;
    Double density_inv = 1. / density;

    // compute velocities
    Double utmp = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16)) * density_inv;
    Double vtmp = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18)) * density_inv;
    Double wtmp = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17)) * density_inv;

    // write to device
    rho(i, j, k) = density;
    u(i, j, k) = utmp;
    v(i, j, k) = vtmp;
    w(i, j, k) = wtmp;

  }
};

struct load_state {

  const ScalarField u, v, w, rho;
  const DistributionField fA, fB;

  const Double weight[19] = { 1. / 3., 1. / 18, 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 18., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1.
      / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36., 1. / 36. };

  load_state(DistributionField fA, DistributionField fB, ScalarField u, ScalarField v, ScalarField w, ScalarField rho) :
      fA(fA), fB(fB), u(u), v(v), w(w), rho(rho) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k) const {

    rho(i, j, k) = 1.0;
    u(i, j, k) = 0.0;
    v(i, j, k) = 0.0;
    w(i, j, k) = 0.0;

    for (int z = 0; z < fA.extent(3); ++z) {

      Double tmp = weight[z];
      fA(i, j, k, z) = tmp;
      fB(i, j, k, z) = tmp;
    }
  }
};

struct is_steady_state {

  const Double tol;
  const DistributionField fA;
  const ScalarField u, v, w, rho;

  is_steady_state(DistributionField fA, ScalarField u, ScalarField v, ScalarField w, ScalarField rho, const Double tol) :
      fA(fA), u(u), v(v), w(w), rho(rho), tol(tol) {
  }

  KOKKOS_INLINE_FUNCTION
  void operator()(const int i, const int j, const int k, int &converged) const {

    // load distributions A
    Double f0 = fA(i, j, k, 0);
    Double f1 = fA(i, j, k, 1);
    Double f2 = fA(i, j, k, 2);
    Double f3 = fA(i, j, k, 3);
    Double f4 = fA(i, j, k, 4);
    Double f5 = fA(i, j, k, 5);
    Double f6 = fA(i, j, k, 6);
    Double f7 = fA(i, j, k, 7);
    Double f8 = fA(i, j, k, 8);
    Double f9 = fA(i, j, k, 9);
    Double f10 = fA(i, j, k, 10);
    Double f11 = fA(i, j, k, 11);
    Double f12 = fA(i, j, k, 12);
    Double f13 = fA(i, j, k, 13);
    Double f14 = fA(i, j, k, 14);
    Double f15 = fA(i, j, k, 15);
    Double f16 = fA(i, j, k, 16);
    Double f17 = fA(i, j, k, 17);
    Double f18 = fA(i, j, k, 18);

    // compute velocities A
    Double umom_A = (f1 + f7 + f9 + f13 + f15 - (f2 + f8 + f10 + f14 + f16));
    Double vmom_A = (f3 + f7 + f11 + f14 + f17 - (f4 + f8 + f12 + f13 + f18));
    Double wmom_A = (f5 + f9 + f11 + f16 + f18 - (f6 + f10 + f12 + f15 + f17));

    // load macroscopic properties B
    Double rho_B = rho(i, j, k);
    Double umom_B = u(i, j, k) * rho_B;
    Double vmom_B = v(i, j, k) * rho_B;
    Double wmom_B = w(i, j, k) * rho_B;

    // convergence criteria
    int converged_momu = fabs(umom_B - umom_A) < max(tol * fabs(umom_A), 1e-12);
    int converged_momv = fabs(vmom_B - vmom_A) < max(tol * fabs(vmom_A), 1e-12);
    int converged_momw = fabs(wmom_B - wmom_A) < max(tol * fabs(wmom_A), 1e-12);
    converged &= converged_momu & converged_momv & converged_momw;
  }
};

