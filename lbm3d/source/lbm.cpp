#include <ostream>
#include <stdlib.h>
#include <memory>

#include "functors.h"
#include "params.h"
#include "output.h"

void update(DistributionField fB, DistributionField fA, ScalarField u, ScalarField v, ScalarField w, ScalarField rho, Params &params, const int step,
    int &converged);

int main(int narg, char *arg[]) {
  Kokkos::initialize(narg, arg);

  Params params;

  std::unique_ptr < Output > o(new Output());

  {
    // params
    params.nx = atoi(arg[1]);
    params.ny = atoi(arg[2]);
    params.nz = atoi(arg[3]);
    params.re = atof(arg[4]);
    params.tol = atof(arg[5]);
    params.max_steps = atoi(arg[6]);
    params.output_rate = atoi(arg[7]);

    params.lid_u = 0.05;
    params.lid_w = 0.0;
    params.lid_mag = sqrt(params.lid_u * params.lid_u + params.lid_w * params.lid_w);
    params.nu = params.lid_mag * double(params.ny - 2) / params.re;
    params.tau = 3. * params.nu + 0.5;

    const int nx = params.nx;
    const int ny = params.ny;
    const int nz = params.nz;

    // allocate memory on device
    DistributionField fA("fA", ny, nx, nz);
    DistributionField fB("fB", ny, nx, nz);
    DistributionField fC = fA;

    ScalarField rho("rho", ny, nx, nz);
    ScalarField u("u", ny, nx, nz);
    ScalarField v("v", ny, nx, nz);
    ScalarField w("w", ny, nx, nz);

    // initialize values on device
    Kokkos::parallel_for("load_state", range_3d( { 0, 0, 0 }, { fA.extent(0), fA.extent(1), fA.extent(2) }), load_state(fA, fB, u, v, w, rho));

    ScalarField::HostMirror h_u = Kokkos::create_mirror_view(u);
    ScalarField::HostMirror h_v = Kokkos::create_mirror_view(v);
    ScalarField::HostMirror h_w = Kokkos::create_mirror_view(w);
    ScalarField::HostMirror h_rho = Kokkos::create_mirror_view(rho);

    // write t = 0 to output file
    o->write_view("output/u", h_u);
    o->write_view("output/v", h_v);
    o->write_view("output/w", h_w);
    o->write_view("output/rho", h_rho);
    o->frame += 1;

    printf("Solving lid driven cavity (Re = %.2e, tau = %.2e, [%i x %i x %i] , RAM (MB) = %.1f...\n", params.re, params.tau, nx - 2, ny - 2, nz - 2,
        double(nx) * double(ny) * double(nz) * 8. * (2. * 19. + 3.) / (1024. * 1024.));

    int converged = 0;
    int step = 0;

    Kokkos::Timer timer;

    while (step < params.max_steps && !converged) {

      // collide stream
      update(fB, fA, u, v, w, rho, params, step, converged);

      // swap pointers
      fC = fA;
      fA = fB;
      fB = fC;

      // write to file
      if ((step + 1) % params.output_rate == 0 || converged) {

        printf("...output step = %i\n", step + 1);

        // deep copy from device to host
        Kokkos::deep_copy(h_u, u);
        Kokkos::deep_copy(h_v, v);
        Kokkos::deep_copy(h_w, w);
        Kokkos::deep_copy(h_rho, rho);

        // write to output file
        o->write_view("output/u", h_u);
        o->write_view("output/v", h_v);
        o->write_view("output/w", h_w);
        o->write_view("output/rho", h_rho);

        o->frame += 1;
      }
      step += 1;
    }

    double time = timer.seconds();

    double site_updates = double(nx - 2) * double(ny - 2) * double(nz - 2) * double(step) / (1000. * 1000.);

    double msus = site_updates / time;

    double bandwidth = msus * 1000. * 1000. * 2. * 19. * 8. / (1024. * 1024. * 1024.);

    double cost = time / double(params.max_steps);

    if (converged) {
      printf("Solution converged to steady state tolerance of %.3e\n", params.tol);
    } else {
      printf("Solution did not converged within %i steps\n", params.max_steps);
    }

    printf("MLUPS = %.1f, GB/s = %.1f, cost per step (s) = %.3e\n", msus, bandwidth, cost);
  }

  Kokkos::finalize();

  exit(EXIT_SUCCESS);
}

void update(DistributionField fB, DistributionField fA, ScalarField u, ScalarField v, ScalarField w, ScalarField rho, Params &params, const int step,
    int &converged) {

  const int D1 = fA.extent(0);
  const int D2 = fA.extent(1);
  const int D3 = fA.extent(2);

  Kokkos::parallel_for("collide_stream", range_3d( { 1, 1, 1 }, { D1 - 1, D2 - 1, D3 - 1 }), collide_stream(fB, fA, params));

  Kokkos::parallel_for("bb_left", range_2d( { 1, 1 }, { D1 - 1, D3 - 1 }), bb_left(fB));

  Kokkos::parallel_for("bb_right", range_2d( { 1, 1 }, { D1 - 1, D3 - 1 }), bb_right(fB));

  Kokkos::parallel_for("bb_front", range_2d( { 1, 1 }, { D1 - 1, D2 - 1 }), bb_front(fB));

  Kokkos::parallel_for("bb_back", range_2d( { 1, 1 }, { D1 - 1, D2 - 1 }), bb_back(fB));

  Kokkos::parallel_for("bb_bottom", range_2d( { 1, 1 }, { D2 - 1, D3 - 1 }), bb_bottom(fB));

  Kokkos::parallel_for("bb_top", range_2d( { 1, 1 }, { D2 - 1, D3 - 1 }), bb_top(fB, params));

  if ((step + 1) % params.output_rate == 0) {

    Kokkos::parallel_for("compute_macroscopic", range_3d( { 1, 1, 1 }, { D1 - 1, D2 - 1, D3 - 1 }), compute_macroscopic(fB, u, v, w, rho));

    Kokkos::parallel_reduce("check_if_steady_state", range_3d( { 1, 1, 1 }, { D1 - 1, D2 - 1, D3 - 1 }), is_steady_state(fA, u, v, w, rho, params.tol),
        Kokkos::BAnd<int>(converged));
  }
}
