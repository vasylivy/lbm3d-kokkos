#ifndef _TYPEDEF_H_
#define _TYPEDEF_H_

#include <Kokkos_Core.hpp>

typedef double Double;

typedef Kokkos::View<Double***> ScalarField;
typedef Kokkos::View<Double***[19]> DistributionField;

typedef Kokkos::MDRangePolicy< Kokkos::Rank<3> > range_3d;
typedef Kokkos::MDRangePolicy< Kokkos::Rank<2> > range_2d;

#endif

