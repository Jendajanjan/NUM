#ifndef FIRSTDERIVATIVES_HPP
#define FIRSTDERIVATIVES_HPP

#include <cmath>
#include "geometry/vector.hpp"

template <typename var>
void firstDerivatives(Vector2<var>& gradW, const var& wA, const var& wB,
		      const var& wL, const var& wR, const Point2d& L,
		      const Point2d& R, const Vector2d& s) {

  Vector2d LR(L, R);

  Vector2d sLR(LR.y, -LR.x);

  double volume = fabs(cross(sLR, s)) / 2.;

  gradW.x = 1. / (2. * volume) * ((wR - wL) * s.x + (wA - wB) * sLR.x);
  gradW.y = 1. / (2. * volume) * ((wR - wL) * s.y + (wA - wB) * sLR.y);
}

#endif
