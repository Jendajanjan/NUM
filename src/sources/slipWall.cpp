#include "slipWall.hpp"

Compressible slipWall(const Compressible& wInside, const Vector2d& s, const Setting& setting) {
  Vector2d n = s / s.length();

  Vector2d rhoU = wInside.rhoU - 2. * dot(wInside.rhoU, n) * n;

  return Compressible(wInside.rho, rhoU, wInside.e);
}
