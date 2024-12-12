#include "compressible.hpp"

double Compressible::kappa;

Compressible (*Compressible::flux)(const Compressible& wl, const Compressible& wr, const Vector2d& s);

double Compressible::p() {
  return (kappa - 1.) * (e - 0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / rho);
}

double Compressible::a() {
  return sqrt(kappa * p() / rho);
}

double Compressible:Ma() {
  return rhoU.length() / rho / a();
}

