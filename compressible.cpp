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

Compressible Compressible::max(const Compressible& a, const Compressible& b) {
  double _rho = std::max(a.rho, b.rho);
  double _rhou = std::max(a.rhoU.x, b.rhoU.x);
  double _rhov = std::max(a.rhoU.y, b.rhoU.y);
  double _e = std::max(a.e, b.e);

  return Compressible(_rho, Vector2d(_rhou, _rhov), _e);
}

Compressible Compressible::min(const Compressible& a, const Compressible& b) {
  double _rho = std::min(a.rho, b.rho);
  double _rhou = std::min(a.rhoU.x, b.rhoU.x);
  double _rhov = std::min(a.rhoU.y, b.rhoU.y);
  double _e = std::min(a.e, b.e);

  return Compressible(_rho, Vector2d(_rhou, _rhov), _e);
}

Compressible Compressible::fabs(const Compressible& a) {
  double _rho = std::fabs(a.rho);
  double _rhou = std::fabs(a.rhoU.x);
  double _rhov = std::fabs(a.rhoU.y);
  double _e = std::fabs(a.e);

  return Compressible(_rho, Vector2d(_rhou, _rhov), _e);
}

