#include "compressible.hpp"

double Compressible::kappa;

Compressible (*Compressible::flux)(const Compressible& wl, const Compressible& wr, const Vector2d& s);

double Compressible::p() const {
  return (kappa - 1.) * (e - 0.5 * (pow(rhoU.x, 2) + pow(rhoU.y, 2)) / rho);
}

double Compressible::a() const {
  return std::sqrt(kappa * p() / rho);
}

double Compressible::Ma() const {
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

Compressible Compressible::sqrt(const Compressible& a) {
  double _rho = std::sqrt(a.rho);
  double _rhou = std::sqrt(a.rhoU.x);
  double _rhov = std::sqrt(a.rhoU.y);
  double _e = std::sqrt(a.e);

  return Compressible(_rho, Vector2d(_rhou, _rhov), _e);
}

Compressible Compressible::Upwind(const Compressible& wl, const Compressible& wr, const Vector2d& s) {
  Vector2d n = s / s.length();

  Vector2d u = (wl.rhoU/wl.rho + wr.rhoU/wr.rho) / 2.;
  double un = dot(u, n);

  Compressible flx;

  if (un >= 0.) flx = wl * un;
  else flx = wr * un;

  double p = (wl.p() + wr.p()) / 2.;

  flx += Compressible(0., p*n, p*un);

  return flx * s.length();
}


