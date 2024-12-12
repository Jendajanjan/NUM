#ifndef COMPRESSIBLE_HPP
#define COMPRESSIBLE_HPP

#include <cmath>
#include <algorithm>
#include "geometry/vector.hpp"

using namespace std;

class Compressible {
public:
  double rho;
  Vector2d rhoU;
  double e;

  Compressible() {};
  Compressible(const double& rho_, const Vector2d& rhoU_, const double& e_):
    rho(rho_), rhoU(rhoU_), e(e_) {};
  ~Compressible() {};

  void zero() {
    rho = 0.; rhoU = Vector2d(0., 0.), e = 0.;
  }

  void ones() {
    rho = 1.; rhoU = Vector2d(1., 1.), e = 1.;
  }

  double p();
  double a();
  double Ma();

  static double kappa;

  static Compressible (*flux)(const Compressible& wl, const Compressible& wr, const Vector2d s);
};

inline Compressible operator+(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho+b.rho, a.rhoU+b.rhoU, a.e+b.e);
}

inline Compressible operator-(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho-b.rho, a.rhoU-b.rhoU, a.e-b.e);
}

inline Compressible operator*(const Compressible& a, const Compressible& b) {
  return Compressible(a.rho*b.rho, a.rhoU*b.rhoU, a.e*b.e);
}

template <typename S>
inline Compressible operator*(const Compressible& a, const S& b) {
  return Compressible(a.rho*b, a.rhoU*b, a.e*b);
}

template <typename S>
inline Compressible operator*(const S& b, const Compressible& a) {
  return Compressible(a.rho*b, a.rhoU*b, a.e*b);
}

template <typename S>
inline Compressible operator/(const Compressible& a, const S& b) {
  return Compressible(a.rho/b, a.rhoU/b, a.e/b);
}

inline Compressible operator+=(Compressible& a, const Compressible& b) {
  a.rho+=b.rho; a.rhoU+=b.rhoU; a.e+=b.e;
  return a;
}

inline Compressible operator-=(Compressible& a, const Compressible& b) {
  a.rho-=b.rho; a.rhoU-=b.rhoU; a.e-=b.e;
  return a;
}

template <typename S>
inline Compressible operator*=(Compressible& a, const S& b) {
  a.rho*=b; a.rhoU*=b; a.e*=b;
  return a;
}

template <typename S>
inline Compressible operator/=(Compressible& a, const S& b) {
  a.rho/=b; a.rhoU/=b; a.e/=b;
  return a;
}

#endif
