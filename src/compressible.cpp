#include "compressible.hpp"

double Compressible::kappa;
double Compressible::R;
double Compressible::cp;
double Compressible::cv;
double Compressible::Pr;

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

double Compressible::T() const {
  return p() / (R * rho);
}

double Compressible::mu() const {
  return (1.45 * pow(T(), 3./2.)) / (T() + 110) * 1e-6;
}

double Compressible::k() const {
  return cp * mu() / Pr;
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

Compressible Compressible::Rusanov(const Compressible& wl, const Compressible& wr, const Vector2d& s) {
  // Vector2d n = s / s.length();

  // double unL = dot(wl.rhoU/wl.rho, n);
  // double unR = dot(wr.rhoU/wr.rho, n);

  // double lambdaL = std::fabs(unL) + wl.a();
  // double lambdaR = std::fabs(unR) + wr.a();

  // double Slambda = std::max(lambdaL, lambdaR);

  // Compressible Fl(wl.rho*unL, wl.rhoU*unL + wl.p()*n, (wl.e + wl.p()) * unL);
  // Compressible Fr(wr.rho*unR, wr.rhoU*unR + wr.p()*n, (wr.e + wr.p()) * unR);

  // Compressible flx = 0.5 * (Fl + Fr) - 0.5 * Slambda * (wr - wl);

  // return flx * s.length();

  
  // Standardne: 1)rotace do n,t, 2)vypcoteni toku, 3) zpetna rotace do x,y
  Vector2d n = s / s.length();
  Vector2d t(-n.y, n.x);

  Compressible WR = wr;
  Compressible WL = wl;

  double rhouR = dot(wr.rhoU, n);
  double rhovR = dot(wr.rhoU, t);
  WR.rhoU = Vector2d(rhouR, rhovR);

  double rhouL = dot(wl.rhoU, n);
  double rhovL = dot(wl.rhoU, t);
  WL.rhoU = Vector2d(rhouL, rhovL);

  double lambdaL = std::fabs(WL.rhoU.x/WL.rho) + WL.a();
  double lambdaR = std::fabs(WR.rhoU.x/WR.rho) + WR.a();

  double Slambda = std::max(lambdaL, lambdaR);

  Compressible Fl(WL.rhoU.x, Vector2d(pow(WL.rhoU.x,2)/WL.rho + WL.p(), WL.rhoU.x*WL.rhoU.y/WL.rho),
		  (WL.e + WL.p()) * WL.rhoU.x/WL.rho);
  Compressible Fr(WR.rhoU.x, Vector2d(pow(WR.rhoU.x,2)/WR.rho + WR.p(), WR.rhoU.x*WR.rhoU.y/WR.rho),
		  (WR.e + WR.p()) * WR.rhoU.x/WR.rho);

  Compressible flx = 0.5 * (Fl + Fr) - 0.5 * Slambda * (WR - WL);

  Vector2d nInv(n.x, -n.y);
  Vector2d tInv(-t.x, t.y);

  double RHOU = dot(flx.rhoU, nInv);
  double RHOV = dot(flx.rhoU, tInv);

  flx.rhoU = Vector2d(RHOU, RHOV);

  return flx * s.length();
}

Compressible Compressible::fluxDissipative(const Vector2<PrimitiveVars>& gradP, const Compressible& wFace, const Vector2d& s) {
  double mu = wFace.mu();
  double k = wFace.k();

  Vector2d u = wFace.rhoU / wFace.rho;

  double divU = gradP.x.u.x + gradP.y.u.y;   // prvni .x - derivace podle x, druhe .x - x-ova slozka rychlosti

  double Txx = 2. * mu * (gradP.x.u.x - 1./3. * divU);
  double Tyy = 2. * mu * (gradP.y.u.y - 1./3. * divU);
  double Txy = mu * (gradP.y.u.x + gradP.x.u.y);

  Vector2d firstRow(Txx, Txy);
  Vector2d secondRow(Txy, Tyy);

  double Theta_x = dot(firstRow, u) + k * gradP.x.T;
  double Theta_y = dot(secondRow, u) + k * gradP.y.T;

  double firstComponent = dot(firstRow, s);
  double secondComponent = dot(secondRow, s);
  double last = dot(Vector2d(Theta_x, Theta_y), s);

  return Compressible(0., Vector2d(firstComponent, secondComponent), last);
}
