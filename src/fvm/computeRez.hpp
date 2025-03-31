#ifndef COMPUTEREZ_HPP
#define COMPUTEREZ_HPP

#include "grid.hpp"
#include "cellfield.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include <omp.h>

template <typename var>
void computeRez(const CellField<var>& w, CellField<var>& rez, const Grid& g) {
  int M = w.M();
  int N = w.N();
  int gh = w.gh();

  int imin = w.Imin(), imax=w.Imax();
  int jmin = w.Jmin(), jmax = w.Jmax();

  CellField<var> psi(g);
  CellField<Vector2<var> > gradW(g);

  grad<var>(w, gradW, g);
  limiter<var>(w, gradW, psi, g);

  // vynulovani rezidui
# pragma omp parallel for
  for (int i=imin; i<imax; i++) {
    for (int j=jmin; j<jmax; j++) {
      rez[i][j].zero();
    }
  }

  // cyklus pres steny ve smeru i
# pragma omp parallel for
  for (int i=0; i<M; i++) {
    for (int j=0; j<N+1; j++) {

      const Face& f = g.faceI(i, j);
      Vector2d rL(g.center(i, j), f.center);
      Vector2d rR(g.center(i, j-1), f.center);
      
      var wl = w[i][j] + psi[i][j]
	     * (gradW[i][j].x * rL.x + gradW[i][j].y * rL.y);
      
      var wr = w[i][j-1] + psi[i][j-1]
	     * (gradW[i][j-1].x * rR.x + gradW[i][j-1].y * rR.y);

      var flx = var::flux(wl, wr, f.s);

      rez[i][j] += flx;
      rez[i][j-1] -= flx;
    }
  }

  // cyklus pres steny ve smeru j
# pragma omp parallel for
  for (int j=0; j<N; j++) {
    for (int i=0; i<M+1; i++) {

      const Face& f = g.faceJ(i, j);
      Vector2d rL(g.center(i-1, j), f.center);
      Vector2d rR(g.center(i, j), f.center);
      
      var wl = w[i-1][j] + psi[i-1][j]
	     * (gradW[i-1][j].x * rL.x + gradW[i-1][j].y * rL.y);
      
      var wr = w[i][j] + psi[i][j]
	     * (gradW[i][j].x * rR.x + gradW[i][j].y * rR.y);

      var flx = var::flux(wl, wr, f.s);

      rez[i-1][j] += flx;
      rez[i][j] -= flx;
    }
  }

  // deleni objemem bunky a nasobeni -1x
# pragma omp parallel for
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      rez[i][j] = -1. * rez[i][j] / g.volume(i, j);
    }
  }
}

#endif
