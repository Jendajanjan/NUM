#ifndef COMPUTEREZ_HPP
#define COMPUTEREZ_HPP

#include "grid.hpp"
#include "cellfield.hpp"
#include <omp.h>

template <typename var>
void computeRez(const CellField<var>& w, CellField<var>& rez, const Grid& g) {
  int M = w.M();
  int N = w.N();
  int gh = w.gh();

  int imin = w.Imin(), imax=w.Imax();
  int jmin = w.Jmin(), jmax = w.Jmax();

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
      const var& wl = w[i][j];
      const var& wr = w[i][j-1];
      const Vector2d& s = g.faceI(i, j).s;

      var flx = var::flux(wl, wr, s);

      rez[i][j] += flx;
      rez[i][j-1] -= flx;
    }
  }

  // cyklus pres steny ve smeru j
# pragma omp parallel for
  for (int j=0; j<N; j++) {
    for (int i=0; i<M+1; i++) {
      const var& wl = w[i-1][j];
      const var& wr = w[i][j];
      const Vector2d& s = g.faceJ(i, j).s;

      var flx = var::flux(wl, wr, s);

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
