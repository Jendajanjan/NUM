#ifndef COMPUTEREZ_HPP
#define COMPUTEREZ_HPP

#include "grid.hpp"
#include "cellfield.hpp"
#include "computeResidualConv.hpp"
#include "computeResidualDiss.hpp"
#include "grad.hpp"
#include "limiter.hpp"
#include "sources/setting.hpp"
#include <omp.h>

template <typename var>
void computeRez(const CellField<var>& w, CellField<var>& rez, const Grid& g, const Setting& setting) {
  int M = w.M();
  int N = w.N();
  
  int imin = w.Imin(), imax=w.Imax();
  int jmin = w.Jmin(), jmax = w.Jmax();

  // vynulovani rezidui
# pragma omp parallel for
  for (int i=imin; i<imax; i++) {
    for (int j=jmin; j<jmax; j++) {
      rez[i][j].zero();
    }
  }

  switch (setting.convection) {
  case 0: break;     // konvekce je vypnuta, nedelame nic
  case 1: computeResidualConv(w, rez, g);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(71);
  }

  switch (setting.diffusion) {
  case 0: break;     // difuze je vypnuta, nedelame nic
  case 1: computeResidualDiss(w, rez, g);
    break;
  default: cout << "No such possibility for convection!" << endl;
    exit(71);
  }

  // deleni objemem bunky
# pragma omp parallel for
  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      rez[i][j] = rez[i][j] / g.volume(i, j);
    }
  }
}

#endif
