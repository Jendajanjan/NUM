#ifndef STEP_HPP
#define STEP_HPP

#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "../fvm/computeRez.hpp"
#include "typedefs.hpp"
#include "setting.hpp"
#include "setGhostCells.hpp"

template <typename var>
void step(CellField<var>& w, CellField<var>& rez, const Grid& g, const double& dt,
	  const map<string, bCondition>& BC, const Setting& setting) {

  CellField<var> wStar(g);

  wStar = w;
  for (int k=0; k<setting.alphaK.size(); k++) {
    setGhostCells(wStar, g, setting, BC);

    computeRez(wStar, rez, g, setting);
      
    wStar = w + setting.alphaK[k] * dt * rez;
  }

  w = wStar;
  
}

#endif
