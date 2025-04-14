#ifndef STEPEXPLICIT_HPP
#define STEPEXPLICIT_HPP

#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "../fvm/computeRez.hpp"
#include "typedefs.hpp"
#include "setting.hpp"
#include "setGhostCells.hpp"
#include "linearSolver.hpp"

template <typename var>
void stepExplicit(CellField<var>& w, CellField<var>& wOld, CellField<var>& rez, const Grid& g,
		  const double& dt, const map<string, bcWithJacobian>& BC,
		  LinearSolver<var>& linSolver, const Setting& setting) {

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
