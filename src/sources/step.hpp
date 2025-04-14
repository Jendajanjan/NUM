#ifndef STEP_HPP
#define STEP_HPP

#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "../fvm/computeRez.hpp"
#include "typedefs.hpp"
#include "setting.hpp"
#include "setGhostCells.hpp"
#include "linearSolver.hpp"

template <typename var>
void (*step)(CellField<var>& w, CellField<var>& wOld, CellField<var>& rez, const Grid& g,
	     const double& dt, const map<string, bcWithJacobian>& BC, LinearSolver<var>& linSolver,
	     const Setting& setting);

#endif
