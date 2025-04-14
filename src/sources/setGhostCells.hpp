#ifndef SETGHOSTCELLS_HPP
#define SETGHOSTCELLS_HPP

#include "../compressible.hpp"
#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "setting.hpp"
#include "typedefs.hpp"

using namespace std;

void setGhostCells(CellField<Compressible>& w, const Grid& g,
		   const Setting& setting, const map<string, bcWithJacobian>& BC);

#endif
