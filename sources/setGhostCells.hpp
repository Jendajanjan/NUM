#ifndef SETGHOSTCELLS_HPP
#define SETGHOSTCELLS_HPP

#include "../compressible.hpp"
#include "../fvm/grid.hpp"
#include "../cellfield.hpp"
#include "setting.hpp"
#inlude "inlet.hpp"
#include "outlet.hpp"
#include "wall.hpp"

void setGhostCells(CellField<Compressible>& w, const Grid& g, const Setting& setting);

#endif
