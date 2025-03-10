#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include "../fvm/cellfield.hpp"
#include "../compressible.hpp"
#include "setting.hpp"

void initialisation(CellField<Compressible>& w, const Setting& setting);

#endif
