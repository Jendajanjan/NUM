#ifndef SAVENORMREZ_HPP
#define SAVENORMREZ_HPP

#include <fstream>
#include <iomanip>
#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "../compressible.hpp"

using namespace std;

void saveNormRez(const CellField<Compressible>& w, const Grid& g, const int& iter);

#endif
