#ifndef SAVERESULTS_HPP
#define SAVERESULTS_HPP

#include <fstream>
#include "../fvm/grid.hpp"
#include "../fvm/cellfield.hpp"
#include "../fvm/nodefield.hpp"
#include "../compressible.hpp"

using namespace std;

void saveResults(const CellField<Compressible>& w, const Grid& g);

#endif
