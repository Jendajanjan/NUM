#ifndef INITIALISATION_HPP
#define INITIALISATION_HPP

#include <cstdlib>
#include "../fvm/cellfield.hpp"
#include "fvm/limiter.hpp"
#include "fvm/grad.hpp"
#include "fvm/zeroLimiter.hpp"
#include "fvm/zeroGrad.hpp"
#include "fvm/gradLSM.hpp"
#include "fvm/barthJespersen.hpp"
#include "fvm/venkatakrishnan.hpp"
#include "../compressible.hpp"
#include "setting.hpp"

using namespace std;

void initialisation(CellField<Compressible>& w, const Setting& setting);

#endif
