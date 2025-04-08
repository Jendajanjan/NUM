#ifndef SLIPWALL_HPP
#define SLIPWALL_HPP

#include <cmath>
#include "../geometry/vector.hpp"
#include "../compressible.hpp"
#include "setting.hpp"

using namespace std;

Compressible slipWall(const Compressible& wInside, const Vector2d& s, const Setting& setting);

#endif
