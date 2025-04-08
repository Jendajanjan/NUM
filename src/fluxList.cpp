#include "fluxList.hpp"

FluxList::FluxList() {
  fluxTypes["Upwind"] = Compressible::Upwind;
  fluxTypes["Rusanov"] = Compressible::Rusanov;
}
