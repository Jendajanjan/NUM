#include "fluxList.hpp"

FluxList::FluxList() {
  fluxTypes["Upwind"] = Compressible::Upwind;
  fluxTypes["Rusanov"] = Compressible::Rusanov;
  fluxTypes["AUSMUP"] =Compressible::AUSMUP;
}

FluxImplicitList::FluxImplicitList() {
  fluxImplicitTypes["Upwind"] = Compressible::UpwindImplicit;
  fluxImplicitTypes["Rusanov"] = Compressible::RusanovImplicit;
}
