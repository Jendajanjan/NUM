#include "initialisation.hpp"

void initialisation(CellField<Compressible>& w, const Setting& setting) {
  system("mkdir -p results");
  system("rm -f results/*");

  Compressible::kappa = setting.kappa;
  Compressible::flux = Compressible::Upwind;

  double rhoInit = 0.9;
  Vector2d uInit(0.5, 0.);
  double pInit = 0.8;

  double eInit = pInit / (Compressible::kappa - 1.)
               + 0.5 * rhoInit * (uInit.x*uInit.x + uInit.y*uInit.y);

  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N(); j++) {
      w[i][j] = Compressible(rhoInit, rhoInit*uInit, eInit);
    }
  }
}

