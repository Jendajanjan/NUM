#include <iostream>

#include "fvm/cellfield.hpp"
#include "fvm/grid.hpp"
#include "fvm/grid_gamm.hpp"
#include "fvm/computeRez.hpp"
#include "sources/initialisation.hpp"
#include "sources/timeStep.hpp"
#include "sources/setGhostCells.hpp"
#include "saving/saveNormRez.hpp"
#include "saving/saveResults.hpp"
#include "compressible.hpp"

using namespace std;

int main() {
  cout << "Welcome in Sim2024!" << endl;

  Setting setting;

  setting.p0 = 1.;
  setting.rho0 = 1.;
  setting.alpha = 0.;
  setting.Ma2is = 0.675;
  setting.CFL = 0.4;
  setting.kappa = 1.4;

  Grid g = Grid_gamm(150, 50, 1);

  CellField<Compressible> w(g);
  CellField<Compressible> wn(g);
  CellField<Compressible> rez(g);

  initialisation(w, setting);

  for (int i=0; i<100000; i++) {
    setGhostCells(w, g, setting);

    double dt = timeStep(w, g, setting);

    computeRez(w, rez, g);

    wn = w + dt * rez;
    w = wn;

    if (i%10 == 0) {
      cout << "iter: " << i << ", dt = " << dt << endl;
      saveNormRez(rez, g, i);
    }
  }

  setGhostCells(w, g, setting);
  saveResults(w, g);

  cout << "Bye Bye!" << endl;

  return 0;
}
