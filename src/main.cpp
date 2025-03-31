#include <iostream>
#include <cstdlib>

#include "fvm/cellfield.hpp"
#include "fvm/grid.hpp"
#include "fvm/grid_gamm.hpp"
#include "fvm/computeRez.hpp"
#include "sources/initialisation.hpp"
#include "sources/timeStep.hpp"
#include "sources/setGhostCells.hpp"
#include "sources/BC.hpp"
#include "saving/saveNormRez.hpp"
#include "saving/saveResults.hpp"
#include "compressible.hpp"

using namespace std;

int main() {
  cout << "Welcome in Sim2024!" << endl;

  Setting setting("starter.txt");

  Grid g = Grid_gamm(setting.mCells, setting.nCells, setting.ghostCells);

  map<string, bCondition> BC;
  for (auto iter=setting.usedBC.begin(); iter!=setting.usedBC.end(); iter++) {
    auto bCond = bcList.find(iter->second);
    if (bCond!=bcList.end()) {
      BC[iter->first] = bcList[iter->second];
    }
    else {
      cout << "Not found \"" << iter->second << "\" boundary condition!" << endl;
      cout << "Possibilities are:" << endl;
      for (bCond = bcList.begin(); bCond!=bcList.end(); bCond++) {
	cout << bCond->first << endl;
      }
      exit(52);
    }
  }

  CellField<Compressible> w(g);
  CellField<Compressible> wn(g);
  CellField<Compressible> wStar(g);
  CellField<Compressible> rez(g);

  initialisation(w, setting);

  for (int i=0; i<setting.stop; i++) {

    double dt = timeStep(w, g, setting);

    wStar = w;
    for (int k=0; k<setting.alphaK.size(); k++) {
      setGhostCells(wStar, g, setting, BC);

      computeRez(wStar, rez, g);

      wStar = w + dt * rez;
    }
    wn = wStar;
    w = wn;

    if (i%10 == 0) {
      cout << "iter: " << i << ", dt = " << dt << endl;
      saveNormRez(rez, g, i);
    }
  }

  setGhostCells(w, g, setting, BC);
  saveResults(w, g);

  cout << "Bye Bye!" << endl;

  return 0;
}
