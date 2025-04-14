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
#include "sources/setGrid.hpp"
#include "sources/step.hpp"
#include "sources/linearSolver.hpp"
#include "compressible.hpp"

using namespace std;

int main(int argc,char **args) {
  
  PetscInitialize( &argc , &args , (char *)0 , 0 );
  
  cout << "Welcome in Sim2024!" << endl;

  Setting setting("starter.txt");

  Grid g;
  setGrid(g, setting);

  map<string, bcWithJacobian> BC;
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

  LinearSolver<Compressible> linSolver(setting.solver, g);

  CellField<Compressible> w(g);
  CellField<Compressible> wOld(g);
  CellField<Compressible> rez(g);

  initialisation(w, setting);
  wOld = w;
  
  for (int i=0; i<setting.stop; i++) {

    double dt = timeStep(w, g, setting);

    step<Compressible>(w, wOld, rez, g, dt, BC, linSolver, setting);

    if (i%10 == 0) {
      cout << "iter: " << i << ", dt = " << dt << endl;
      saveNormRez(rez, g, i);
    }
  }

  setGhostCells(w, g, setting, BC);
  saveResults(w, g);

  cout << "Bye Bye!" << endl;

  linSolver.free();

  PetscFinalize();

  return 0;
}
