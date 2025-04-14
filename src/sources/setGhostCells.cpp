#include "setGhostCells.hpp"

void setGhostCells(CellField<Compressible>& w, const Grid& g,
		   const Setting& setting, const map<string, bcWithJacobian>& BC) {
  int M = w.M();
  int N = w.N();
  int gh = w.gh();

  // nastaveni pomocnych bunek na leve a prave hranici
  for (int j=0; j<N; j++) {
    // leva hranice
    Face f = g.faceJ(0, j);

    Compressible wInside = w[0][j];

    auto it = BC.find(f.name);
    Compressible wOutside = it->second.first(wInside, f.s, setting);

    for (int k=1; k<=gh; k++) {
      w[-k][j] = wOutside;
    }

    // prava hranice
    f = g.faceJ(M, j);

    wInside = w[M-1][j];

    it = BC.find(f.name);
    wOutside = it->second.first(wInside, f.s, setting);

    for (int k=1; k<=gh; k++) {
      w[M-1+k][j] = wOutside;
    }
  }

  // dolni a horni hranice
  int imin = -gh;
  int imax = M+gh;

  for (int i=imin; i<imax; i++) {
    // dolni hranice
    Face f = g.faceI(i, 0);

    Compressible wInside = w[i][0];

    auto it = BC.find(f.name);
    Compressible wOutside = it->second.first(wInside, f.s, setting);

    for (int k=1; k<=gh; k++) {
      w[i][-k] = wOutside;
    }


    // horni hranice
    f = g.faceI(i, N);

    wInside = w[i][N-1];
    
    it = BC.find(f.name);
    wOutside = it->second.first(wInside, f.s, setting);

    for (int k=1; k<=gh; k++) {
      w[i][N-1+k] = wOutside;
    }
  }
}
