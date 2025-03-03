#include "setGhostCells.hpp"

void setGhostCells(CellField<Compressible>& w, const Grid& g, const Setting& setting) {
  int M = w.M();
  int N = w.N();
  int gh = w.gh();

  // nastaveni pomocnych bunek na leve a prave hranici
  for (int j=0; j<N; j++) {
    // leva hranice
    Vector2d s = g.faceJ(0, j).s;

    Compressible wInside = w[0][j];
    Compressible wOutside = inlet(wInside, s, setting);

    for (int k=1; k<=gh; k++) {
      w[-k][j] = wOutside;
    }

    // prava hranice
    s = g.faceJ(M, j).s;

    wInside = w[M-1][j];
    wOutside = outlet(wInside, s, setting);

    for (int k=1; k<=gh; k++) {
      w[M-1+k][j] = wOutside;
    }
  }

  // dolni a horni hranice
  int imin = -gh;
  int imax = M+gh;

  for (int i=imin; i<imax; i++) {
    // dolni hranice
    Vector2d s = g.faceI(i, 0).s;

    Compressible wInside = w[i][0];
    Compressible wOutside = wall(wInside, s, setting);

    for (int k=1; k<=gh; k++) {
      w[i][-k] = wOutside;
    }


    // horni hranice
    s = g.faceI(i, N).s;

    wInside = w[i][N-1];
    wOutside = wall(wInside, s, setting);

    for (int k=1; k<=gh; k++) {
      w[i][N-1+k] = wOutside;
    }
  }
}
