#include "grid_gamm.hpp"

Grid_gamm::Grid_gamm(int M, int N, int ght) {
  Mnodes = M + 1;
  Nnodes = N + 1;
  Mvolumes = M;
  Nvolumes = N;
  ghost = ght;

  nodes.allocate(-ghost, Mnodes+ghost, -ghost, Nnodes+ghost);
  facesI.allocate(0, Mvolumes, 0, Nnodes);
  facesJ.allocate(0, Mnodes, 0, Nvolumes);
  centers.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  volumes.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);

  // generovani souradnic vrcholu
  double dx = 3. / Mvolumes;
  for (int i=0; i<Mnodes; i++) {
    double x = i * dx;
    double y1 = 1.;
    double y0 = 0.;
    if (x >= 1. && x <= 2.) y0 = sqrt(1.69 - pow(x-1.5, 2)) - 1.2;
    double dy = (y1 - y0) / Nvolumes;
    for (int j=0; j<Nnodes; j++) {
      double y = j * dy;
      nodes[i][j].vertex = Point2d(x, y);
    }
  }

  // generovani souradnic pomocnych uzlu
  // uzly ve smeru i
  for (int j=0; j<Nnodes; j++) {
    // leva hranice
    Vector2d s(vertex(1, j), vertex(0, j));
    for (int k=1; k<=ghost; k++) {
      nodes[-k][j].verex = vertex(0, j) + k * s;
    }

    // prava hranice
    s = Vector2d(vertex(Mnodes-2, j), vertex(Mnodes-1, j));
    for (int k=1; k<=ghost; k++) {
      nodes[Mnodes-1+k][j].vertex = vertex(Mnodes-1, j) + k * s;
    }
  }

  // uzly ve smeru j
  for (int i=-ghost; i<Mnodes+ghost; i++) {
    // spodni hranice
    Vector2d s(vertex(i, 1), vertex(i, 0));
    for (int k=1; k<=ghost; k++) {
      nodes[i][-k] = vertex(i, 0) + k * s;
    }

    // horni hranice
    s = Vector2d(vertex(i, Nnodes-2), vertex(i, Nnodes-1));
    for (int k=0; k<=ghost; k++) {
      nodes[i][Nnodes-1+k] = vertex(i, Nnodes-1) + k * s;
    }
  }

  // nastveni sten, objemu a stredu bunek
  update();
}

      
