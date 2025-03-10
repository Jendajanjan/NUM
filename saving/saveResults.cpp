#include "saveResults.hpp"

void saveResults(const CellField<Compressible>& w, const Grid& g) {
  ofstream fout("results/results.vtk");

  fout << "# vtk DataFile Version 3.1" << endl;
  fout << "Results produced by Sim2024" << endl;
  fout << "ASCII" << endl;
  fout << "DATASET STRUCTURED_GRID" << endl;

  fout << "DIMENSIONS " << g.Mnd() << " " << g.Nnd() << " 1" << endl;

  int totalSize = g.Mnd() * g.Nnd();
  
  fout << "POINTS " << totalSize << " DOUBLE" << endl;

  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << g.x(i, j) << " " << g.y(i, j) << " 0.0" << endl;
    }
  }

  fout << endl;

  NodeField<Compressible> wNode(g);

  for (int i=0; i<wNode.M(); i++) {
    for (int j=0; j<wNode.N(); j++) {
      Compressible& wNij = wNode[i][j];
      
      double volume = g.volume(i, j) + g.volume(i-1, j)
	            + g.volume(i-1, j-1) + g.volume(i, j-1);

      wNij = w[i][j] * g.volume(i, j)
	   + w[i-1][j] * g.volume(i-1, j)
	   + w[i-1][j-1] * g.volume(i-1, j-1)
	   + w[i][j-1] * g.volume(i, j-1);

      wNij /= volume;
    }
  }

  fout << "POINT_DATA " << totalSize << endl;
  fout << "SCALARS rho DOUBLE" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << wNode[i][j].rho << endl;
    }
  }

  fout << "VECTORS u DOUBLE" << endl;
  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << wNode[i][j].rhoU.x/w[i][j].rho
	   << " " << wNode[i][j].rhoU.y/w[i][j].rho
	   << " 0.0" << endl;
    }
  }

  fout << "SCALARS e DOUBLE" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << wNode[i][j].e << endl;
    }
  }

  fout << "SCALARS p DOUBLE" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << wNode[i][j].p() << endl;
    }
  }

  fout << "SCALARS Ma DOUBLE" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << wNode[i][j].Ma() << endl;
    }
  }

  fout.close();
}
