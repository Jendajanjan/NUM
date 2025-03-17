#include "saveResults.hpp"

void saveResults(const CellField<Compressible>& w, const Grid& g) {

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

  NodeField<double> rho(g);
  NodeField<Vector2d> u(g);
  NodeField<double> e(g);
  NodeField<double> p(g);
  NodeField<double> Ma(g);

  for (int i=0; i<wNode.M(); i++) {
    for (int j=0; j<wNode.N(); j++) {
      rho[i][j] = wNode[i][j].rho;
      u[i][j] = wNode[i][j].rhoU / wNode[i][j].rho;
      e[i][j] = wNode[i][j].e;
      p[i][j] = wNode[i][j].p();
      Ma[i][j] = wNode[i][j].Ma();
    }
  }

  vtk_output v("results/results.vtk");
  v.save(g);
  v.save(rho, "rho");
  v.save(u, "u");
  v.save(e, "e");
  v.save(p, "p");
  v.save(Ma, "Ma");
}
