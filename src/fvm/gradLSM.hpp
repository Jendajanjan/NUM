#ifndef GRADLSM_HPP
#define GRADLSM_HPP

#include <vector>
#include "grid.hpp"
#include "cellfield.hpp"
#include <omp.h>

using namespace std;

template <typename var>
void gradLSM(const CellField<var>& w, CellField<Vector2<var> >& gradW,
	     const Grid& g) {

#pragma omp parallel for
  for (int i=w.Imin(); i<w.Imax(); i++) {
    for (int j=w.Jmin(); j<w.Jmax(); j++) {
      gradW[i][j].x.zero();
      gradW[i][j].y.zero();
    }
  }

#pragma omp parallel for
  for (int i=0; i<w.M(); i++) {
    for (int j=0; j<w.N(); j++) {

      vector<var> wk(8);
      vector<Point2d> centers_k(8);

      const var& w_ij = w[i][j];
      const Point2d& center_ij = g.center(i, j);

      int k=0;
      for (int p=i-1; p<=i+1; p++) {
	for (int r=j-1; r<=j+1; r++) {
	  if (p!=i || r!= j) {
	    wk[k] = w[p][r];
	    centers_k[k] = g.center(p, r);
	    k++;
	  }
	}
      }

      double Ixx = 0.;
      double Ixy = 0.;
      double Iyy = 0.;
      var Jx; Jx.zero();
      var Jy; Jy.zero();

      for (k=0; k<8; k++) {
	double Ix = centers_k[k].x - center_ij.x;
	double Iy = centers_k[k].y - center_ij.y;

	Ixx += Ix * Ix;
	Ixy += Ix * Iy;
	Iyy += Iy * Iy;

	Jx += (wk[k] - w_ij) * Ix;
	Jy += (wk[k] - w_ij) * Iy;
      }

      double D = Ixx * Iyy - Ixy * Ixy;

      gradW[i][j].x = (Jx*Iyy - Jy*Ixy) / D;
      gradW[i][j].y = (Jy*Ixx - Jx*Ixy) / D;
    }
  }
}

#endif
