#ifndef SETNODEFIELD_HPP
#define SETNODEFIELD_HPP

#include "grid.hpp"
#include "cellfield.hpp"
#include "nodefield.hpp"
#include <omp.h>

template <typename var>
void setNodeField(const CellField<var>& w, NodeField<var>& wNode, const Grid& g) {
#pragma omp parallel for
  for (int i=0; i<wNode.M(); i++) {
    for (int j=0; j<wNode.N(); j++) {
      var& wNd = wNode[i][j];
      wNd.zero();

      vector<var> wK(4);
      wK[0] = w[i][j];
      wK[1] = w[i-1][j];
      wK[2] = w[i-1][j-1];
      wK[3] = w[i][j-1];

      vector<double>& alpha = g.alpha(i, j);

      for (int k=0; k<4; k++) {
	wNd += wK[k] * alpha[k];
      }
    }
  }
}

#endif
