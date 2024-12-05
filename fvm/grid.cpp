#include "grid.hpp"

void Face::updateI(Grid& g, int i, int j) {
  const Point2d& a = g.vertex(i, j);
  const Point2d& b = g.vertex(i+1, j);

  center = (a + b) / 2.;
  s = Vector2d(b.y-a.y, a.x-b.x);
}

void Face::updateJ(Grid& g, int i, int j) {
  const Point2d& a = g.vertex(i, j);
  const Point2d& b = g.vertex(i, j+1);

  center = (a + b) / 2.;
  s = Vector2d(b.y-a.y, a.x-b.x);
}

double Grid::x(const int& i, const int& j) const {
  return nodes[i][j].vertex.x;
}

double Grid::y(const int& i, const int& j) const {
  return nodes[i][j].vertex.y;
}

Point2d& Grid::vertex(const int& i, const int& j) const {
  return nodes[i][j].vertex;
}

Node& Grid::node(const int& i, const int& j) const {
  return nodes[i][j];
}

Face& Grid::faceI(const int& i, const int& j) const {
  return facesI[i][j];
}

Face& Grid::faceJ(const int& i, const int& j) const {
  return facesJ[i][j];
}

Point2d& Grid::center(const int& i, const int& j) const {
  return centers[i][j];
}

double Grid::volume(const int& i, const int& j) const {
  return volumes[i][j];
}

void Grid::updateVolumes() {
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      const Point2d& A = vertex(i, j);
      const Point2d& B = vertex(i+1, j);
      const Point2d& C = vertex(i+1, j+1);
      const Point2d& D = vertex(i, j+1);

      Vector2d f(A, C);
      Vector2d e(B, D);

      volumes[i][j] = fabs(cross(f, e)) / 2.;
    }
  }
}

void Grid::updateCenters() {
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      centers[i][j] = (vertex(i, j) + vertex(i+1, j) + vertex(i+1, j+1) +
		       vertex(i, j+1)) / 4.;
    }
  }
}

void Grid::update() {
  for (int i=0; i<Mvolumes; i++) {
    for (int j=0; i<Nnodes; j++) {
      facesI[i][j].updateI(*this, i, j);
    }
  }

  for (int i=0; i<Mnodes; i++) {
    for (int j=0; j<Nvolumes; j++) {
      facesJ[i][j].updateJ(*this, i, j);
    }
  }

  updateCenters();
  updateVolumes();
}
