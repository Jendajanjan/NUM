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

vector<double>& Grid::alpha(const int& i, const int& j) const {
  return nodes[i][j].alpha;
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
  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nnodes+ghost; j++) {
      facesI[i][j].updateI(*this, i, j);
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      facesJ[i][j].updateJ(*this, i, j);
    }
  }

  updateCenters();
  updateVolumes();
}

void Grid::computeAlphaWeight(Grid& g) {
  for (int i=0; i<g.Mnd(); i++) {
    for (int j=0; j<g.Nnd(); j++) {
      double volume = g.volume(i, j) + g.volume(i-1, j)
	            + g.volume(i-1, j-1) + g.volume(i, j-1);

      Node& nd = g.node(i, j);
      nd.alpha.resize(4);

      nd.alpha[0] = g.volume(i, j) / volume;
      nd.alpha[1] = g.volume(i-1, j) / volume;
      nd.alpha[2] = g.volume(i-1, j-1) / volume;
      nd.alpha[3] = g.volume(i, j-1) / volume;
    }
  }
}

void Grid::computeAlphaLSM(Grid& g) {
  for (int i=0; i<g.Mnd(); i++) {
    for (int j=0; j<g.Nnd(); j++) {
      Node& nd = g.node(i, j);
      nd.alpha.resize(4);

      vector<Point2d> centers(4);
      centers[0] = g.center(i, j);
      centers[1] = g.center(i-1, j);
      centers[2] = g.center(i-1, j-1);
      centers[3] = g.center(i, j-1);

      const Point2d& V = nd.vertex;

      double Rx=0., Ry=0., Ixx=0., Ixy=0., Iyy=0.;

      for (int k=0; k<4; k++) {
	double Jx = centers[k].x - V.x;
	double Jy = centers[k].y - V.y;

	Rx += Jx;
	Ry += Jy;

	Ixx += Jx * Jx;
	Ixy += Jx * Jy;
	Iyy += Jy * Jy;
      }

      double D = Ixx * Iyy - Ixy * Ixy;
      double lx = (Ry * Ixy - Rx * Iyy) / D;
      double ly = (Rx * Ixy - Ry * Ixx) / D;

      for (int k=0; k<4; k++) {
	nd.alpha[k] = (1. + lx * (centers[k].x - V.x) + ly * (centers[k].y - V.y))
	            / (4. + lx * Rx + ly * Ry);
      }
    }
  }
}

void (*Grid::computeAlpha)(Grid& g);

map<string, coefficients> Grid::coeffList = {pair<string, coefficients>("Weight", computeAlphaWeight),
					      pair<string, coefficients>("LSM", computeAlphaLSM)};

