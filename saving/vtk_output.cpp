#include "vtk_output.hpp"

vtk_output::vtk_output(const string& name) {
  fout.open(name.c_str());
  if (!fout) {
    cout << "File \"" << name << "\" could not be opened!" << endl;
    exit(54);
  }

  fout << "# vtk DataFile Version 3.1" << endl;
  fout << "Results produced by Sim2024" << endl;
  fout << "ASCII" << endl;
  fout << "DATASET STRUCTURED_GRID" << endl;

  fout << scientific << setprecision(8);
}

vtk_output::~vtk_output() {
  fout.close();
}

void vtk_output::save(const Grid& g) {
  fout << "DIMENSIONS " << g.Mnd() << " " << g.Nnd() << " 1" << endl;

  int totalSize = g.Mnd() * g.Nnd();
  
  fout << "POINTS " << totalSize << " DOUBLE" << endl;

  for (int j=0; j<g.Nnd(); j++) {
    for (int i=0; i<g.Mnd(); i++) {
      fout << g.x(i, j) << " " << g.y(i, j) << " 0.0" << endl;
    }
  }

  fout << endl;

  fout << "POINT_DATA " << totalSize << endl;
}

void vtk_output::save(const NodeField<double>& p, const string& name) {
  fout << "SCALARS " << name << " DOUBLE" << endl;
  fout << "LOOKUP_TABLE default" << endl;
  for (int j=0; j<p.N(); j++) {
    for (int i=0; i<p.M(); i++) {
      fout << p[i][j]<< endl;
    }
  }
}

void vtk_output::save(const NodeField<Vector2d>& d, const string& name) {
  fout << "VECTORS " << name << " DOUBLE" << endl;
  for (int j=0; j<d.N(); j++) {
    for (int i=0; i<d.M(); i++) {
      fout << d[i][j].x << " " << d[i][j].y << " 0.0" << endl;
    }
  }
}
