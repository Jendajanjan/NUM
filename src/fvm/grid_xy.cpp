#include "grid_xy.hpp"

void Grid_xy::loadXYFiles(const string& name1, const string& name2,
		   vector<double>& xCoord, vector<double>& yCoord) {

  ifstream fin(name1.c_str());
  if (!fin) {
    cout << "File \"" << name1 << "\" with x coordinates can not be opened!" << endl;
    exit(20);
  }

  int nb;
  fin >> nb;

  xCoord.resize(nb);
  for (int i=0; i<nb; i++) {
    fin >> xCoord[i];
  }

  fin.close();

  fin.open(name2.c_str());
  
  if (!fin) {
    cout << "File \"" << name2 << "\" with y coordinates can not be opened!" << endl;
    exit(20);
  }

  fin >> nb;

  yCoord.resize(nb);
  for (int i=0; i<nb; i++) {
    fin >> yCoord[i];
  }

  fin.close();
}

Grid_xy::Grid_xy(const string name1, const string name2, const int gh, const string& type) {
  vector<double> xCoord, yCoord;

  loadXYFiles(name1, name2, xCoord, yCoord);
    
  //m, n pocet bunek ve smeru i, j
  Mnodes = xCoord.size();
  Nnodes = yCoord.size();
  Mvolumes = Mnodes - 1;
  Nvolumes = Nnodes -1;
  ghost = gh;

  nodes.allocate(-ghost, Mnodes+ghost, -ghost, Nnodes+ghost);
  centers.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  volumes.allocate(-ghost, Mvolumes+ghost, -ghost, Nvolumes+ghost);
  facesI.allocate(-ghost, Mvolumes+ghost, -ghost, Nnodes+ghost);
  facesJ.allocate(-ghost, Mnodes+ghost, -ghost, Nvolumes+ghost);

  //generovani souradnic
  
  for (int i=0; i<Mnodes; i++) {
    for (int j=0; j<Nnodes; j++) {
      nodes[i][j].vertex = Point2d(xCoord[i], yCoord[j]);
    }
  }

  // vytvoreni pomocnych bunek
  for (int j=0; j<Nnodes; j++) {
    for (int k=1; k<=ghost; k++) {
      //leva hranice
      Vector2d s(vertex(1, j), vertex(0, j));
       nodes[-k][j].vertex = vertex(0, j) + k * s;

      //prava hranice
      s = Vector2d(vertex(Mnodes-2, j), vertex(Mnodes-1, j));
      nodes[Mnodes-1+k][j].vertex = vertex(Mnodes-1, j) + k * s;
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int k=1; k<=ghost; k++) {
      //dolni hranice
      Vector2d s(vertex(i, 1), vertex(i, 0));
      nodes[i][-k].vertex = vertex(i, 0) + k * s;

      //horni hranice
      s = Vector2d(vertex(i, Nnodes-2), vertex(i, Nnodes-1));
      nodes[i][Nnodes-1+k].vertex = vertex(i, Nnodes-1) + k * s;
    }
  }

  // nastveni sten, objemu a stredu bunek
  update();

  map<string, coefficients>::iterator it;

  it = coeffList.find(type);

  if (it != coeffList.end()) {
    computeAlpha = it->second;
  }
  else {
    cout << "No a such possibility for a node weight computation!" << endl;
    cout << "Possibilities are: ";
    for (it=coeffList.begin(); it!=coeffList.end(); it++) {
      cout << it->first << ", ";
    }
    cout << endl;
    exit(21);
  }

  computeAlpha(*this);

  for (int i=-ghost; i<Mvolumes+ghost; i++) {
    for (int j=-ghost; j<Nnodes+ghost; j++) {
      facesI[i][j].name = "internal";

      if (j == 0) {
	if (i < Mnodes/2) {
	  facesI[i][j].name = "Lower";
	}
	else {
	  facesI[i][j].name = "Wall";
	}
      }
      
      if (j == (Nnodes-1)) {
	facesI[i][j].name = "Upper";
      }
    }
  }

  for (int i=-ghost; i<Mnodes+ghost; i++) {
    for (int j=-ghost; j<Nvolumes+ghost; j++) {
      facesJ[i][j].name = "internal";

      if (i == 0) facesJ[i][j].name = "Inlet";
      if (i == (Mnodes-1)) facesJ[i][j].name = "Outlet";
    }
  }
}
