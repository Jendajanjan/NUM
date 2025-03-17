#ifndef VTK_OUTPUT_HPP
#define VTK_OUTPUT_HPP

#include <fstream>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include "../fvm/grid.hpp"
#include "../fvm/nodefield.hpp"

using namespace std;

class vtk_output {
  ofstream fout;
  
public:
  vtk_output(const string& name);
  ~vtk_output();

  void save(const Grid& g);
  void save(const NodeField<double>& p, const string& name);
  void save(const NodeField<Vector2d>& d, const string& name);
};

#endif
