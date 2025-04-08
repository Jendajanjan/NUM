#ifndef FLUXLIST_HPP_
#define FLUXLIST_HPP_

#include <iostream>
#include <map>
#include <set>
#include <string>
#include <cstdlib>
#include <utility>
#include "compressible.hpp"

using namespace std;

typedef Compressible (*fluxType)(const Compressible& wl, const Compressible& wr, const Vector2d& s);

class FluxList {
public:
  map<string, fluxType> fluxTypes;

  FluxList();
  ~FluxList() {};

  fluxType operator[](const string& name) {
    map<string, fluxType>::iterator it;
    it = fluxTypes.find(name);
    if (it != fluxTypes.end()) {
      return fluxTypes[name];
    }
    else {
      cout << "There is no \"" << name << "\" flux type for compressible variables!" << endl;
      cout << "Possible choice is: " << endl;
      for (it=fluxTypes.begin(); it!=fluxTypes.end(); it++) {
	cout << it->first << endl;
      }
      exit(0);
    }
  }
};

#endif
