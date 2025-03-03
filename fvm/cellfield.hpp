#ifndef CELLFIELD_HPP
#define CELLFIELD_HPP

#include <iostream>
#include <cstdlib>
#include "grid.hpp"

using namespace std;

template<typename var>
class CellField {
  int Mvol;
  int Nvol;
  int ghost;

  Field2<var> data;
public:
  CellField(const Grid& g) {
    Mvol = g.Mvol();
    Nvol = g.Nvol();
    ghost = g.gh();

    data.allocate(-ghost, Mvol+ghost, -ghost, Nvol+ghost);
  };

  CellFiled(int M, int N, int gh): Mvol(M), Nvol(N), ghost(gh) {
    data.allocate(-ghost, Mvol+ghost, -ghost, Nvol+ghost);
  }

  ~CellField() {};

  inline var* operator[](int i) {
    return data[i];
  }

  int M() const {return Mvol;};
  int N() const {return Nvol;};
  int gh() const {return ghost;};
  int Imin() const {return data.Imin();};
  int Imax() const {return data.Imax();};
  int Jmin() const {return data.Jmin();};
  int Jmax() const {return data.Jmax();};
};

template <typename var>
inline CellField<var> opretator+(const CellField<var>& a, const CellField<var>& b) {
  if (a.M() != b.M() || a.N() != b.n() || a.gh() != b.gh()) {
    cout << "Nelze scitat pole CellField ruznych velikosti!" << endl;
    exit(10);
  }

  CellField<var> c(a.M(), a.N(), a.gh());
  for (int i=Imin(); i<Imax(); i++) {
    for (int j=Jmin(); j<Jmax(); j++) {
      c[i][j] = a[i][j] + b[i][j];
    }
  }

  return c;
}

template <typename var>
inline CellField<var> opretator-(const CellField<var>& a, const CellField<var>& b) {
  if (a.M() != b.M() || a.N() != b.n() || a.gh() != b.gh()) {
    cout << "Nelze odecitat pole CellField ruznych velikosti!" << endl;
    exit(10);
  }

  CellField<var> c(a.M(), a.N(), a.gh());
  for (int i=Imin(); i<Imax(); i++) {
    for (int j=Jmin(); j<Jmax(); j++) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }

  return c;
}

template <typename var, typename S>
inline CellField operator*(const CellField<var>& a, const S& b) {
  CellField c(a.M(), a.N(), a.gh());
  
  for (int i=Imin(); i<Imax(); i++) {
    for (int j=Jmin(); j<Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}


template <typename var, typename S>
inline CellField operator*(const S& b, const CellField<var>& a) {
  CellField c(a.M(), a.N(), a.gh());
  
  for (int i=Imin(); i<Imax(); i++) {
    for (int j=Jmin(); j<Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}

template <typename var, typename S>
inline CellField operator/(const CellField<var>& a, const S& b) {
  CellField c(a.M(), a.N(), a.gh());
  
  for (int i=Imin(); i<Imax(); i++) {
    for (int j=Jmin(); j<Jmax(); j++) {
      c[i][j] = a[i][j] / b;
    }
  }

  return c;
}

template <typename var>
inline CellField<var> opretator+=(CellField<var>& a, const CellField<var>& b) {
  if (a.M() != b.M() || a.N() != b.n() || a.gh() != b.gh()) {
    cout << "Nelze scitat (+=) pole CellField ruznych velikosti!" << endl;
    exit(10);
  }

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] += b[i][j];
    }
  }

  return a;
}

template <typename var>
inline CellField<var> opretator-=(CellField<var>& a, const CellField<var>& b) {
  if (a.M() != b.M() || a.N() != b.n() || a.gh() != b.gh()) {
    cout << "Nelze odecitat (-=) pole CellField ruznych velikosti!" << endl;
    exit(10);
  }

  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] -= b[i][j];
    }
  }

  return a;
}

template <typename var, typename S>
inline CellField operator*=(CellField<var>& a, const S& b) {
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] *= b;
    }
  }

  return a;
}

template <typename var, typename S>
inline CellField operator/=(CellField<var>& a, const S& b) {
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] /= b;
    }
  }

  return a;
}

#endif
