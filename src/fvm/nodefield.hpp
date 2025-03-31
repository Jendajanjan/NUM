#ifndef NODEFIELD_HPP
#define NODEFIELD_HPP

#include <iostream>
#include <cstdlib>
#include "grid.hpp"

using namespace std;

template<typename var>
class NodeField {
  int Mnd;
  int Nnd;
  int ghost;

  Field2<var> data;
public:
  NodeField(const Grid& g) {
    Mnd = g.Mnd();
    Nnd = g.Nnd();
    ghost = g.gh();

    data.allocate(-ghost, Mnd+ghost, -ghost, Nnd+ghost);
  };

  NodeField(int M, int N, int gh): Mnd(M), Nnd(N), ghost(gh) {
    data.allocate(-ghost, Mnd+ghost, -ghost, Nnd+ghost);
  }

  ~NodeField() {};

  inline var* operator[](int i) const {
    return data[i];
  }

  int M() const {return Mnd;};
  int N() const {return Nnd;};
  int gh() const {return ghost;};
  int Imin() const {return data.Imin();};
  int Imax() const {return data.Imax();};
  int Jmin() const {return data.Jmin();};
  int Jmax() const {return data.Jmax();};
};

template <typename var>
inline NodeField<var> operator+(const NodeField<var>& a, const NodeField<var>& b) {
  if (a.M() != b.M() || a.N() != b.N() || a.gh() != b.gh()) {
    cout << "Nelze scitat pole NodeField ruznych velikosti!" << endl;
    exit(10);
  }

  NodeField<var> c(a.M(), a.N(), a.gh());
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] + b[i][j];
    }
  }

  return c;
}

template <typename var>
inline NodeField<var> operator-(const NodeField<var>& a, const NodeField<var>& b) {
  if (a.M() != b.M() || a.N() != b.N() || a.gh() != b.gh()) {
    cout << "Nelze odecitat pole NodeField ruznych velikosti!" << endl;
    exit(10);
  }

  NodeField<var> c(a.M(), a.N(), a.gh());
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] - b[i][j];
    }
  }

  return c;
}

template <typename var, typename S>
inline NodeField<var> operator*(const NodeField<var>& a, const S& b) {
  NodeField<var> c(a.M(), a.N(), a.gh());
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}


template <typename var, typename S>
inline NodeField<var> operator*(const S& b, const NodeField<var>& a) {
  NodeField<var> c(a.M(), a.N(), a.gh());
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] * b;
    }
  }

  return c;
}

template <typename var, typename S>
inline NodeField<var> operator/(const NodeField<var>& a, const S& b) {
  NodeField<var> c(a.M(), a.N(), a.gh());
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      c[i][j] = a[i][j] / b;
    }
  }

  return c;
}

template <typename var>
inline NodeField<var> operator+=(NodeField<var>& a, const NodeField<var>& b) {
  if (a.M() != b.M() || a.N() != b.N() || a.gh() != b.gh()) {
    cout << "Nelze scitat (+=) pole NodeField ruznych velikosti!" << endl;
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
inline NodeField<var> operator-=(NodeField<var>& a, const NodeField<var>& b) {
  if (a.M() != b.M() || a.N() != b.N() || a.gh() != b.gh()) {
    cout << "Nelze odecitat (-=) pole NodeField ruznych velikosti!" << endl;
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
inline NodeField<var> operator*=(NodeField<var>& a, const S& b) {
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] *= b;
    }
  }

  return a;
}

template <typename var, typename S>
inline NodeField<var> operator/=(NodeField<var>& a, const S& b) {
  
  for (int i=a.Imin(); i<a.Imax(); i++) {
    for (int j=a.Jmin(); j<a.Jmax(); j++) {
      a[i][j] /= b;
    }
  }

  return a;
}

#endif
