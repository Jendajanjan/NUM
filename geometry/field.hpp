#ifndef FIELD_HPP
#define FIELD_HPP

#include <iostream>
#include <cstdlib>

using namespace std;

template <typename T>
class Field2 {
  int imin;
  int imax;
  int jmin;
  int jmax;

  T* data;
  bool allocated;

public:
  Field2(): imin(0), imax(0), jmin(0), jmax(0) {allocated = false;};

  Field2(const int& _imin, const int& _imax,
	 const int& _jmin, const int& _jmax): imin(_imin), imax(_imax),
					      jmin(_jmin), jmax(_jmax) {
    allocated = true;

    int totalSize = (imax - imin) * (jmax - jmin);
    data = new T[totalSize];
  };

  Field2(const Field2& fld) {
    allocated = false;
    allocate(fld.imin, fld.imax, fld.jmin, fld.jmax);
    
    int totalSize = (imax - imin) * (jmax - jmin);  
    for (int i=0; i<totalSize; i++) {
      data[i] = fld.data[i];
    }
  };

  ~Field2() {
    if (allocated) {
      delete [] data;
    }
    allocated = false;
  };

  void allocate(const int& _imin, const int& _imax,
		const int& _jmin, const int& _jmax) {
    if (allocated) {
      cout << "Pole2 uz je naalokovane, pouzijte funkci reallocate!" << endl;
      exit(1);
    }

    allocated = true;

    imin = _imin;
    imax = _imax;
    jmin = _jmin;
    jmax = _jmax;

    int totalSize = (imax - imin) * (jmax - jmin);
    data = new T[totalSize];
  };

  void reallocate(const int& _imin, const int& _imax,
		  const int& _jmin, const int& _jmax) {
    if (!allocated) {
      cout << "Pole2 neni naalokovane, pouzijte funkci allocate!" << endl;
      exit(1);
    }

    delete [] data;

    allocated = false;

    allocate(_imin, _imax, _jmin, _jmax);
  };

  int Imin() const {return imin;};
  int Imax() const {return imax;};
  int Jmin() const {return jmin;};
  int Jmax() const {return jmax;};

  T* operator[](int i) const {
    int jSize = jmax - jmin;
    return data + (i - imin) * jSize - jmin;
  };

  Field2 operator=(const Field2& fld) {
    if (allocated) {
      reallocate(fld.imin, fld.imax, fld.jmin, fld.jmax);
    }
    else {
      allocate(fld.imin, fld.imax, fld.jmin, fld.jmax);
    }

    int totalSize = (imax - imin) * (jmax - jmin);
    for (int i=0; i<totalSize; i++) {
      data[i] = fld.data[i];
    }

    return *this;
  };
};
  
#endif
