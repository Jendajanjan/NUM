#include "saveNormRez.hpp"

void saveNormRez(const CellField<Compressible>& w, const Grid& g, const int& iter) {
  Compressible L2, Linf;

  L2.zero();
  Linf.zero();

  int M = w.M(), N = w.N();

  for (int i=0; i<M; i++) {
    for (int j=0; j<N; j++) {
      const Compressible& wij = w[i][j];

      Compressible R = Compressible::fabs(wij);

      Linf = Compressible::max(Linf, R);
      L2 += R * R * g.volume(i, j);
    }
  }

  L2 = Compressible::sqrt(L2);

  ofstream fout("results/L2.txt", ios::app);
  fout << iter << " " << L2.rho << " " << L2.rhoU.x << " " << L2.rhoU.y
       << " " << L2.e << endl;

  fout.close();

  fout.open("results/Linf.txt", ios::app);
  fout << iter << " " << Linf.rho << " " << Linf.rhoU.x << " " << Linf.rhoU.y
       << " " << Linf.e << endl;

  fout.close();
}
