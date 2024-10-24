#include <iostream>
#include "pyqmat.h"
using namespace std;

QMAT::QMAT(const std::string &filename, const std::string &maname) {
  openmeshfile(&input, this, filename, maname);
}

void QMAT::simplify_slab(int n)
{
  simplifySlab(this, &input, n);
}

void QMAT::export_ply(const std::string &fname)
{
  ExportPly(fname, &input);
}

void QMAT::export_ma(const std::string &fname)
{
  Export(fname, &input);
}

std::vector<double> QMAT::hausdorff()
{
  int nv = input.pVertexList.size();
  std::vector<double> hausdoff(nv);
  ComputeHausdorffDistance(*this, input);
  for (int i = 0; i < nv; i++)
  {
    auto it = input.pVertexList[i];
    hausdoff[i] = it->slab_hausdorff_dist;
  }
  return hausdoff;
}
