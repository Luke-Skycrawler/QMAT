#include <iostream>
#include "pyqmat.h"
using namespace std;

QMAT::QMAT(const std::string &filename, const std::string &maname) : input(make_shared<Mesh>())
{
  openmeshfile(input.get(), this, filename, maname);
}

void QMAT::simplify_slab(int n)
{
  simplifySlab(this, input.get(), n);
}

void QMAT::export_ply(const std::string &fname)
{
  ExportPly(fname, input.get());
}

void QMAT::export_ma(const std::string &fname)
{
  Export(fname, input.get());
}

std::vector<double> QMAT::hausdorff()
{
  int nv = input->pVertexList.size();
  std::vector<double> hausdoff(nv);
  ComputeHausdorffDistance(*this, *input);
  for (int i = 0; i < nv; i++)
  {
    auto it = input->pVertexList[i];
    hausdoff[i] = it->slab_hausdorff_dist;
  }
  return hausdoff;
}

QMAT::QMAT(const QMAT &a) : input(a.input), SlabMesh(a) {}
double QMAT::get_diagonal()
{
  return input->bb_diagonal_length;
}


QMAT::QMAT(const QMAT &a, const std::string &maname) : input(a.input)
{
  initialize_slab(input.get(), this, maname);
}