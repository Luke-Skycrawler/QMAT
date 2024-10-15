#include "pyqmat.h"
// #include "polyscope/polyscope.h"
// #include "polyscope/surface_mesh.h"

using namespace std;
int main(int argc, char ** argv) {
  if (4 > argc) {
    std::cerr << "Usage: " << argv[0]
              << " <surface_mesh.off> <medial_mesh.ma> <num_target_spheres>"
              << std::endl;
    return 1;
  }
  std::string filename = argv[1];
  std::string maname = argv[2];
  unsigned num_spheres = atoi(argv[3]);
  printf("reading off file %s\n", filename.c_str());
    QMAT qmat(filename, maname);
    qmat.simplifySlab(num_spheres);
    qmat.ExportMA("export_half");
    qmat.ExportPly("export_half");
//   ps::init();
//   int nv = pinput ->pVertexList.size();
//   int nf = pinput -> pFaceList.size();
//   Eigen::MatrixXd V(nv, 3);
//   Eigen::MatrixXi F(nf, 3); 
//   std::vector<double> hausdoff(nv);
//   ComputeHausdorffDistance(slabMesh, input);
//   for (int i = 0; i < nv; i++) {
//     auto it = input.pVertexList[i];
//     V.row(i) = Eigen::Vector3d(it->point()[0], it->point()[1], it->point()[2]);
//     hausdoff[i] = it-> slab_hausdorff_dist;
//   }
//   for (int i = 0; i < nf; i++) {
//     auto it = input.pFaceList[i] -> facet_begin();
//     int i0 = it -> vertex() -> id;
//     int i1 = it -> next() -> vertex() -> id;
//     int i2 = it -> next() -> next() -> vertex() -> id;
//     F.row(i) = Eigen::Vector3i(i0, i1, i2);
//   }

//   auto input_mesh = ps::registerSurfaceMesh("input mesh", V, F);
//   input_mesh -> addVertexScalarQuantity("hausdorff", hausdoff);
//   ps::show();
  return 0;
}