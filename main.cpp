#include <iostream>
#include "subgraph.h"
#include "point_adder.h"
#include <assert.h>
#include <set>
#include "util.h"
using namespace std;
namespace ps = polyscope;

std::vector<std::vector<unsigned>> test_reentrant()
{
  string input_name = "../data/spider.off";
  // string fine = "../data/spider_v100.ma";
  string fine = "../data/spider.ma";
  Mesh input;
  SlabMesh slabmesh;
  openmeshfile(&input, &slabmesh, input_name, fine);
  simplifySlab(&slabmesh, &input, 100);
  slabmesh.ExportPly("../output/spider_to100", &input);

  slabmesh.initMergeList();
  slabmesh.initCollapseQueue();
  slabmesh.Simplify(75);

  slabmesh.ExportPly("../output/spider_100to25", &input);
  std::vector<std::vector<unsigned>> L;
  for (int i = 0; i < slabmesh.vertices.size(); i++)
  {
    if (slabmesh.vertices[i].first)
    {
      auto &list{slabmesh.vertices[i].second->merged_vertices};
      cout << i << ": { ";
      for (auto j : list)
      {
        cout << j << ",";
      }
      cout << "} " << endl;
      L.push_back(list);
    }
  }
  return L;
}
void test_add_sphere()
{
  string input_name = "../data/spider.off";
  //string fine = "../data/spider_v100.ma";
  string fine = "../data/spider_v100.ma";
  string coarse = "../data/spider_v25.ma";
  Mesh input;
  SlabMesh slab_coarse, slab_fine;
  openmeshfile(&input, &slab_fine, input_name, fine);
  simplifySlab(&slab_fine, &input, 100);
  slab_fine.AdjustStorage();

  initialize_slab(&input, &slab_coarse, coarse);
  Sphere new_sphere = Sphere{Vector3d(-0.2176294, 0.06370435, 0.09288197), 0.04448237};
  PointAdder adder(input.bb_diagonal_length, slab_coarse, slab_fine);
  // adder.generate_collapsed_list();
  adder.set_collapsed_list(test_reentrant());
  adder.add_new_node(new_sphere);
  adder.export_ply("../output/spider_v26");
}

void test_add_node2()
{
  string input_name = "../data/spider.off";
  string fine = "../data/spider.ma";
  Mesh input;
  SlabMesh fine_mesh;
  openmeshfile(&input, &fine_mesh, input_name, fine);
  simplifySlab(&fine_mesh, &input, 100);
  fine_mesh.ExportPly("../output/spider_to100", &input);

  SlabMesh coarse_mesh = fine_mesh;

  coarse_mesh.initMergeList();
  coarse_mesh.initCollapseQueue();
  coarse_mesh.Simplify(75);

  coarse_mesh.ExportPly("../output/spider_100to25", &input);
  std::vector<std::vector<unsigned>> L;
  for (int i = 0; i < coarse_mesh.vertices.size(); i++)
  {
    if (coarse_mesh.vertices[i].first)
    {
      auto &list{coarse_mesh.vertices[i].second->merged_vertices};
      cout << i << ": { ";
      for (auto j : list)
      {
        cout << j << ",";
      }
      cout << "} " << endl;
      L.push_back(list);
    }
  }
  Sphere new_sphere = Sphere{Vector3d(-0.2176294, 0.06370435, 0.09288197), 0.04448237};

  PointAdder adder(input.bb_diagonal_length, coarse_mesh, fine_mesh);
  adder.set_collapsed_list(L);
  adder.add_new_node(new_sphere);
  adder.export_ply("../output/spider_v26");
}

int main()
{
  // test_reentrant();
   test_add_sphere();
  //test_add_node2();
  return 0;
}
// int main(int argc, char** argv) {
//   if (4 > argc) {
//     std::cerr << "Usage: " << argv[0]
//               << " <surface_mesh.off> <medial_mesh.ma> <num_target_spheres>"
//               << std::endl;
//     return 1;
//   }
//   std::string filename = argv[1];
//   std::string maname = argv[2];
//   unsigned num_spheres = atoi(argv[3]);
//   printf("reading off file %s\n", filename.c_str());

//   Mesh input;
//   Mesh* pinput = &input;
//   SlabMesh slabMesh;
//   SlabMesh* pslabMesh = &slabMesh;
//   openmeshfile(pinput, pslabMesh, filename, maname);
//   printf("done openmeshfile\n");
//   simplifySlab(pslabMesh, pinput, num_spheres);
//   printf("done simplifyslab\n");
//   pslabMesh->ExportPly("export_half", pinput);
//   pslabMesh->Export("export_half", pinput);
//   printf("done export\n");

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

//   // V = slabMesh.V();
//   // V *= 2.0;
//   // slabMesh.displace_vertices(V);
//   // slabMesh.Export("2x", pinput);
//   ps::show();

//   return 0;
// }
