#include <iostream>
#include <vector>
#include "src/NonManifoldMesh.h"
#include "src/SlabMesh.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "subgraph.h"
#include <assert.h>
#include <set>
using namespace std;
namespace ps = polyscope;

void ComputeHausdorffDistance(SlabMesh &slab_mesh, Mesh &input)
{

  // ma_qem_mesh.maxhausdorff_distance = 0.;
  slab_mesh.maxhausdorff_distance = 0;
  double sumhausdorff_distance = 0;
  for (unsigned i = 0; i < input.pVertexList.size(); i++)
  {
    double min_dis = DBL_MAX;
    unsigned min_index = -1;
    Vector3d bou_ver(input.pVertexList[i]->point()[0], input.pVertexList[i]->point()[1], input.pVertexList[i]->point()[2]);
    bou_ver /= input.bb_diagonal_length;

    for (unsigned j = 0; j < slab_mesh.numVertices; j++)
    {
      Sphere ma_ver = slab_mesh.vertices[j].second->sphere;
      double temp_length = abs((bou_ver - ma_ver.center).Length() - ma_ver.radius);
      // if (temp_length >= 0 && temp_length < min_dis)
      if (temp_length < min_dis)
      {
        min_dis = temp_length;
        min_index = j;
      }
    }

    if (min_index != -1)
    {
      double temp_near_dis = slab_mesh.NearestPoint(bou_ver, min_index);
      min_dis = min(temp_near_dis, min_dis);

      sumhausdorff_distance += min_dis;

      // ma_qem_mesh.vertices[min_index].second->bplist.push_back(i);
      // ma_qem_mesh.maxhausdorff_distance = max(ma_qem_mesh.maxhausdorff_distance, min_dis);

      slab_mesh.vertices[min_index].second->bplist.insert(i);
      slab_mesh.maxhausdorff_distance = max(slab_mesh.maxhausdorff_distance, min_dis);

      // input.pVertexList[i]->vqem_hausdorff_dist = min_dis / input.bb_diagonal_length;
      input.pVertexList[i]->vqem_hausdorff_dist = min_dis;
      input.pVertexList[i]->vqem_hansdorff_index = min_index;

      // input.pVertexList[i]->slab_hausdorff_dist = min_dis / input.bb_diagonal_length;
      input.pVertexList[i]->slab_hausdorff_dist = min_dis;
      input.pVertexList[i]->slab_hansdorff_index = min_index;
    }
  }

  // ma_qem_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
  slab_mesh.meanhausdorff_distance = sumhausdorff_distance / input.pVertexList.size();
}

void LoadInputNMM(Mesh *input, SlabMesh *slabMesh, std::string maname, bool init_merged_list = true)
{
  std::ifstream mastream(maname.c_str());
  NonManifoldMesh newinputnmm;
  newinputnmm.numVertices = 0;
  newinputnmm.numEdges = 0;
  newinputnmm.numFaces = 0;
  int nv, ne, nf;
  mastream >> nv >> ne >> nf;

  // slab mesh
  slabMesh->numVertices = 0;
  slabMesh->numEdges = 0;
  slabMesh->numFaces = 0;
  slabMesh->bound_weight = 0.1;

  double len[4];
  len[0] = input->m_max[0] - input->m_min[0];
  len[1] = input->m_max[1] - input->m_min[1];
  len[2] = input->m_max[2] - input->m_min[2];
  len[3] = sqrt(len[0] * len[0] + len[1] * len[1] + len[2] * len[2]);
  newinputnmm.diameter = len[3];

  for (unsigned i = 0; i < input->pVertexList.size(); i++)
    newinputnmm.BoundaryPoints.push_back(SamplePoint(
        input->pVertexList[i]->point()[0], input->pVertexList[i]->point()[1],
        input->pVertexList[i]->point()[2]));

  for (unsigned i = 0; i < nv; i++)
  {
    char ch;
    double x, y, z, r;
    mastream >> ch >> x >> y >> z >> r;

    // handle the slab mesh
    Bool_SlabVertexPointer bsvp2;
    bsvp2.first = true;
    bsvp2.second = new SlabVertex;
    (*bsvp2.second).sphere.center[0] = x / input->bb_diagonal_length;
    (*bsvp2.second).sphere.center[1] = y / input->bb_diagonal_length;
    (*bsvp2.second).sphere.center[2] = z / input->bb_diagonal_length;
    (*bsvp2.second).sphere.radius = r / input->bb_diagonal_length;
    unsigned index = slabMesh->vertices.size();
    (*bsvp2.second).index = index;
    if (init_merged_list)
      (*bsvp2.second).merged_vertices = {index};
    else
    {
      (*bsvp2.second).merged_vertices = {};
    }
    slabMesh->vertices.push_back(bsvp2);
    slabMesh->numVertices++;
  }

  for (unsigned i = 0; i < ne; i++)
  {
    char ch;
    unsigned ver[2];
    mastream >> ch;
    mastream >> ver[0];
    mastream >> ver[1];

    // handle the slab mesh
    Bool_SlabEdgePointer bsep2;
    bsep2.first = true;
    bsep2.second = new SlabEdge;
    (*bsep2.second).vertices_.first = ver[0];
    (*bsep2.second).vertices_.second = ver[1];
    (*slabMesh->vertices[(*bsep2.second).vertices_.first].second)
        .edges_.insert(slabMesh->edges.size());
    (*slabMesh->vertices[(*bsep2.second).vertices_.second].second)
        .edges_.insert(slabMesh->edges.size());
    (*bsep2.second).index = slabMesh->edges.size();
    slabMesh->edges.push_back(bsep2);
    slabMesh->numEdges++;
  }

  for (unsigned i = 0; i < nf; i++)
  {
    char ch;
    unsigned vid[3];
    unsigned eid[3];
    mastream >> ch >> vid[0] >> vid[1] >> vid[2];

    // handle the slab mesh
    Bool_SlabFacePointer bsfp2;
    bsfp2.first = true;
    bsfp2.second = new SlabFace;
    (*bsfp2.second).vertices_.insert(vid[0]);
    (*bsfp2.second).vertices_.insert(vid[1]);
    (*bsfp2.second).vertices_.insert(vid[2]);
    if (slabMesh->Edge(vid[0], vid[1], eid[0]))
      (*bsfp2.second).edges_.insert(eid[0]);
    if (slabMesh->Edge(vid[0], vid[2], eid[1]))
      (*bsfp2.second).edges_.insert(eid[1]);
    if (slabMesh->Edge(vid[1], vid[2], eid[2]))
      (*bsfp2.second).edges_.insert(eid[2]);
    (*bsfp2.second).index = slabMesh->faces.size();
    slabMesh->vertices[vid[0]].second->faces_.insert(slabMesh->faces.size());
    // slab_mesh.vertices[vid[0]].second->related_face += 2;
    slabMesh->vertices[vid[1]].second->faces_.insert(slabMesh->faces.size());
    // slab_mesh.vertices[vid[1]].second->related_face += 2;
    slabMesh->vertices[vid[2]].second->faces_.insert(slabMesh->faces.size());
    // slab_mesh.vertices[vid[2]].second->related_face += 2;
    slabMesh->edges[eid[0]].second->faces_.insert(slabMesh->faces.size());
    slabMesh->edges[eid[1]].second->faces_.insert(slabMesh->faces.size());
    slabMesh->edges[eid[2]].second->faces_.insert(slabMesh->faces.size());
    slabMesh->faces.push_back(bsfp2);
    slabMesh->numFaces++;
  }

  // newinputnmm.ComputeFacesNormal();
  // newinputnmm.ComputeFacesCentroid();
  // newinputnmm.ComputeFacesSimpleTriangles();
  // newinputnmm.ComputeEdgesCone();
  // input_nmm = newinputnmm;

  slabMesh->iniNumVertices = slabMesh->numVertices;
  slabMesh->iniNumEdges = slabMesh->numEdges;
  slabMesh->iniNumFaces = slabMesh->numFaces;

  slabMesh->CleanIsolatedVertices();
  slabMesh->computebb();
  slabMesh->ComputeFacesCentroid();
  slabMesh->ComputeFacesNormal();
  slabMesh->ComputeVerticesNormal();
  slabMesh->ComputeEdgesCone();
  slabMesh->ComputeFacesSimpleTriangles();
  slabMesh->DistinguishVertexType();
}

bool importMA(Mesh *input, SlabMesh *slabMesh, std::string maname)
{
  // std::string filename = filename + ".ma";
  // if (!std::filesystem::exists(filename)) {
  //     std::cerr << "Related .ma file is missing." << std::endl;
  //     return false;
  // }

  std::cout << "Loading ma file " << maname << std::endl;
  // bool success = false;

  // m_pThreeDimensionalShape->input_nmm.meshname = filename;
  // m_pThreeDimensionalShape->input_nmm.domain =
  // m_pThreeDimensionalShape->input.domain;
  // m_pThreeDimensionalShape->input_nmm.pmesh =
  // &(m_pThreeDimensionalShape->input);
  // m_pThreeDimensionalShape->slab_mesh.pmesh =
  // &(m_pThreeDimensionalShape->input);
  slabMesh->type = 1;
  slabMesh->bound_weight = 1.0;
  // m_pThreeDimensionalShape->slab_mesh.type = 1;
  // m_pThreeDimensionalShape->slab_mesh.bound_weight = 1.0;
  // m_pThreeDimensionalShape->LoadInputNMM(filename);
  LoadInputNMM(input, slabMesh, maname);
  std::cout << "import MA done." << std::endl;
  // success = true;

  // if (success) {
  //     m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
  // }

  return true;
}

void LoadSlabMesh(SlabMesh *slabMesh)
{
  slabMesh->clear();
  // long startt = clock();
  // handle each face
  for (unsigned i = 0; i < slabMesh->vertices.size(); i++)
  {
    if (!slabMesh->vertices[i].first)
      continue;

    SlabVertex sv = *slabMesh->vertices[i].second;
    std::set<unsigned> fset = sv.faces_;
    Vector4d C1(sv.sphere.center.X(), sv.sphere.center.Y(),
                sv.sphere.center.Z(), sv.sphere.radius);

    for (set<unsigned>::iterator si = fset.begin(); si != fset.end(); si++)
    {
      SlabFace sf = *slabMesh->faces[*si].second;

      if (sf.valid_st == false || sf.st[0].normal == Vector3d(0., 0., 0.) ||
          sf.st[1].normal == Vector3d(0., 0., 0.))
        continue;

      Vector4d normal1(sf.st[0].normal.X(), sf.st[0].normal.Y(),
                       sf.st[0].normal.Z(), 1.0);
      Vector4d normal2(sf.st[1].normal.X(), sf.st[1].normal.Y(),
                       sf.st[1].normal.Z(), 1.0);

      // compute the matrix of A
      Matrix4d temp_A1, temp_A2;
      temp_A1.MakeTensorProduct(normal1, normal1);
      temp_A2.MakeTensorProduct(normal2, normal2);
      temp_A1 *= 2.0;
      temp_A2 *= 2.0;

      // compute the matrix of b
      double normal_mul_point1 = normal1.Dot(C1);
      double normal_mul_point2 = normal2.Dot(C1);
      Wm4::Vector4d temp_b1 = normal1 * 2 * normal_mul_point1;
      Wm4::Vector4d temp_b2 = normal2 * 2 * normal_mul_point2;

      // compute c
      double temp_c1 = normal_mul_point1 * normal_mul_point1;
      double temp_c2 = normal_mul_point2 * normal_mul_point2;

      slabMesh->vertices[i].second->slab_A += temp_A1;
      slabMesh->vertices[i].second->slab_A += temp_A2;
      slabMesh->vertices[i].second->slab_b += temp_b1;
      slabMesh->vertices[i].second->slab_b += temp_b2;
      slabMesh->vertices[i].second->slab_c += temp_c1;
      slabMesh->vertices[i].second->slab_c += temp_c2;

      slabMesh->vertices[i].second->related_face += 2;
    }
  }

  switch (slabMesh->preserve_boundary_method)
  {
  case 1:
    slabMesh->PreservBoundaryMethodOne();
    break;
  case 2:
    // slab_mesh.PreservBoundaryMethodTwo();
    break;
  case 3:
    slabMesh->PreservBoundaryMethodThree();
    break;
  default:
    slabMesh->PreservBoundaryMethodFour();
    break;
  }

  slabMesh->initCollapseQueue();
  // long endt = clock();
}

void openmeshfile(Mesh *input, SlabMesh *slabMesh, std::string filename,
                  std::string maname)
{
  // QString filename = QFileDialog::getOpenFileName(this, tr("Select a 3D model
  // to open"), NULL, tr("3D model(*.off)"));
  if (!filename.empty())
  {
    // QDir qd(filename);
    std::string prefix = filename.substr(0, filename.size() - 4);
    // ThreeDimensionalShape * pThreeDimensionalShape = new
    // ThreeDimensionalShape;
    bool suc = false;
    std::ifstream stream(filename);

    if (stream)
    {
      stream >> *input;
      // compute the properties of the input mesh
      input->computebb();           // bounding box
      input->GenerateList();        // generate vertex and triangle list
      input->GenerateRandomColor(); // color of vertex and triangle
      input->compute_normals();     // normal of vertex and triangle
      // pThreeDimensionalShape->input_nmm.meshname = filename.;
      cout << input->pVertexList.size() << " " << input->pFaceList.size() << endl;
      std::ifstream streampol(filename);
      Polyhedron pol;
      streampol >> pol;

      // set the non manifold mesh
      // Mesh_domain * pdom;
      // pdom = new Mesh_domain(pol);
      suc = true;
    }
    if (suc)
    {
      // m_pThreeDimensionalShape = pThreeDimensionalShape;
      // //ui.actionShow_Edge->setChecked(true);
      //
      // ui.actionShow_Face->setChecked(true);
      // ui.actionReverse_Orientation->setChecked(true);

      // importVP(prefix);
      slabMesh->pmesh = input;
      bool re = importMA(input, slabMesh, maname);
      if (re == false)
        return;

      float k = 0.00001;
      slabMesh->k = k;
      // initialize();
      slabMesh->preserve_boundary_method = 0;
      slabMesh->hyperbolic_weight_type = 3;
      slabMesh->compute_hausdorff = true;
      slabMesh->boundary_compute_scale = 0;
      slabMesh->prevent_inversion = false;

      LoadSlabMesh(slabMesh);
      // long ti = m_pThreeDimensionalShape->LoadSlabMesh();
      // slab_initial = true;

      // // set back
      // ui.actionSet_k_value->setChecked(false);
      //
      // // show the input mesh in the dialog.
      // m_pGLWidget->set3DShape(m_pThreeDimensionalShape);
      // statusBar()->showMessage(filename + tr(" is loaded successfully.") );
      // setWindowTitle( tr("Medial Axis Simplification 3D - ") + filename );

      // m_isSimplified = false;
      std::cout << "openmeshfile done." << std::endl;
    }
    else
    {
    }
  }
  else
  {
    std::cout << "Filename is empty !" << std::endl;
  }
}

void initialize_slab(Mesh *input, SlabMesh *slabMesh, std::string maname)
{
  slabMesh->pmesh = input;
  bool re = importMA(input, slabMesh, maname);
  if (re == false)
    return;

  float k = 0.00001;
  slabMesh->k = k;
  // initialize();
  slabMesh->preserve_boundary_method = 0;
  slabMesh->hyperbolic_weight_type = 3;
  slabMesh->compute_hausdorff = true;
  slabMesh->boundary_compute_scale = 0;
  slabMesh->prevent_inversion = false;

  LoadSlabMesh(slabMesh);
}
void simplifySlab(SlabMesh *slabMesh, Mesh *mesh, unsigned num_spheres)
{
  slabMesh->CleanIsolatedVertices();
  // int threhold = min(10000, (int)(slabMesh->numVertices / 2));
  int threhold = num_spheres;

  // bool ok = true;
  // int simplifyNum = min(10000, (int)(slabMesh->numVertices / 2));
  // if (ok)
  //     threhold = simplifyNum;
  // else
  //     return;

  // if(slabMesh == NULL)
  //     return;

  // long start_time = clock();
  // slabMesh->initCollapseQueue();
  // slabMesh->initBoundaryCollapseQueue();
  slabMesh->compute_hausdorff = true;
  slabMesh->Simplify(slabMesh->numVertices - threhold);
  // long end_time = clock();
  //
  // std::string res;
  // std::stringstream ss;
  // ss << end_time - start_time;
  // ss >> res;

  // slabMesh->ComputeFacesNormal();
  // slabMesh->ComputeVerticesNormal();
  // slabMesh->ComputeEdgesCone();
  // slabMesh->ComputeFacesSimpleTriangles();

  std::cout << "Simplify done." << std::endl;
}

struct PointAdder
{
  std::vector<set<unsigned>> collapsed_list;
  std::vector<unsigned> included_in;
  double diagonal;
  SlabMesh &fine, coarse;
  PointAdder(double diag, SlabMesh &coarse, SlabMesh &fine) : diagonal(diag), coarse(coarse), fine(fine) {}
  uint nearest_node(Vector3d p, SlabMesh &slabmesh)
  {
    uint min_idx = 0;
    double minl = DBL_MAX;
    for (unsigned i = 0; i < slabmesh.numVertices; i++)
    {
      Vector3d pi = slabmesh.vertices[i].second->sphere.center;
      double l = (pi - p).Length();
      if (l < minl)
      {
        minl = l;
        min_idx = i;
      }
    }
    return min_idx;
  }

  void set_collapsed_list(const std::vector<std::vector<unsigned>> &merged_list)
  {
    collapsed_list.resize(coarse.numVertices);
    included_in.resize(fine.numVertices);
    for (int i = 0; i < coarse.numVertices; i++)
    {
      for (auto j : merged_list[i])
      {

        collapsed_list[i].insert(j);
        included_in[j] = i;
      }
    }
  }
  void generate_collapsed_list()
  {
    collapsed_list.resize(coarse.numVertices);
    included_in.resize(fine.numVertices);

    for (int i = 0; i < fine.numVertices; i++)
    {
      Vector3d pi = fine.vertices[i].second->sphere.center;
      // int min_idx = -1;
      // double minl = DBL_MAX;
      int min_idx = nearest_node(pi, coarse);
      collapsed_list[min_idx].insert(i);
      included_in[i] = min_idx;
    }
    for (int i = 0; i < coarse.numVertices; i++)
    {
      assert(!collapsed_list[i].empty());
    }
  }

  void add_new_node(Sphere &new_sphere)
  {
    new_sphere.center /= diagonal;
    new_sphere.radius /= diagonal;
    uint ik = nearest_node(new_sphere.center, fine);
    int sigma = included_in[ik];

    cout << "ik = " << ik << "sigma = " << sigma;
    assert(collapsed_list[sigma].size() > 1);
    if (!collapsed_list[sigma].size() > 1)
    {
      exit(1);
    }

    Sphere s_sigma = coarse.vertices[sigma].second->sphere;
    SubGraph c_sg(collapsed_list[sigma], fine);
    uint i_sigma = nearest_node(s_sigma.center, fine);
    int id = coarse.numVertices++;
    collapsed_list.push_back({});
    if (c_sg.taps.count(ik) && c_sg.critical(ik))
    {
      included_in[ik] = id;
      collapsed_list[id] = {ik};
      collapsed_list[sigma].erase(ik);
    }
    else if (c_sg.taps.count(i_sigma) && c_sg.critical(i_sigma))
    {
      collapsed_list[id] = collapsed_list[sigma];
      collapsed_list[id].erase(i_sigma);
      for (auto i : collapsed_list[sigma])
      {
        included_in[i] = id;
      }
      collapsed_list[sigma] = {i_sigma};
    }
    else
    {
      c_sg.split(ik, id, included_in, collapsed_list[id], collapsed_list[sigma]);
    }
    // collapsed_list[sigma].erase(ik);

    // included_in[ik] = id;
    Bool_SlabVertexPointer bsvp;
    bsvp.first = true;
    bsvp.second = new SlabVertex;
    bsvp.second->sphere.center = new_sphere.center;
    bsvp.second->sphere.radius = new_sphere.radius;
    bsvp.second->index = id;
    coarse.vertices.push_back(bsvp);
  }

  void export_ply(const string &filename)
  {
    string ply_name = filename + ".ply";
    std::ofstream fout(ply_name);

    int tot_prims = 0;
    for (int i = 0; i < fine.edges.size(); i++)
    {
      int uu = fine.edges[i].second->vertices_.first;
      int vv = fine.edges[i].second->vertices_.second;
      int u = included_in[uu];
      int v = included_in[vv];
      if (u != v)
      {
        // fout << "2 " << u << " " << v << std::endl;
        tot_prims++;
      }
    }
    for (int i = 0; i < fine.faces.size(); i++)
    {
      int idx[3];

      int j = 0;
      for (auto it = fine.faces[i].second->vertices_.begin(); it != fine.faces[i].second->vertices_.end(); it++)
      {
        idx[j++] = *it;
      }

      int u = included_in[idx[0]];
      int v = included_in[idx[1]];
      int w = included_in[idx[2]];
      if (u == v && v == w)
      {
        continue;
      }
      else if (u != v && v != w && w != u)
      {
        // fout << "3 " << u << " " << v << " " << w << std::endl;
        tot_prims++;
      }
      else
      {
        if (u == v)
        {
          // fout << "2" << u << " " << w << std::endl;
          tot_prims++;
        }
        if (v == w)
        {
          // fout << "2" << u << " " << v << std::endl;
          tot_prims++;
        }
        if (u == w)
        {
          // fout << "2" << u << " " << v << std::endl;
          tot_prims++;
        }
      }
    }

    fout << "ply" << std::endl;
    fout << "format ascii 1.0" << std::endl;
    fout << "element vertex " << coarse.numVertices << std::endl;
    fout << "property float x" << endl;
    fout << "property float y" << endl;
    fout << "property float z" << endl;
    fout << "property float r" << endl;
    fout << "element face " << tot_prims << endl;
    fout << "property list uchar uint vertex_indices" << endl;
    fout << "end_header" << endl;

    for (unsigned i = 0; i < coarse.numVertices; i++)
    {
      auto &vertices{coarse.vertices};
      fout << setiosflags(ios::fixed) << setprecision(15) << (vertices[i].second->sphere.center * diagonal) << " " << (vertices[i].second->sphere.radius * diagonal) << std::endl;
    }

    for (int i = 0; i < fine.edges.size(); i++)
    {
      int uu = fine.edges[i].second->vertices_.first;
      int vv = fine.edges[i].second->vertices_.second;
      int u = included_in[uu];
      int v = included_in[vv];
      if (u != v)
      {
        fout << "2 " << u << " " << v << std::endl;
      }
    }
    for (int i = 0; i < fine.faces.size(); i++)
    {
      int idx[3];

      int j = 0;
      for (auto it = fine.faces[i].second->vertices_.begin(); it != fine.faces[i].second->vertices_.end(); it++)
      {
        idx[j++] = *it;
      }

      int u = included_in[idx[0]];
      int v = included_in[idx[1]];
      int w = included_in[idx[2]];
      if (u == v && v == w)
      {
        continue;
      }
      else if (u != v && v != w && w != u)
      {
        fout << "3 " << u << " " << v << " " << w << std::endl;
      }
      else
      {
        if (u == v)
        {
          fout << "2" << u << " " << w << std::endl;
        }
        if (v == w)
        {
          fout << "2" << u << " " << v << std::endl;
        }
        if (u == w)
        {
          fout << "2" << u << " " << v << std::endl;
        }
      }
    }
  }
};

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
  string fine = "../data/spider_v100.ma";
  string coarse = "../data/spider_v25.ma";
  Mesh input;
  SlabMesh slab_coarse, slab_fine;
  openmeshfile(&input, &slab_fine, input_name, fine);

  initialize_slab(&input, &slab_coarse, coarse);
  Sphere new_sphere = Sphere{Vector3d(-0.2176294, 0.06370435, 0.09288197), 0.04448237};
  PointAdder adder(input.bb_diagonal_length, slab_coarse, slab_fine);
  // adder.generate_collapsed_list();
  adder.set_collapsed_list(test_reentrant());
  adder.add_new_node(new_sphere);
  adder.export_ply("../output/spider_v26");
}

int main()
{
  // test_reentrant();
  test_add_sphere();
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
