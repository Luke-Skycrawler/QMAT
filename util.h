#pragma once
#include "src/NonManifoldMesh.h"
#include "src/SlabMesh.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include <string>

void ComputeHausdorffDistance(SlabMesh &slab_mesh, Mesh &input);
void LoadInputNMM(Mesh *input, SlabMesh *slabMesh, std::string maname, bool init_merged_list = true);
bool importMA(Mesh *input, SlabMesh *slabMesh, std::string maname);
void LoadSlabMesh(SlabMesh *slabMesh);
void openmeshfile(Mesh *input, SlabMesh *slabMesh, std::string filename, std::string maname);
void initialize_slab(Mesh *input, SlabMesh *slabMesh, std::string maname);
void simplifySlab(SlabMesh *slabMesh, Mesh *mesh, unsigned num_spheres);
