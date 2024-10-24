#pragma once
#include <iostream>
#include <vector>
#include "src/SlabMesh.h"
#include <set>

#include "subgraph.h"
#include <assert.h>

using uint = unsigned;
struct PointAdder
{
    std::vector<set<unsigned>> collapsed_list;
    std::vector<unsigned> included_in;
    double diagonal;
    SlabMesh &fine, coarse;
    PointAdder(double diag, SlabMesh &coarse, SlabMesh &fine) : diagonal(diag), coarse(coarse), fine(fine) {}

    double eval(set<uint> &merged, Sphere &s);
    double evaluate(uint t, SubGraph &G, uint x, Sphere s_x, Sphere s_sigma);
    uint nearest_node(Vector3d p, SlabMesh &slabmesh);
    void set_collapsed_list(const std::vector<std::vector<unsigned>> &merged_list);

    void generate_collapsed_list();
    void add_new_node(Sphere &new_sphere);
    void export_ply(const std::string &filename);
};