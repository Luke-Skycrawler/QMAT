#pragma once
#include <tuple>
#include <vector>
#include "util.h"

struct QMAT : public SlabMesh
{
    QMAT(const std::string &file, const std::string &maname);
    void simplify_slab(int n);
    Mesh input;
    std::vector<double> hausdorff();
    void export_ply(const std::string &fname);
    void export_ma(const std::string &fname);
};
