#pragma once
#include <tuple>
#include <vector>
#include "util.h"
#include <memory>

struct QMAT : public SlabMesh
{
    QMAT(const std::string &file, const std::string &maname);
    QMAT(const QMAT &a);
    QMAT(const QMAT &a, const std::string &maname);
    void simplify_slab(int n);
    std::shared_ptr<Mesh> input;
    std::vector<double> hausdorff();
    void export_ply(const std::string &fname);
    void export_ma(const std::string &fname);
    double get_diagonal();
};
