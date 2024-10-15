#include <string>
#include "src/NonManifoldMesh.h"
#include "src/SlabMesh.h"
#include <tuple>
#include <vector>
struct QMAT{
    QMAT(const std::string &filename, const std::string& maname);

    void ComputeHausdorffDistance();
    void simplifySlab(unsigned num_spheres);
    SlabMesh slab_mesh;
    Mesh input; 
    void ExportPly(const std::string & fname);
    void ExportMA(const std::string & fname);
    std::vector<double> export_hausdorff_distance();
};
