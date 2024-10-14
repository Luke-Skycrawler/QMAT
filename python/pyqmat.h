#include <string>
#include "src/NonManifoldMesh.h"
#include "src/SlabMesh.h"
struct QMAT{
    QMAT(const std::string &filename, const std::string& maname);

    void ComputeHausdorffDistance();
    void simplifySlab(unsigned num_spheres);
    SlabMesh slab_mesh;
    Mesh input; 
};
