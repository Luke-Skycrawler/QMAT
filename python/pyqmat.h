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

struct QMATinh : public SlabMesh
{
    QMATinh(const std::string &file, const std::string &maname);
    void simplify(int n);
    Mesh input;
    std::vector<double> hausdorff();
    void clean_up();
    void export_ply(const std::string &fname);
    void export_ma(const std::string &fname);
};
