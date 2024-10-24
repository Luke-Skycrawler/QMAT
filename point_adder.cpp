#include "point_adder.h"
using namespace std;

double PointAdder::eval(set<uint> &merged, Sphere &s)
{
    Matrix4d A;
    Vector4d b(0, 0, 0, 0), r;
    double c;
    A.MakeZero();
    r = Vector4d(s.center.X(), s.center.Y(), s.center.Z(), s.radius);
    for (auto i : merged)
    {
        auto &v{fine.vertices[i].second};
        A += v->slab_A + v->add_A;
        b += v->slab_b + v->add_b;
        c += v->slab_c + v->add_c;
    }
    // 1/2 xT A x - bT x + c
    double cost = 0.5 * (r * A).Dot(r) - b.Dot(r) + c;
    return cost;
}
// evaluate the cost of splitting subgraph G to {t} and G\{t}
double PointAdder::evaluate(uint t, SubGraph &G, uint x, Sphere s_x, Sphere s_sigma)
{
    auto set1 = G.vertices_set;
    set1.erase(t);
    set<uint> set2{t};
    Sphere s1, s2;
    if (t == x)
    {
        s1 = s_sigma;
        s2 = s_x;
    }
    else
    {
        s1 = s_x;
        s2 = s_sigma;
    }
    return eval(set1, s1) + eval(set2, s2);
}
uint PointAdder::nearest_node(Vector3d p, SlabMesh &slabmesh)
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

void PointAdder::set_collapsed_list(const std::vector<std::vector<unsigned>> &merged_list)
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
void PointAdder::generate_collapsed_list()
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

void PointAdder::add_new_node(Sphere &new_sphere)
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
    SubGraph G(collapsed_list[sigma], fine);
    uint i_sigma = nearest_node(s_sigma.center, fine);
    int id = coarse.numVertices++;
    collapsed_list.push_back({});

    // G.split(ik, id, included_in, collapsed_list[id], collapsed_list[sigma]);
    auto candidates = G.non_critical_taps();
    double min_cost = DBL_MAX;
    uint min_t = 0;
    for (auto t : candidates)
    {
        double cost = evaluate(t, G, ik, new_sphere, s_sigma);
        if (cost < min_cost)
        {
            min_cost = cost;
            min_t = t;
        }
    }

    if (min_t == ik)
    {
        included_in[ik] = id;
        collapsed_list[id] = {ik};
        collapsed_list[sigma].erase(ik);
    }
    else
    {
        collapsed_list[id] = collapsed_list[sigma];
        collapsed_list[id].erase(min_t);
        for (auto i : collapsed_list[id])
        {
            included_in[i] = id;
        }
        collapsed_list[sigma] = {min_t};
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

void PointAdder::export_ply(const string &filename)
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
