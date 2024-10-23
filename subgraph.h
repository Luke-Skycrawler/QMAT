#include <vector>
#include <set>
#include "src/SlabMesh.h"
#include <map>

using uint = unsigned;
struct SubGraph
{
    std::vector<uint> faces;
    std::set<std::pair<uint, uint>> edges;
    std::set<uint> vertices_set, taps;
    std::map<uint, std::vector<uint>> adjacency;

    bool critical(uint t);
    bool connected(uint u, uint v, uint r);
    uint split(uint mx, uint coase_id, std::vector<uint> &included_in, std::set<uint> &collapsed_list_x, std::set<uint> &collapsed_list_c);
    SubGraph(std::set<uint> &merged_list, SlabMesh &fine);
    std::vector<uint> non_critical_taps();
    inline void make_tap(uint u)
    {
        taps.insert(u);
    }
    inline void make_edge(int u, int v)
    {
        int sum = u + v;
        edges.insert({min(u, v), sum - min(u, v)});
        adjacency[u].push_back(v);
        adjacency[v].push_back(u);
    }
    inline void make_face(int u, int v, int w)
    {
        faces.push_back(u);
        faces.push_back(v);
        faces.push_back(w);
        adjacency[u].push_back(v);
        adjacency[u].push_back(w);
        adjacency[v].push_back(u);
        adjacency[v].push_back(w);
        adjacency[w].push_back(u);
        adjacency[w].push_back(v);
    }
};