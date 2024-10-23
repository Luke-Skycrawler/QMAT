#include "subgraph.h"
#include <queue>
#include <assert.h>
using namespace std;

// construct merged subgraph for coarse mesh vertex u
// from merged vertices list
SubGraph::SubGraph(set<uint> &merged, SlabMesh &fine) : vertices_set(merged)
{
    std::set<uint> vset;
    for (int i = 0; i < fine.edges.size(); i++)
    {
        auto &e{fine.edges[i].second->vertices_};
        int u = e.first;
        int v = e.second;

        bool has_u = vertices_set.find(u) != vertices_set.end();
        bool has_v = vertices_set.find(v) != vertices_set.end();
        if (has_u && has_v)
        {
            // inside the collapsed subgraph
            make_edge(u, v);
        }
        else
        {
            if (has_u)
            {
                make_tap(u);
            }
            if (has_v)
            {
                make_tap(v);
            }
        }
    }
    for (int i = 0; i < fine.faces.size(); i++)
    {
        int idx[3];
        auto &f{fine.faces[i].second->vertices_};
        int j = 0;
        for (auto it = f.begin(); it != f.end(); it++)
        {
            idx[j++] = *it;
        }
        bool has_u = vertices_set.find(idx[0]) != vertices_set.end();
        bool has_v = vertices_set.find(idx[1]) != vertices_set.end();
        bool has_w = vertices_set.find(idx[2]) != vertices_set.end();

        if (has_u && has_v && has_w)
        {
            make_face(idx[0], idx[1], idx[2]);
        }
        else
        {
            if (has_u)
            {
                make_tap(idx[0]);
            }
            if (has_v)
            {
                make_tap(idx[1]);
            }
            if (has_w)
            {
                make_tap(idx[2]);
            }
            if (has_u && has_v)
            {
                make_edge(idx[0], idx[1]);
            }
            if (has_v && has_w)
            {
                make_edge(idx[1], idx[2]);
            }
            if (has_w && has_u)
            {
                make_edge(idx[2], idx[0]);
            }
        }
    }
}

// test if the subgraph is connected after removing the tap vertice t
// by connected we mean there still exists a path between any pair of the rest of the taps
bool SubGraph::critical(uint t)
{
    if (taps.size() <= 2)
    {
        return false;
    }
    for (auto i : taps)
    {
        if (i == t)
            continue;
        for (auto j : taps)
        {
            if (j == t || j == i)
                continue;
            if (connected(i, j, t))
            {
            }
            else
            {
                return true;
            }
        }
    }
    return false;
}

// if it is able to connect tap u and v after removing tap r
bool SubGraph::connected(uint u, uint v, uint r)
{
    queue<uint> q;
    set<uint> visited;
    assert(u != v && u != r && v != r);
    assert(taps.count(u) && taps.count(v) && taps.count(r));
    q.push(u);
    while (!q.empty())
    {
        uint i = q.front();
        if (i == v)
        {
            return true;
        }
        if (i == r || visited.count(i))
        {
            continue;
        }

        // new vertex reached, enqueue all its neighbors
        visited.insert(i);
        for (auto j : adjacency[i])
        {
            if (!visited.count(j))
            {
                q.push(j);
            }
        }
    }
}

// given a point x to be split from collapsed set c
// split the subgraph into two parts that x in one part and c gets the rest

// @returns the splited point u
// if u == x, it means x is a sigle-point subgraph
// otherwise, c becomes a single-point subgraph {u}, and x gets the rest of the elements
// fixme: for the single point subgraph, it might has some dangling edges that are not connected to the outside as a tap. These vertices should also be added to the one-point subgraph
uint SubGraph::split(uint x, uint coarse_id, vector<uint> &included_in, set<uint> &merged_x, set<uint> &merged_c)
{
    uint ret = 1<<30;
    for (auto u : taps)
    {
        if (critical(u))
            continue;

        // u is non-critical, returns
        ret = u;
        break;
    }
    if (ret == x) {
        included_in[x] = coarse_id;
        merged_x = {x};
        merged_c.erase(x);
    }
    else {
        merged_x = merged_c;
        merged_x.erase(ret);
        for (auto i: merged_x) {
            included_in[i] = coarse_id;
        }
        merged_c = {x};
    }
    if (ret == 1 << 30) {
        // must exists one non-critical tap
        assert(false);
    }
    return 0;
}