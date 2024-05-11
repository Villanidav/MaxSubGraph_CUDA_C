#include <iostream>
#include <vector>

using namespace std;

// vtx_set: selected label class
// g: selected graph
int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<float> >& g) {
    // selects node from graph given a label, choosing an adjacent node with the maximum degree

    int max_deg = -1;
    int vtx = 0;
    for (int c_vtx : vtx_set) {
        int deg = 0;
        for (float i : g[c_vtx]) {
            if (i != 0) {
                deg += 1;
            }
        }

        if (deg > max_deg) {
            max_deg = deg;
            vtx = c_vtx;
        }
    }
    return vtx;
}
