#include <vector>
using namespace std;
// vtx: selected node
// g: selected graph
// edge: bond type
std::vector<int> hood(int vtx, const std::vector<std::vector<float>>& g, float edge) {
    // Return the neighbors of a specified node, with the specified bond type.
    std::vector<int> friends;
    for (std::size_t i = 0; i < g.size(); ++i) {
        if (g[i][vtx] == edge && static_cast<std::size_t>(vtx) != i) {
            friends.push_back(i);
        }
    }

    return friends;
}
