//#include <string>
#include <vector>

using namespace std; // Add this line

std::vector<std::pair<std::string, int>> gen_rotations(const std::string& s);

int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<int>>& g);

std::vector<int> hood(int vtx, const std::vector<std::vector<int>>& g, int edge);