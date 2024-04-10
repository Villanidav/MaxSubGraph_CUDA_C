#include <vector>
#include <algorithm>

class LabelClass {


    // Constructor
    public :     
    std::vector<int> g;
    std::vector<int> h;
    int adj;
    int label;
    std::vector<int> rings_g;
    
    LabelClass(const std::vector<int> elems_g, const std::vector<int> elems_h, const std::vector<int> rings, int adjj , int labell ) 
        {
            g = elems_g;
            h = elems_h;
            rings_g = rings;
            adj = adjj;
            label = labell;
        };

     // Remove method
    void remove(int graph, int elem) {
        if (graph == 0) {
            auto it = std::find(g.begin(), g.end(), elem);
            if (it != g.end()) {
                int idx = std::distance(g.begin(), it);
                g.erase(it);
                rings_g.erase(rings_g.begin() + idx);
            }
        } else {
            auto it = std::find(h.begin(), h.end(), elem);
            if (it != h.end()) {
                h.erase(it);
            }
        }
    }

     // Get ring match data method
    std::vector<int> get_ring_match_data(const std::vector<int>& elems) {
        std::vector<int> res;
        for (int i : elems) {
            auto it = std::find(g.begin(), g.end(), i);
            if (it != g.end()) {
                int idx = std::distance(g.begin(), it);
                res.push_back(rings_g[idx]);
            }
        }
        return res;
    }


};

