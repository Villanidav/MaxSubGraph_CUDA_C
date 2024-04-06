#include <vector>
#include <algorithm>

class LabelClass {
    private:
        std::vector<int> g;
        std::vector<int> h;
        int adj;
        int label;
        std::vector<int> rings_g;
    
    public:
    // Constructor
    LabelClass(const std::vector<int>& elems_g, const std::vector<int>& elems_h, const std::vector<int>& rings, int adj = 0, int label = 0) 
        : g(elems_g), h(elems_h), adj(adj), label(label), rings_g(rings) {}
    // Getters and Setters
    const std::vector<int>& getG() const { return g; }
    const std::vector<int>& getH() const { return h; }
    int getAdj() const { return adj; }
    int getLabel() const { return label; }
    const std::vector<int>& getRingsG() const { return rings_g; }

    void setG(const std::vector<int>& elems_g) { g = elems_g; }
    void setH(const std::vector<int>& elems_h) { h = elems_h; }
    void setAdj(int adj_value) { adj = adj_value; }
    void setLabel(int label_value) { label = label_value; }
    void setRingsG(const std::vector<int>& rings) { rings_g = rings; }


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

