#include <vector>
#include <algorithm>
#include <iostream>

class LabelClass {


    // Constructor
    public :     
    std::vector<int> g;
    std::vector<int> h;
    int adj;
    std::string label;
    std::vector<std::vector<int> > rings_g;
    
    LabelClass(const std::vector<int> elems_g, const std::vector<int> elems_h, const std::vector<std::vector<int> > rings, int adjj , std::string labell ) 
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
            if(g.size() == 1){
                g.clear();
                rings_g.clear();
            }else{
                auto it = std::find(g.begin(), g.end(), elem);
                int posizione = 0;
                for(int elem : g){
                    posizione++;
                }
                if (it != g.end()) {
                    int idx = std::distance(g.begin(), it);
                    g.erase(g.begin() + idx);
                    rings_g.erase(rings_g.begin() + idx);
                    int posizione = 0;
                    for(int elem : g){
                        posizione++;
                    }
                }else{
                    int posizione = 0;
                }
            }


        } else {
            auto it = std::find(h.begin(), h.end(), elem);
            if (it != h.end()) {
                h.erase(it);
            }
        }
    }

     // Get ring match data method
    std::vector<std::vector<int> > get_ring_match_data( std::vector<int>& elems) {
        std::vector<std::vector<int> > res = {};
        std::vector<int> idxList = {};
        int c;

        if( !elems.empty() ) {

            for (int i : elems) {
                c=0;
                if( !g.empty() ){
                    for ( int k : g ) {
                        if ( k == i ) idxList.push_back(c);
                        c++;
                    }
                }
            }
        }


        if ( !rings_g.empty() && !idxList.empty() ) {
            for ( int j : idxList ) {
                res.push_back(rings_g.at(j));
            }
        }

        return res;
    }

    // Overload the operator==
    bool operator==(const LabelClass& other) const {
        // Compare relevant attributes of LabelClass for equality
        // For example:
        return (this->g == other.g && this->h == other.h && this->label == other.label);
    }


};