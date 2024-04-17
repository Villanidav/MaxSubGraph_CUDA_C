#include <vector>
#include <algorithm>
#include "testFra.hpp"
#include "test.hpp"
using namespace std;

void search_mcs(std::vector<std::vector<int> > g0, std::vector<std::vector<int> > g1, std::vector<LabelClass>& label_classes, std::vector<float> edge_labels, std::vector<std::pair<int, int> > m){
    //incumbent has to be GLOBAL
    std::vector<std::pair<int, int> > incumbent;
    int bound = m.size() + calc_bound(label_classes);
    
    if(m.size() > incumbent.size()){
        incumbent = m;
    }
    if(incumbent.size() >= bound){
        return;
    }

    std::vector<LabelClass*> label_class_pointers;
    label_class_pointers.reserve(label_classes.size()); // Reserve space for the pointers

    for (LabelClass& item : label_classes) {
    label_class_pointers.push_back(&item); // Add the address of each element to the new vector
    }

    //da completare
    LabelClass* label_class = select_label(label_class_pointers, m.size());

    // label_class.size() = 0 || label_class = null
    if(label_class == nullptr){
        return;
    }

    int v = select_vertex(label_class->g, g0);
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = label_class -> get_ring_match_data(vector).at(0);

    for( int w : label_class -> h){
        if( (v_ring_atoms.size() != 0) && 
        (std::find(v_ring_atoms.begin(),v_ring_atoms.end(), -1) != v_ring_atoms.end() ) ||
        ( std::find(v_ring_atoms.begin(),v_ring_atoms.end(), w) == v_ring_atoms.end()) ){
            continue;
        }

        std::vector<LabelClass> l_draft;

        for(LabelClass label : label_classes){
            for(float edge_l : edge_labels){
                std::vector<int> v_conn;
                std::vector<int> w_conn;
                std::vector<std::vector<int> > v_c_rings;
                for(int vtx : hood(v,g0,edge_l)){
                    if(std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                        v_conn.push_back(vtx);             
                    }
                }
                v_c_rings = label.get_ring_match_data(v_conn);
                for(int vtx : hood(w,g1,edge_l)){
                    if(std::find(label.h.begin(),label.h.end(),vtx) != label.g.end() ){
                        w_conn.push_back(vtx);             
                    }
                }
                int adj;
                if(v_conn.size() != 0 && w_conn.size() != 0){
                    if(edge_l != 0.0 || label.adj == 1){
                        adj = 1;
                    }else{
                        adj = 0;
                    }
                    LabelClass tmp(v_conn,w_conn,v_c_rings,adj=adj, label.label);
                    l_draft.push_back(tmp);

                }
            }
        }
        std::pair<int, int> nuova_coppia = std::make_pair(v, w);
        m.push_back(nuova_coppia);
        search_mcs(g0, g1, l_draft, edge_labels, m );
    }

    auto it = std::find(label_classes.begin(), label_classes.end(), label_class);
    if (it != label_classes.end()) {
        label_classes.erase(it);
    }
    label_class->remove(0, v);

    
    if(label_class->g.size() > 0){
        LabelClass tmp = *label_class;
        label_classes.push_back(tmp);
    }
    search_mcs(g0, g1, label_classes, edge_labels, m);
}