//
// Created by davide on 5/3/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>


 std::vector<std::vector<double>> g0;
 std::vector<std::vector<double>> g1;
 std::vector<double> edge_labels;

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

void gpu_mc_split(const std::vector<std::vector<double>>& g00, const std::vector<std::vector<double>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes) {
    g0 = g00;
    g1 = g11;
    edge_labels = gen_bond_labels(g0, g1);
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);
    std::vector<LabelClass* > label_class_pointers;
    label_class_pointers.reserve(initial_label_classes.size()+1);

for (LabelClass& item : initial_label_classes) {
        label_class_pointers.push_back(&item);
    }
    std::pair<int, int> pair;
    std::vector<std::pair<int, int>> m ;
    std::vector<std::pair<int, int>> m_best ;
    std::queue<LabelClass> labelQueue;

    int n_thread = 0;
    int bound = 0;

    while( !initial_label_classes.empty() ) {
        LabelClass label_class = *select_label( label_class_pointers , m.size() );

        initial_label_classes.erase(std::find(initial_label_classes.begin(), initial_label_classes.end(), label_class));

        bound = m.size() + calc_bound(initial_label_classes);
        if ( bound <= m_best.size() ) continue;
        /*something related to threads*/
        pair = pair_vertex(label_class, g0);
        if( pair.second == -1  ) continue;
        m.push_back(pair);

        std::vector<LabelClass> l_draft;

        if ( m.size() > m_best.size() ) m_best = m;
        for(LabelClass label : initial_label_classes){
            for(double edge_l : edge_labels){
                std::vector<int> v_conn;
                std::vector<int> w_conn;
                std::vector<std::vector<int> > v_c_rings;
                for(int vtx : hood(pair.first,g0,edge_l)){
                    if(std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                        v_conn.push_back(vtx);
                    }
                }
                v_c_rings = label.get_ring_match_data(v_conn);
                for(int vtx : hood(pair.second,g1,edge_l)){
                    if(std::find(label.h.begin(),label.h.end(),vtx) != label.h.end() ){
                        w_conn.push_back(vtx);
                    }
                }
                int adj;
                if(!v_conn.empty() && !w_conn.empty()){
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

        initial_label_classes = l_draft;

    }
    cout<< "M_BEST[";
    for ( std::pair<int,int> p : m_best) {
        cout<< "["<< p.first << ", "<< p.second << "]";
    }
    cout<< "]"<<endl;

}
