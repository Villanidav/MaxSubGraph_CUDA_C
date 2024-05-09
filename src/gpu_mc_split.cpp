//
// Created by davide on 5/3/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;

 std::vector<std::vector<double>> g0;
 std::vector<std::vector<double>> g1;
 std::vector<double> edge_labels;
vector<pair<int,int>> m_best;

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

bool matchable( int v, int w, LabelClass lc ) {
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = {};
    v_ring_atoms = lc.get_ring_match_data(vector).at(0);

    if( !v_ring_atoms.empty() ) {
        for(int x : v_ring_atoms){

            if( x == -1 ) {
                return false;
            }
            if(x == w ){
                return true;
            }
        }
        return false;
    }
    return true;
}

vector<LabelClass> genNewLabels(int v, int w, const vector<LabelClass>& lcs) {
    vector<LabelClass> l_draft;
    for(LabelClass label : lcs){

        for(double edge_l : edge_labels){
            std::vector<int> v_conn;
            std::vector<int> w_conn;
            std::vector<std::vector<int> > v_c_rings;

            for(int vtx : hood(v,g0,edge_l)){
                if( std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                    v_conn.push_back(vtx);
                }
            }
            v_c_rings = label.get_ring_match_data(v_conn);

            for(int vtx : hood(w,g1,edge_l)){
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
    return l_draft;
}

//Allocare lo spazio
struct queue_elem{
    vector<LabelClass> labels;
    vector<pair<int,int>> m_local;
};

vector<queue_elem> Q;


bool solve_mcs() {
    queue_elem elem =  Q.back();
    Q.pop_back();
    vector<LabelClass> lcs = elem.labels;
    vector<pair<int,int >> m_local = elem.m_local;
    int bound,v;
    bound = m_local.size() + calc_bound(lcs);
    if ( bound <= m_best.size() ){ if( !Q.empty() ){ return true; } return false;}
    vector<LabelClass*> label_class_pointers;
    label_class_pointers.reserve(lcs.size()+1);
    for (LabelClass& item : lcs) {
        label_class_pointers.push_back(&item);
    }
    vector<LabelClass> temp_lcs;


    for ( LabelClass lc : lcs ) {

        //LabelClass lc = *select_label(label_class_pointers, m_local.size());

        for( int v : lc.g )  {
            v = select_vertex(lc.g, g0);

            for ( int w : lc.h ) {
                if ( !matchable(v,w,lc) ) continue;

                temp_lcs = genNewLabels(v,w,lcs);
                pair<int,int> m_temp;
                m_local.push_back(m_temp);
                queue_elem qel;
                qel.labels = temp_lcs;
                qel.m_local = m_local;
                Q.push_back(qel);
                if ( m_local.size() > m_best.size() ) m_best = m_local;
                m_local.pop_back();

            }
            //lc.remove(0, v);
        }
        //lcs.erase(std::find(lcs.begin(), lcs.end(), lc));
    }
    return true;
}







vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<double>>& g00, const std::vector<std::vector<double>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes) {
    g0 = g00;
    g1 = g11;
    edge_labels = gen_bond_labels(g0, g1);

    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);
    std::vector<LabelClass* > label_class_pointers;
    int v;
    label_class_pointers.reserve(initial_label_classes.size()+1);
    for (LabelClass& item : initial_label_classes) {
        label_class_pointers.push_back(&item);
    }
    vector<pair<int,int>> m_local={{}};
    vector<LabelClass> temp_lcs;

   for( LabelClass lc : initial_label_classes )  {
        //LabelClass lc = *select_label(label_class_pointers, 0 );

        for(int v : lc.g) {
             v = select_vertex(lc.g, g0);

            for ( int w : lc.h ) {
                if ( !matchable(v,w,lc) ) continue;

                temp_lcs = genNewLabels(v,w,initial_label_classes);


                m_local.at(0).first = v;
                m_local.at(0).second = w;
                queue_elem elem;
                elem.labels = temp_lcs; elem.m_local=m_local;
                Q.push_back(elem);
                //manca while dei threads
            }
            //lc.remove(0, v);
        }
        //initial_label_classes.erase(std::find(initial_label_classes.begin(), initial_label_classes.end(), lc));
    }
    //cout<<"\nQsize:"<< Q.size();
    bool flag = true;
    while( flag ) {
        flag = solve_mcs();
    }
    return m_best;
}
