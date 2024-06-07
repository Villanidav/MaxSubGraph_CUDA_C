//
// Created by davide on 5/3/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

using namespace std;
    std::vector<std::vector<float>> g0;
    std::vector<std::vector<float>> g1;
    std::vector<float> edge_labelss;
    vector<pair<int,int>> m_best;


    struct queue_elem{
        vector<LabelClass> labels;
        vector<pair<int,int>> m_local;
    };
    queue<queue_elem> Q;

void printLabelClasss(LabelClass lb) {
    if( true) {
        cout<< lb.label << " [ ";
        cout<< " G("<< lb.g.size() << "): ";
        if(!lb.g.empty()) {for ( int i : lb.g ) cout<<"["<<i<<"]";}
        cout<< " H("<< lb.h.size() << "): ";
        if(!lb.h.empty()) {for ( int i : lb.h ) cout<<"["<<i<<"]";}
        cout<< " RINGS("<< lb.rings_g.size() << "): [";
        for( vector<int> i : lb.rings_g ){cout<<"("<<i.size()<<")"<<"["; for( int j: i ) cout<<j<<", ";  cout<<" ]";}
        cout<<"]";
        cout<< " edge : " <<lb.adj<<" " ;
        cout<< lb.label << " ] "<<endl;
    }
}


LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);


bool matchableee( int v, int w, LabelClass lc ) {
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = {};
    v_ring_atoms = (lc.get_ring_match_data(vector)).at(0);


    if( !v_ring_atoms.empty() ) {
        for(int x : v_ring_atoms){
            //cout<<"  "<<v<< " result =  "<< x;
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

vector<LabelClass> genNewLabelsss(int v, int w, const vector<LabelClass>& lcs) {
    vector<LabelClass> l_draft;
    for(LabelClass label : lcs){

        for(float edge_l : edge_labelss){
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
                LabelClass tmp(v_conn,w_conn,v_c_rings,adj, label.label);
                l_draft.push_back(tmp);
            }
        }
    }
    return l_draft;
}




bool solve_mcs(int size) {

 for( int i = 0 ; i < size ; i++){
        queue_elem elem =  Q.front();
        Q.pop();

    
        vector<LabelClass> lcs = elem.labels;
        vector<pair<int,int >> m_local = elem.m_local;

        std::vector<LabelClass*> label_class_pointers;
        label_class_pointers.reserve(lcs.size());
        for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}

        LabelClass *lcc = select_label(label_class_pointers, m_local.size());
        if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){return true;}

        queue_elem qel;
        LabelClass lc = *lcc;
        pair<int,int> m_temp;

        for( int v : lc.g )  {
            for ( int w : lc.h ) {
                if ( !matchableee(v,w,lc) ) continue;
                m_temp.first = v;
                m_temp.second = w;
                m_local.push_back(m_temp);
                qel.labels = genNewLabelsss(v,w,lcs);
                qel.m_local = m_local;
                Q.push(qel);
                if ( m_local.size() > m_best.size() ) {
                    m_best = m_local;
                    cout<<"\nsize : "<< m_best.size()<<endl;
                    for( pair<int,int> p : m_best ) cout<<"["<< p.first << " , "<< p.second <<" ]";
                }
                m_local.pop_back();
            }
        }
     }
    

   
   return Q.empty();
}


int calcSize(vector<LabelClass> lcs) {
    int result=0;
    for( LabelClass lc : lcs) {
        result = result + (lc.g.size()+lc.h.size());
        for( vector<int> i : lc.rings_g ) {
            result += i.size();
        }
        result += 3;
    }
    return result;
}



vector<pair<int,int>> cpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes){
    g0 = g00;
    g1 = g11;
    edge_labelss = gen_bond_labels(g0, g1);
    int min = std::min(l0.size(), l1.size());
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);
    int size_of_label_classes = calcSize(initial_label_classes);


    vector<pair<int,int>> m_local={{}};
    queue_elem elem;
    int v=0,w=0;
    for( LabelClass lc : initial_label_classes ) {
        v = select_vertex(lc.g,g0);
        w = select_vertex(lc.h,g1);
        //for( int w : lc.h ){
            if ( !matchableee(v,w,lc) ) continue;
            m_local.at(0).first = v;
            m_local.at(0).second = w;
            elem.labels = genNewLabelsss(v,w,initial_label_classes);

            elem.m_local=m_local;
            Q.push(elem);
            //todo
        //}
    }

    int flag;
    flag = Q.size();
    int i = 0;
    do{
   
    solve_mcs(flag);
    flag = Q.size();
    cout<<"\nDIM : "<< flag;
    /*     for( queue_elem e : Q ){
        cout<<"\n\n";
        for( LabelClass l : e.labels){
            printLabelClasss(l);
        }
    }*/
    
    i++;
     }while(flag > 0);
    //for( int i = 0 ; i < 4 ; i++ ) solve_mcs();
    cout << "\n\n\n\n\n\n\n";
    return m_best;
}

