//
// Created by francesco on 24/04/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <algorithm>
std::vector<std::pair<int, int>> incumbent;

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

void printLabelClass(LabelClass lb) {
    if( true) {
        cout<< lb.label << " [ ";
        cout<< " G("<< lb.g.size() << "): ";
        if(!lb.g.empty()) {for ( int i : lb.g ) cout<<"["<<i<<"]";}
        cout<< " H("<< lb.h.size() << "): ";
        if(!lb.h.empty()) {for ( int i : lb.h ) cout<<"["<<i<<"]";}
        cout<< " RINGS("<< lb.rings_g.size() << "): [";
        for( vector<int> i : lb.rings_g ){cout<<"["; for( int j: i ) cout<<j<<", ";  cout<<" ]";}
        cout<<"]";
        cout<< " edge : " <<lb.adj<<" " ;
        cout<< lb.label << " ] "<<endl;
    }
}


bool matchablee( int v, int w, LabelClass lc ) {
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


void search_mcs(std::vector<std::vector<double> > g0, std::vector<std::vector<double> > g1,
    std::vector<LabelClass>& label_classes, std::vector<double> edge_labels,
    const std::vector<std::pair<int, int> >& m){


    //incumbent has to be GLOBAL
    int bound;
    bound = m.size() + calc_bound(label_classes);

    if(m.size() > incumbent.size()){
        incumbent = m;
    }
    if(incumbent.size() >= bound){
        return;
    }

    std::vector<LabelClass*> label_class_pointers;
    if( !label_classes.empty() ) label_class_pointers.reserve(label_classes.size()+1); // Reserve space for the pointers


    for (LabelClass& item : label_classes) {
        label_class_pointers.push_back(&item); // Add the address of each element to the new vector
    }


    LabelClass* single_label_class_pointer;


    if ( !m.empty() )  single_label_class_pointer = select_label(label_class_pointers, m.size());
    else single_label_class_pointer = select_label(label_class_pointers, 0);


    // label_class.size() = 0 || label_class = null
    if( m.size()>0 && !single_label_class_pointer  ){
        return;
    }
        LabelClass label_class = *single_label_class_pointer;




        int v = select_vertex(label_class.g, g0);




        std::vector<int> vector;
        vector.push_back(v);



        std::vector<int>  v_ring_atoms = {};
        std::vector<std::vector<int> > returned_data = {};
        returned_data = label_class.get_ring_match_data(vector);

        if ( !returned_data.empty()) {
            v_ring_atoms = label_class.get_ring_match_data(vector).at(0);
        }


        for( int w : label_class.h){
            /*bool flag = false;
            bool flag_m = false;
            if( !v_ring_atoms.empty() ) {
                for(int x : v_ring_atoms){

                    if( x == -1 ) {
                        flag_m = true;
                    }
                    if(x == w ){
                        flag = true;
                        break;
                    }

                }
                if ( !flag || flag_m ) continue;
            }*/
            if ( !matchablee(v,w,label_class) ) continue;



            std::vector<LabelClass> l_draft;

            for(LabelClass label : label_classes){

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

            std::pair<int, int> nuova_coppia = std::make_pair(v, w);

            std::vector<std::pair<int, int> > h = m;

            h.push_back(nuova_coppia);
            cout<<"\n\n new LABEL CLASSES "<< nuova_coppia.first << " , " << nuova_coppia.second<<endl;
            for ( LabelClass lb : l_draft ) printLabelClass(lb);
            search_mcs(g0, g1, l_draft, edge_labels, h );
        }


        label_classes.erase(std::find(label_classes.begin(), label_classes.end(), label_class));
        label_class.remove(0, v);

        if( !label_class.g.empty() ){
            label_classes.push_back(label_class);
        }
        search_mcs(g0, g1, label_classes, edge_labels, m);

}



std::vector<std::pair<int, int>> getIncumbent(){
    return incumbent;
}




std::vector<std::pair<int, int>> mc_split(const std::vector<std::vector<double>> g0, const std::vector<std::vector<double>> g1,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes) {

    incumbent.clear(); // Clear incumbent result
    // Generate label data
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);

    cout<<"MATRIX"<<endl;
    for( std::vector<double> i : g0 ) {
        for( double j : i )
            cout<<"["<<j<<"]";
        cout<<"\n";
    }

    std::vector<double> edge_labels = gen_bond_labels(g0, g1);


    // Search maximum common connected subgraph
    search_mcs(g0, g1, initial_label_classes, edge_labels, {});
    cout<< "INCUMBENT \n[";
    for ( std::pair<int,int> p : incumbent) {
        cout<< "["<< p.first << ", "<< p.second << "]";
    }
    cout<< "]"<<endl;
    return incumbent;
}

