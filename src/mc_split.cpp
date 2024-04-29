//
// Created by francesco on 24/04/24.
//

#include "test.hpp"
#include <vector>
#include <string>
#include <algorithm>
std::vector<std::pair<int, int>> incumbent;

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

void search_mcs(std::vector<std::vector<double> > g0, std::vector<std::vector<double> > g1, std::vector<LabelClass>& label_classes, std::vector<double> edge_labels, std::vector<std::pair<int, int> > m){

    cout << "\n INIZIO RICORSIONE:   "<<endl ;
    cout << "\n m.size() = " << m.size() << "  incumbent.size() = " << incumbent.size() << "  bound = " << m.size() + calc_bound(label_classes) << endl;
    //incumbent has to be GLOBAL
    int bound;
    if( !m.empty() )  bound = m.size() + calc_bound(label_classes);
    else bound = 0 + calc_bound(label_classes);

    if(m.size() > incumbent.size()){
        incumbent = m;
    }
    if(incumbent.size() >= bound){
       // cout << "\n fine \n   "<<endl ;
        return;
    }
    std::vector<LabelClass*> label_class_pointers;
    if( !label_classes.empty() ) label_class_pointers.reserve(label_classes.size()); // Reserve space for the pointers


    for (LabelClass& item : label_classes) {
        label_class_pointers.push_back(&item); // Add the address of each element to the new vector
    }

    LabelClass* label_class_pointer;
   

    //da completare
    if ( !m.empty() )  label_class_pointer = select_label(label_class_pointers, m.size());
    else label_class_pointer = select_label(label_class_pointers, 0);



    // label_class.size() = 0 || label_class = null
    if( label_class_pointer == nullptr ){
        return;
    }else{
        LabelClass label_class = *label_class_pointer;


        

    //    std::cout << "\n label_class selected is: "<<label_class.label<<endl;
        

        int v = select_vertex(label_class.g, g0);
        
        
    //    std::cout << "\nv selezionato: "<<v<<endl;
        
        std::vector<int> vector;
        vector.push_back(v);



        std::vector<int>  v_ring_atoms = {};
        std::vector<std::vector<int> > returned_data = {};
        returned_data = label_class.get_ring_match_data(vector);

        if ( !returned_data.empty()) {
            v_ring_atoms = label_class.get_ring_match_data(vector).at(0);
        }

        
        for( int w : label_class.h){

            if(   !v_ring_atoms.empty()    &&
                ( std::find(v_ring_atoms.begin(), v_ring_atoms.end(), -1) == v_ring_atoms.end() ||
                std::find(v_ring_atoms.begin(),v_ring_atoms.end(), w) == v_ring_atoms.end() ) ){
                //cout << "continue : "<<endl ;
                continue;
            }

            std::vector<LabelClass> l_draft;

            for(LabelClass label : label_classes){
                for(double edge_l : edge_labels){
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
            
            for(std::pair<int, int> pair : h){
                cout << "\n[ " << pair.first << " ," << pair.second << "]" <<endl;
            } 
            search_mcs(g0, g1, l_draft, edge_labels, h );
        }
       // cout << "\n salto, fine primo w:   "<<endl ;


        auto it = std::find(label_classes.begin(), label_classes.end(), label_class);
        if (it != label_classes.end()) {
            int index = std::distance(label_classes.begin(), it);
            label_classes.erase(it);
        }
        label_class.remove(0, v);

        if( !label_class.g.empty() ){
            label_classes.push_back(label_class);
        }
        search_mcs(g0, g1, label_classes, edge_labels, m);
    }   
}



std::vector<std::pair<int, int>> getIncumbent(){
    return incumbent;
}




std::vector<std::pair<int, int>> mc_split(const std::vector<std::vector<double>> g0, const std::vector<std::vector<double>> g1,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes) {

    incumbent.clear(); // Clear incumbent result
    cout << "\n step up classes 2: \n " << endl;
    // Generate label data
    std::vector<LabelClass> initial_label_classes = gen_initial_labels(l0, l1, ring_classes);
/*
    for( LabelClass l : initial_label_classes ) {
        cout<< "\nlabel : " << l.label <<endl;
        cout<< "val in g" << endl;
        for( int g : l.g)
            cout<< " " << g << " ";
        cout<< "\nval in h" << endl;
        for( int g : l.h)
            cout<< " " << g << " ";
        cout<< "\nrings.g [" ;
        for ( vector<int> r : l.rings_g) {
            cout<< "["  ;
            for ( int a : r)
                cout<<" " <<a << " ";
            cout<< "]"  ;
        }

        cout<< "]" << endl;
    }

   cout << "\n fino a qui porta tutto \nfunzione mc_split riga 167";

*/
    std::ofstream outfile("mc_split_output.txt");
    std::streambuf* coutbuf = std::cout.rdbuf();
    std::cout.rdbuf(outfile.rdbuf());


    std::vector<double> edge_labels = gen_bond_labels(g0, g1);
    cout << "\nBOND LABELS: "  ;
    for ( double d : edge_labels)
        cout << " " <<d  ;

    // Search maximum common connected subgraph
    search_mcs(g0, g1, initial_label_classes, edge_labels, {});
    std::cout.rdbuf(coutbuf);
    return incumbent;
}