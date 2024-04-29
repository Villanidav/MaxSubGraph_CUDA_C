#include <iostream>
#include "test.hpp"

#include <vector>
#include <string>
#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

//questa è la funzione nostra
std::vector<std::vector<int>> gen_ring_classes_1(const RDKit::RWMol& mol0, const RDKit::RWMol& mol1) {
    
    std::vector<std::string> l0, l1;
    for (const auto& atom : mol0.atoms()) {
        l0.push_back(atom->getSymbol());
    }
    for (const auto& atom : mol1.atoms()) {
        l1.push_back(atom->getSymbol());
    }

    std::vector<std::vector<int>> ring_info_m0, ring_info_m1;
    ring_info_m0 = mol0.getRingInfo()->atomRings();

   
    ring_info_m1 = mol1.getRingInfo()->atomRings();
  

    std::vector<std::vector<int> >ring_comp_m0; // Initialize with -1
    ring_comp_m0.resize(l0.size());
    for ( std::vector<int> r : ring_info_m0 ) {
   
        if( !r.empty() ) {
            for ( int a : r ) {
                if( a < ring_comp_m0.size() ) {
                    ring_comp_m0.at(a) = {-1};
                }
            }
        }
    }

    
    for (const std::vector<int>& r0 : ring_info_m0) {
        std::string r0_label;
        if( !r0.empty() ) {
            for (int atomIdx : r0) {
                r0_label += l0[atomIdx];
            }
        }
     
        std::string r0_label_rev = r0_label;
        std::reverse(r0_label_rev.begin(), r0_label_rev.end());

        for (const std::vector<int>& r1 : ring_info_m1) {
            if( !r1.empty() ) {
                if (r0.size() == r1.size()) {
                    std::string r1_label;
                    for (int atomIdx : r1) {
                        r1_label += l1[atomIdx];
                    }

                    std::vector<std::pair<std::string, int>> rotations = gen_rotations(r0_label);
                    std::vector<std::pair<std::string, int>> inv_rotations = gen_rotations(r0_label_rev);

                    std::vector<std::pair<std::string, int>> r0_rots;
                    std::vector<std::pair<std::string, int>> r0_rev_rots;

                    for (const auto& rot : rotations) {
                        if (r1_label == rot.first) {
                            r0_rots.push_back(rot);
                        }
                    }
                    for (const auto& rot : inv_rotations) {
                        if (r1_label == rot.first) {
                            r0_rev_rots.push_back(rot);
                        }
                    }

                    for(std::pair<std::string, int> comb : r0_rev_rots ){
                        cout << "\n" << comb.first << " " << comb.second << "\n";
                    }
                    for(std::pair<std::string, int> comb : r0_rots ){
                        cout << "\n" << comb.first << " " << comb.second << "\n";
                    }
                    
                    int i = 0;
                    for (const auto& rot : r0_rots) {
                        for (size_t idx = 0; idx < r0.size() ; ++idx) {
                            int targetIdx = idx - rot.second;
                            if(targetIdx < 0){
                                targetIdx = r1.size() + targetIdx;
                            }
                            cout << "\n"<<r0[idx] << "       ";
                            if (ring_comp_m0[r0[idx]][0] == -1) {
                                
                                ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                                cout <<"\n" << i << "direc 0: la posizione idx dell'anello r1 considerata: "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            } else {
                                ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                                cout <<"\n" << i << "direc 1: la posizione idx dell'anello r1 considerata: "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";;
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            }
                        }
                    }

                    for (const auto& rot : r0_rev_rots) {
                        for (size_t idx = 0; idx < r0.size() ; ++idx) {
                            cout << "\n"<< r0[idx] << "       ";
                            int targetIdx = r1.size() - idx - 1 - rot.second;
                            if(targetIdx < 0){
                                targetIdx = r1.size() + targetIdx;
                            }
                            if (ring_comp_m0[r0[idx]][0] == -1) {
                                ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                                cout <<"\n" << i << "rev 0: la posizione idx dell'anello r1 considerata:  "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";;
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            } else if(find(ring_comp_m0.at(r0[idx]).begin(), ring_comp_m0.at(r0[idx]).end(), r1[targetIdx] ) == ring_comp_m0.at(r0[idx]).end() ) {
                                ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                                cout <<"\n" << i << " rev 1 :la posizione idx dell'anello r1 considerata:  "<<targetIdx << "ed è : "<< r1[targetIdx] <<"\n";
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            }
                        }
                    }
                }
            }

        }
    }

    return ring_comp_m0;
}

std::vector<std::vector<int>> s_smiles_mcs( std::string& smile0,  std::string& smile1 ) {

    RWMol mol0 = *SmilesToMol(smile0);
    RWMol mol1 = *SmilesToMol(smile1);

    return gen_ring_classes_1(mol0,mol1);
    
}



//funzione di prova
std::vector<std::vector<int>> prova_gen_ring_classes(std::vector<std::string> l0, std::vector<std::string>  l1, const std::vector<std::vector<int> > ring_info_m0, const std::vector<std::vector<int> > ring_info_m1 ) {
    
    std::vector<std::vector<int> >ring_comp_m0; // Initialize with -1
    ring_comp_m0.resize(l0.size());
    for ( std::vector<int> r : ring_info_m0 ) {
   
        if( !r.empty() ) {
            for ( int a : r ) {
                if( a < ring_comp_m0.size() ) {
                    ring_comp_m0.at(a) = {-1};
                }
            }
        }
    }

    
    for (const std::vector<int>& r0 : ring_info_m0) {
        std::string r0_label;
        if( !r0.empty() ) {
            for (int atomIdx : r0) {
                r0_label += l0[atomIdx];
            }
        }
     
        std::string r0_label_rev = r0_label;
        std::reverse(r0_label_rev.begin(), r0_label_rev.end());

        for (const std::vector<int>& r1 : ring_info_m1) {
            if( !r1.empty() ) {
                if (r0.size() == r1.size()) {
                    std::string r1_label;
                    for (int atomIdx : r1) {
                        r1_label += l1[atomIdx];
                    }

                    std::vector<std::pair<std::string, int>> rotations = gen_rotations(r0_label);
                    std::vector<std::pair<std::string, int>> inv_rotations = gen_rotations(r0_label_rev);

                    std::vector<std::pair<std::string, int>> r0_rots;
                    std::vector<std::pair<std::string, int>> r0_rev_rots;

                    for (const auto& rot : rotations) {
                        if (r1_label == rot.first) {
                            r0_rots.push_back(rot);
                        }
                    }
                    for (const auto& rot : inv_rotations) {
                        if (r1_label == rot.first) {
                            r0_rev_rots.push_back(rot);
                        }
                    }

                    for(std::pair<std::string, int> comb : r0_rev_rots ){
                        cout << "\n" << comb.first << " " << comb.second << "\n";
                    }
                    for(std::pair<std::string, int> comb : r0_rots ){
                        cout << "\n" << comb.first << " " << comb.second << "\n";
                    }
                    
                    int i = 0;
                    for (const auto& rot : r0_rots) {
                        for (size_t idx = 0; idx < r0.size() ; ++idx) {
                            int targetIdx = idx - rot.second;
                            if(targetIdx < 0){
                                targetIdx = r1.size() + targetIdx;
                            }
                            cout << "\n"<<r0[idx] << "       ";
                            if (ring_comp_m0[r0[idx]][0] == -1) {
                                
                                ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                                cout <<"\n" << i << "direc 0: la posizione idx dell'anello r1 considerata: "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            } else {
                                ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                                cout <<"\n" << i << "direc 1: la posizione idx dell'anello r1 considerata: "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";;
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            }
                        }
                    }

                    for (const auto& rot : r0_rev_rots) {
                        for (size_t idx = 0; idx < r0.size() ; ++idx) {
                            cout << "\n"<< r0[idx] << "       ";
                            int targetIdx = r1.size() - idx - 1 - rot.second;
                            if(targetIdx < 0){
                                targetIdx = r1.size() + targetIdx;
                            }
                            if (ring_comp_m0[r0[idx]][0] == -1) {
                                ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                                cout <<"\n" << i << "rev 0: la posizione idx dell'anello r1 considerata:  "<<targetIdx  << "ed è : "<< r1[targetIdx] <<"\n";;
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            } else if(find(ring_comp_m0.at(r0[idx]).begin(), ring_comp_m0.at(r0[idx]).end(), r1[targetIdx] ) == ring_comp_m0.at(r0[idx]).end() ) {
                                ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                                cout <<"\n" << i << " rev 1 :la posizione idx dell'anello r1 considerata:  "<<targetIdx << "ed è : "<< r1[targetIdx] <<"\n";
                                for(std::vector<int> vect : ring_comp_m0 ){
                                    cout << "[" ;
                                    for(int pos : vect){
                                        cout << pos << " ";
                                    }
                                    cout << "]" ;
                                }
                                i++;
                            }
                        }
                    }
                }
            }

        }
    }

    return ring_comp_m0;
}



int main() {
/* TESTING CREAZIONE LABEL
    std::vector<int> g;
    std::vector<int> h;
    int adj;
    int label;
    std::vector<std::vector<int> > rings_g;
    rings_g.resize(3);

    g.push_back(8);
    g.push_back(10);
    g.push_back(2);
    
    h.push_back(45);
    h.push_back(40);
    h.push_back(25);
    
    rings_g.at(0).push_back(45);
    rings_g.at(0).push_back(40);

    rings_g.at(1).push_back(25);
    

    LabelClass label(g,h,rings_g,0,0);
    
    std::cout << "\nStampo gli elementi di g:\n ";
    for (int i : label.g)
    {
        std::cout << i << " ";
    }

    std::cout << "\nStampo i match tra l'elemento con l'indice 8 all'interno della molecola ( il primo all'interno di g) ";
    std::vector<int> prova;
    prova.push_back(8);
    prova.push_back(10);
    std::vector<std::vector<int> > rings_ritornati = label.get_ring_match_data(prova);

    
    for (int i = 0; i < rings_ritornati.size(); i++)
    {
        std::cout <<  " \n";
        for(int j: rings_ritornati.at(i)){
            std::cout << j << " ";
        }
    }
    std::cout << " \n";

    label.remove(0,10);

    rings_ritornati = label.get_ring_match_data(prova);

    
    for (int i = 0; i < rings_ritornati.size(); i++)
    {
        std::cout << i << " \n";
        for(int j: rings_ritornati.at(i)){
            std::cout << j << " ";
        }
    }
    
*/
/*TESTING GEN RING CLASSES ANELLO ARTIFICIALE
    const std::vector<std::string> l0({"C","C","O","C","H","H","O"});
    const std::vector<std::string> l1({"C","H","F","O","H","C","B","R"});
    const std::vector<std::vector<int> > ring_info_m0({{0,1,2},{4,3,5}});
    const std::vector<std::vector<int> > ring_info_m1({{4,1,0},{1,2,3},{1,4,5},{5,0,3}});
    std::vector<std::vector<int> > genereted_rings = prova_gen_ring_classes(l0,l1,ring_info_m0,ring_info_m1);
    int pos = 0;
    for(std::vector<int> vect : genereted_rings ){
        cout << "\n: posizione: " << pos << " \n";
        pos++;
        for(int i : vect){
            cout << i << " ";     
        };
    };


*/

/*GENERATE RING CLASSES CON SMIILES*/

string smile0 = "CN(c1ccc(cc1)c2nnn(CC(=O)Nc3ccc4nc(oc4c3)c5ccccc5Cl)n2)c6cc(C)c(N)cn6";
string smile1 = "O=C(Cn1nnc(n1)c2ccc(Nc3ccccn3)cc2)Nc4ccc5nc(oc5c4)c6ccccc6";


std::vector<std::vector<int> > genereted_rings = s_smiles_mcs(smile0,smile1);
    int pos = 0;
    for(std::vector<int> vect : genereted_rings ){
        cout << "\n: posizione: " << pos << " \n";
        pos++;
        for(int i : vect){
            cout << i << " ";     
        };
    };

}