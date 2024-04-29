#include <vector>
#include <string>
#include "test.hpp" // Include your header file where gen_rotations function is declared

std::vector<std::vector<int>> gen_ring_classes(const RDKit::RWMol& mol0, const RDKit::RWMol& mol1) {
    std::vector<std::string> l0, l1;
    for (const auto& atom : mol0.atoms()) {
        l0.push_back(atom->getSymbol());
    }
    for (const auto& atom : mol1.atoms()) {
        l1.push_back(atom->getSymbol());
    }

    std::vector<std::vector<int>> ring_info_m0, ring_info_m1;
    ring_info_m0 = mol0.getRingInfo()->atomRings();
    cout << "ring info m0 " <<endl;
    for ( std::vector<int> r : ring_info_m0 ) {
        cout << "[ " ;
        for( int a : r )
            cout << a ;
        cout << "[ " <<endl;
    }

    ring_info_m1 = mol1.getRingInfo()->atomRings();
    cout << "ring info m1 " <<endl;
    for ( std::vector<int> r : ring_info_m1 ) {
        cout << "[ " <<endl;
        for( int a : r )
            cout << a ;
        cout << "] " <<endl;
    }

    cout << "DEBAGGGGGGGGG111" <<endl;

    std::vector<std::vector<int> >ring_comp_m0;
    ring_comp_m0.resize(l0.size());
    for ( std::vector<int> r : ring_info_m0 ) {
        cout << "\ndentro " <<endl;
        if( !r.empty() ) {
            for ( int a : r ) {
                if( a < ring_comp_m0.size() ) {
                    cout << "\n secondo for " <<endl;
                    ring_comp_m0.at(a) = {-1};
                }
            }
        }
    }

    cout << "\n RING COMP M0 : \n " << endl;
    cout << "[" ;
    if( !ring_comp_m0.empty() ) {
        for ( std::vector<int> r : ring_comp_m0 ) {
            cout << "[" ;
            if( !r.empty() ) {
                for ( int idx : r)
                {cout << "" <<idx << "" ;}
            }
            cout << "]," ;
        }
    }cout << "]" ;

    cout << "r0 label\n" <<endl;
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