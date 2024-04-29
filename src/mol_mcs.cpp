//
// Created by davide on 4/18/24.
// Modified by Francesco on 24/04
//
#include "test.hpp"

#include <iostream>
#include <vector>
#include <algorithm>



vector<vector<double>> getAdjacencyMatrix(const RWMol& mol) {
    int numAtoms = mol.getNumAtoms();
    vector adjacencyMatrix(numAtoms, vector<double>(numAtoms, 0.0));

    // Iterare su tutti i legami nella molecola e aggiornare la matrice di adiacenza con il peso dei legami
    // Iterare su tutti i legami nella molecola e aggiornare la matrice di adiacenza con il peso dei legami
    for (const Bond *bond : mol.bonds()) {
        const Atom *beginAtom = bond->getBeginAtom();
        const Atom *endAtom = bond->getEndAtom();
        Bond::BondType bondType = bond->getBondType();

        // Ottenere gli indici degli atomi connessi
        int beginAtomIdx = beginAtom->getIdx();
        int endAtomIdx = endAtom->getIdx();

        // Impostare il peso del legame in base al tipo di legame
        double bondWeight = 1.0; // Peso predefinito per il legame singolo
        if (bondType == Bond::DOUBLE) {
            bondWeight = 2.0;
        } else if (bondType == Bond::TRIPLE) {
            bondWeight = 3.0;
        } else if (bondType == Bond::AROMATIC) {
            bondWeight = 1.5; // Peso per legame aromatico
        }
        else if (bondType == Bond::QUADRUPLE) {
            bondWeight = 4.0; // Peso per legame aromatico
        }
        else if (bondType == Bond::QUINTUPLE) {
            bondWeight = 5.0; // Peso per legame aromatico
        }
        else if (bondType == Bond::HEXTUPLE) {
            bondWeight = 6.0; // Peso per legame aromatico
        }

        // Aggiornare la matrice di adiacenza
        adjacencyMatrix[beginAtomIdx][endAtomIdx] = bondWeight;
        adjacencyMatrix[endAtomIdx][beginAtomIdx] = bondWeight; // Poiché la matrice di adiacenza è simmetrica per i grafi non orientati
    }

    return adjacencyMatrix;
}



//ho cambiato la return di questa funzione per fini di test
ROMol mol_mcs(const RDKit::RWMol &mol0, const RDKit::RWMol &mol1, int bond_match, int ring_match, int return_map) {
    std::vector mols = {mol0, mol1};

    std::vector<std::string> l0, l1;
    for (const auto &atom : mol0.atoms()) {
        l0.push_back(atom->getSymbol());
    }
    for (const auto &atom : mol1.atoms()) {
        l1.push_back(atom->getSymbol());
    }

    std::pair<std::vector<std::string>, std::vector<std::string>> label_ring_data = {l0, l1};

    //stampa per prove
    cout << "\n first molecule: " <<endl;

    for ( string idx : l0 )
        cout << idx;

    cout << "\n second molecule: " <<endl;
    for ( string idx : l1 )
        cout  << idx;

    std::vector<std::vector<double>> g0,g1;

    if (bond_match) {
        g0 = getAdjacencyMatrix(mol0);
        g1 = getAdjacencyMatrix(mol1);
    } else {
        g0 = getAdjacencyMatrix(mol0);
        g1 = getAdjacencyMatrix(mol1);
    }

    //print matrix
    std::cout<<"\n first adiacent Matrix";
    for (vector<double> line : g0){
        std::cout<<"\n";
        for(double i : line){
            std::cout<<" "<<i;
        }
    }

    std::cout<<"\n second adiacent Matrix";
    for (vector<double> line : g1){
        std::cout<<"\n";
        for(double i : line){
            std::cout<<" "<<i;
        }
    }

    cout << " INIZIO COSA ANELLI  " << endl;

    if (ring_match) {
        std::pair< std::vector<std::vector<int> > , std::vector<std::vector<int>> > ring_info ;

        for ( const std::vector<int>& var : mol0.getRingInfo()->atomRings() )
            ring_info.first.push_back(var);

        for ( const std::vector<int>& var : mol1.getRingInfo()->atomRings() )
            ring_info.second.push_back(var);


        //cout << "  SIZE : " << ring_info.size() <<endl;
        if( !label_ring_data.first.at(0).empty() ) {
            for (auto & ring : ring_info.first) {
                for ( int atm_idx : ring ) {
                    if (  // Check if molecule index is valid
                        !label_ring_data.first.empty() && // Check if sub-vector is not empty
                        atm_idx < label_ring_data.first.size()) { // Check if atom index is valid

                        if (!label_ring_data.first[atm_idx].empty() && label_ring_data.first[atm_idx].back() != 'R') {
                            label_ring_data.first[atm_idx].push_back('R');
                        }
                        }
                }
            }
        }
        if( !label_ring_data.second.at(0).empty() ) {
            for (auto & ring : ring_info.second) {
                for ( int atm_idx : ring ) {
                    if (  // Check if molecule index is valid
                        !label_ring_data.second.empty() && // Check if sub-vector is not empty
                        atm_idx < label_ring_data.second.size()) { // Check if atom index is valid

                        if (!label_ring_data.second[atm_idx].empty() && label_ring_data.second[atm_idx].back() != 'R') {
                            label_ring_data.second[atm_idx].push_back('R');
                        }
                        }
                }
            }
        }
    }

    cout << " LABEL RING DATA \n";
    for ( string r : label_ring_data.first)
        cout<< " "<< r;
    cout << "\n 2 \n";
    for ( string r : label_ring_data.second)
        cout<< " "<< r;


    cout << "\n FINE COSA ANELLI  " << endl;

    //todo risolvere questione stringa inversa

    std::vector<std::vector<int>> ring_classes = gen_ring_classes(mol0, mol1);

    cout << "\n ring classes : \n " << endl;
    cout << "[" ;
    if( !ring_classes.empty() ){
        for ( std::vector<int> r : ring_classes ) {
            cout << "[" ;
            if( !r.empty() ) {
                for ( int idx : r)
                {cout << "'" <<idx << "'" ;}
                cout << "]," ;
            }
        }

    }    cout << "]" ;
    cout << "\n step up classes : \n " << endl;
    //ho cambiato la seguente riga dal codice Python
    std::vector<std::pair<int, int>> mapping = mc_split(g0, g1, label_ring_data.first, label_ring_data.second, ring_classes);

    cout << "\n step ritorno classes : \n " << endl;
    std::sort(mapping.begin(), mapping.end());

    std::vector<int> mapped_atom_idxs_g0;
    for (const auto &pair : mapping) {
        mapped_atom_idxs_g0.push_back(pair.first);
    }

    std::vector<std::string> mcs_labels;
    for (int idx_g0 : mapped_atom_idxs_g0) {
        mcs_labels.push_back(l0[idx_g0]);
    }

    vector<vector<double>> mcs_matrix = g0;
    for (int idx = g0.size() - 1; idx >= 0; --idx) {
        if (std::find(mapped_atom_idxs_g0.begin(), mapped_atom_idxs_g0.end(), idx) == mapped_atom_idxs_g0.end()) {
            mcs_matrix.erase(mcs_matrix.begin() + idx);
            for (auto &row : mcs_matrix) {
                row.erase(row.begin() + idx);
            }
        }
    }

    ROMol mcs = g2mol(mcs_labels, mcs_matrix);

    /*if (return_map) {
         //return incumbent;
    }*/

    return mcs;


}