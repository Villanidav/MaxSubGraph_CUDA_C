//
// Created by davide on 4/18/24.
// Modified by Francesco on 24/04
//
#include "test.hpp"

#include <iostream>
#include <vector>
#include <algorithm>
#include <GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RingInfo.h>
#include <GraphMol/Bond.h>

using namespace RDKit;

vector<vector<double>> getAdjacencyMatrix(const RWMol& mol) {
    int numAtoms = mol.getNumAtoms();
    vector<vector<double>> adjacencyMatrix(numAtoms, vector<double>(numAtoms, 0.0));

    // Iterare su tutti i legami nella molecola e aggiornare la matrice di adiacenza con il peso dei legami
    for (RWMol::BondIterator bondIt = mol.beginBonds(); bondIt != mol.endBonds(); ++bondIt) {
        const Bond *bond = *bondIt;
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

        // Aggiornare la matrice di adiacenza
        adjacencyMatrix[beginAtomIdx][endAtomIdx] = bondWeight;
        adjacencyMatrix[endAtomIdx][beginAtomIdx] = bondWeight; // Poiché la matrice di adiacenza è simmetrica per i grafi non orientati
    }

    return adjacencyMatrix;
}



//ho cambiato la return di questa funzione per fini di test
void mol_mcs(const RDKit::RWMol &mol0, const RDKit::RWMol &mol1, int bond_match=1, int ring_match=1, int return_map=0) {
    std::vector<RDKit::RWMol> mols = {mol0, mol1};

    std::vector<std::string> l0, l1;
    for (const auto &atom : mol0.atoms()) {
        l0.push_back(atom->getSymbol());
    }
    for (const auto &atom : mol1.atoms()) {
        l1.push_back(atom->getSymbol());
    }

    std::vector<std::vector<std::string>> label_ring_data = {l0, l1};


    std::vector<std::vector<double>> go,g1;
    if (bond_match) {
        g0 = getAdjacencyMatrix(mol0);
        g1 = getAdjacencyMatrix(mol1);
    } else {
        g0 = getAdjacencyMatrix(mol0);
        g1 = getAdjacencyMatrix(mol1);
    }

    if (ring_match) {
        std::vector<std::vector<int>> ring_info = {mol0.getRingInfo()->atomRings(), mol1.getRingInfo()->atomRings()};
        for (size_t mol_idx = 0; mol_idx < 2; ++mol_idx) {
            for (const auto &ring : ring_info[mol_idx]) {
                for (int atom_idx : ring) {
                    if (label_ring_data[mol_idx][atom_idx].back() != 'R') {
                        label_ring_data[mol_idx][atom_idx] += 'R';
                    }
                }
            }
        }
    }

    std::vector<std::vector<int>> ring_classes = gen_ring_classes(mol0, mol1);
    //ho cambiato la seguente riga dal codice Python
    std::vector<std::pair<int, int>> mapping = mc_split(g0, g1, l0, l1, ring_classes);

    /*
    std::sort(mapping.begin(), mapping.end());

    std::vector<int> mapped_atom_idxs_g0;
    for (const auto &pair : mapping) {
        mapped_atom_idxs_g0.push_back(pair.first);
    }

    std::vector<std::string> mcs_labels;
    for (int idx_g0 : mapped_atom_idxs_g0) {
        mcs_labels.push_back(l0[idx_g0]);
    }

    RDKit::INT_MATRIX mcs_matrix = g0;
    for (int idx = g0.size() - 1; idx >= 0; --idx) {
        if (std::find(mapped_atom_idxs_g0.begin(), mapped_atom_idxs_g0.end(), idx) == mapped_atom_idxs_g0.end()) {
            mcs_matrix.erase(mcs_matrix.begin() + idx);
            for (auto &row : mcs_matrix) {
                row.erase(row.begin() + idx);
            }
        }
    }

    RDKit::RWMol mcs = g2mol(mcs_labels, mcs_matrix);
    if (return_map) {
        // return incumbent
    }

    return mcs;
     */


    return
}

