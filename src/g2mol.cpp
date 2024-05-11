//
// Created by davide on 4/24/24.
//
#include "test.hpp"
#include <RDGeneral/test.h>

// Function to convert a graph (adjacency matrix and labels) to an RDKit molecule
RDKit::ROMol g2mol( std::vector<std::string>& labels,  std::vector<std::vector<float>>& adj) {
    // Create an empty editable molecule
    RWMol mol;

    // Map to store node index in adj matrix to molecule atom index
    std::unordered_map<int, int> node_to_idx;

    // Add atoms to the molecule
    for (int i = 0; i < labels.size(); ++i) {
        // Create the atom object (ownership remains here)
        Atom *atom(new Atom(labels[i]));

        // Add the atom to the molecule and store its index
        int mol_idx = mol.addAtom(atom);
        node_to_idx[i] = mol_idx;
    }

    // Add bonds between adjacent atoms
    for (size_t ix = 0; ix < adj.size(); ++ix) {
        for (size_t iy = ix + 1; iy < adj[ix].size(); ++iy) {  // Only iterate upper triangle

            float bond_type = adj[ix][iy];
            if (bond_type == 0.0) {
                continue;
            }

            RDKit::Bond::BondType rdkit_bond_type;
            if (bond_type == 1.5) {
                rdkit_bond_type = RDKit::Bond::AROMATIC;
            } else if (bond_type == 1.0) {
                rdkit_bond_type = RDKit::Bond::SINGLE;
            } else if (bond_type == 2.0) {
                rdkit_bond_type = RDKit::Bond::DOUBLE;
            } else if (bond_type == 3.0) {
                rdkit_bond_type = RDKit::Bond::TRIPLE;
            } else if (bond_type == 4.0) {
                rdkit_bond_type = RDKit::Bond::QUADRUPLE;
            } else if (bond_type == 5.0) {
                rdkit_bond_type = RDKit::Bond::QUINTUPLE;
            } else if (bond_type == 6.0) {
                rdkit_bond_type = RDKit::Bond::HEXTUPLE;
            } else {
                // Handle unsupported bond types (throw exception or log warning)
                throw std::runtime_error("Unsupported bond type: " + std::to_string(bond_type));
            }

            mol.addBond(node_to_idx[ix], node_to_idx[iy], rdkit_bond_type);
        }
    }

    // Convert RWMol to a final, non-editable molecule object
    return ROMol(mol);

}
