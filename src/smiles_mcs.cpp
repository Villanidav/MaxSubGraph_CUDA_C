#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/test.h>
#include "test.hpp"
using namespace RDKit;


ROMol smiles_mcs( std::string& smile0,  std::string& smile1, int bond_match , int ring_match ) {

    RWMol mol0 = *SmilesToMol(smile0);
    RWMol mol1 = *SmilesToMol(smile1);

    return mol_mcs(mol0, mol1, 1,1,0);
}
//['/home/davide/miniconda3/envs/my-rdkit-env/lib/python3.12/site-packages/rdkit']