#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/test.h>
using namespace RDKit;

RWMol smiles_mcs( std::string& smile0,  std::string& smile1, int bond_match = 1, int ring_match = 1) {

    RWMol mol0 = *SmilesToMol(smile0);
    RWMol mol1 = *SmilesToMol(smile1);

    return mol0;
}
//['/home/davide/miniconda3/envs/my-rdkit-env/lib/python3.12/site-packages/rdkit']