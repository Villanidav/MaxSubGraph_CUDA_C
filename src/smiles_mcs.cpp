#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using namespace RDKit;

int smiles_mcs( std::string& s0,  std::string& s1, int bond_match = 1, int ring_match = 1) {
    /*RWMol mol0 = RDKit::SmilesToMol(s0);
    RWMol mol1 = RDKit::SmilesToMol(s1);
    return 1;*/
    std::unique_ptr<RWMol> mol0(SmilesToMol(s0));
    std::unique_ptr<RWMol> mol1(SmilesToMol(s1));

    return 1;
}
//['/home/davide/miniconda3/envs/my-rdkit-env/lib/python3.12/site-packages/rdkit']
