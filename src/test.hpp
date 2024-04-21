//#include <string>
#include <vector>
#include "Label.cpp"

#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>


using namespace std; // Add this line

std::vector<std::vector<int> > gen_rings_classes(
    /*RDKit::ROMol mol0,RDKit::ROMol mol1*/
    const std::vector<std::string> l0,
    const std::vector<std::string> l1,
    const std::vector<std::vector<int> > ring_info_m0,
    const std::vector<std::vector<int> > ring_info_m1
);

std::vector<std::pair<std::string, int> > gen_rotations(const std::string& s);

int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<int> >& g);

std::vector<int> hood(int vtx, const std::vector<std::vector<int> >& g, int edge);

LabelClass* select_label(std::vector<LabelClass*>& label_classes, int map_size);

int calc_bound(const std::vector<LabelClass>& label_classes);

std::vector<float> gen_bond_labels(const std::vector<std::vector<float > >& g0, const std::vector<std::vector<float> >& g1);

void smiles_mcs( std::string& s0,  std::string& s1, int bond_match = 1, int ring_match = 1);

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string> l0, const std::vector<std::string> l1,     std::vector<std::vector<int> > ring_classes);