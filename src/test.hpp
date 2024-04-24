//#include <string>
#include <vector>
#include "Label.cpp"

#include <rdkit/GraphMol/GraphMol.h>
#include <GraphMol/FileParsers/MolSupplier.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <Numerics/Matrix.h>


using namespace std; // Add this line
using namespace RDKit; // Add this line

std::vector<std::vector<int> > gen_rings_classes(
    /*RDKit::ROMol mol0,RDKit::ROMol mol1*/
    const std::vector<std::string> l0,
    const std::vector<std::string> l1,
    const std::vector<std::vector<int> > ring_info_m0,
    const std::vector<std::vector<int> > ring_info_m1
);

std::vector<std::pair<std::string, int> > gen_rotations(const std::string& s);

int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<double> >& g);


std::vector<int> hood(int vtx, const std::vector<std::vector<double> >& g, double edge);



int calc_bound( std::vector<LabelClass>& label_classes);

std::vector<double> gen_bond_labels(const std::vector<std::vector<double> >& g0, const std::vector<std::vector<double> >& g1);

ROMol smiles_mcs( std::string& smile0,  std::string& smile1, int bond_match = 1, int ring_match = 1);

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string>& l0, const std::vector<std::string>& l1,     std::vector<std::vector<int> >& ring_classes);

std::vector<std::vector<int>> gen_ring_classes(const RDKit::RWMol& mol0, const RDKit::RWMol& mol1);


std::vector<std::pair<int, int>> mc_split(const std::vector<std::vector<double>> g0, const std::vector<std::vector<double>> g1,
                                          const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                          std::vector<std::vector<int> >& ring_classes);

void search_mcs(std::vector<std::vector<double> > g0, std::vector<std::vector<double> > g1, std::vector<LabelClass>& label_classes, std::vector<double> edge_labels, std::vector<std::pair<int, int> > m);

RDKit::ROMol g2mol( std::vector<std::string>& labels,  std::vector<std::vector<double>>& adj) ;

ROMol mol_mcs(const RDKit::RWMol &mol0, const RDKit::RWMol &mol1, int bond_match=1, int ring_match=1, int return_map=0) ;