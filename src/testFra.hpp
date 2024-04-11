#include <vector>
std::vector<std::vector<int> > gen_rings_classes(
    /*RDKit::ROMol mol0,RDKit::ROMol mol1*/
    const std::vector<std::string> l0,
    const std::vector<std::string> l1,
    const std::vector<std::vector<int> > ring_info_m0,
    const std::vector<std::vector<int> > ring_info_m1
);
//#include <rdkit/RDKit.h>