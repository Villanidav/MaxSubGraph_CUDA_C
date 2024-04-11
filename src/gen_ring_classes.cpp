#include <vector>
#include "test.hpp"

using namespace std;

std::vector<std::vector<int> > gen_rings_classes(
    /*RDKit::ROMol mol0,RDKit::ROMol mol1*/
    const std::vector<std::string> l0,
    const std::vector<std::string> l1,
    const std::vector<std::vector<int> > ring_info_m0,
    const std::vector<std::vector<int> > ring_info_m1
){
    /*const std::vector<std::string> l0; const std::vector<std::string> l1;
    const std::vector<std::vector<int> > ring_info_m0;
    const std::vector<std::vector<int> > ring_info_m1;*/
    std::vector<std::vector<int> > ring_comp_m0;
    ring_comp_m0.resize(l0.size());
    

    // Iterate through rings and assign membership
    for (const std::vector<int>& r0 : ring_info_m0) {
        for (int atomIdx : r0) {
            ring_comp_m0[atomIdx] = std::vector<int>{-1};  // Store ring size as membership info
        }
    }
    

    for (const std::vector<int>& r0 : ring_info_m0) {
        std::string r0_label = "";
        for (int atomIdx : r0) {
            r0_label += l0.at(atomIdx);
        }
        std::string r0_label_rev = r0_label;  // Make a copy of the original string
        std::reverse(r0_label_rev.begin(), r0_label_rev.end());  // Reverse the copied string

        for (const std::vector<int>& r1 : ring_info_m0) {
            if(r0.size() == r1.size()){
                std::string r1_label = "";
                for (int atomIdx : r1) {
                    r1_label += l1.at(atomIdx);
                }

                std::vector<std::pair<std::string, int> > rotations = gen_rotations(r0_label);
                std::vector<std::pair<std::string, int> > inv_rotations = gen_rotations(r0_label_rev);
                
                std::vector<std::pair<std::string, int> > r0_rots;
                std::vector<std::pair<std::string, int> > r0_rev_rots;

                for(std::pair<std::string, int> rot : rotations){
                    if(! r1_label.compare(rot.first)){
                        r0_rots.push_back(rot);
                    }
                }
                for(std::pair<std::string, int> rot : inv_rotations){
                    if(! r1_label.compare(rot.first)){
                        r0_rev_rots.push_back(rot);
                    }
                }

                for(std::pair<std::string, int> rot : r0_rots){
                    for (int idx = 0; idx < r0.size(); ++idx) {
                        if(ring_comp_m0.at(r0.at(idx)).at(0) == -1){
                            ring_comp_m0.at(r0.at(idx)).at(0) = r1.at(idx - rot.second);
                        }else{
                            ring_comp_m0.at(r0.at(idx)).push_back(r1.at(idx - rot.second));
                        }
                    }
                }
                for(std::pair<std::string, int> rot : r0_rev_rots){
                    for (int idx = 0; idx < r0.size(); ++idx) {
                        int rev_idx = r1.size() - idx - 1;
                        if(ring_comp_m0.at(r0.at(idx)).at(0) == -1){
                            ring_comp_m0.at(r0.at(idx)).at(0) = r1.at(rev_idx - rot.second);
                        }else{
                            ring_comp_m0.at(r0.at(idx)).push_back(r1.at(rev_idx - rot.second));
                        }
                    }
                }
                return ring_comp_m0;
            }
        }

    }
    


}