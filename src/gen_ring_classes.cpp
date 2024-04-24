#include <vector>
#include <string>
#include "test.hpp" // Include your header file where gen_rotations function is declared

std::vector<std::vector<int>> gen_ring_classes(const RDKit::RWMol& mol0, const RDKit::RWMol& mol1) {
    std::vector<std::string> l0, l1;
    for (const auto& atom : mol0.atoms()) {
        l0.push_back(atom->getSymbol());
    }
    for (const auto& atom : mol1.atoms()) {
        l1.push_back(atom->getSymbol());
    }

    std::vector<std::vector<int>> ring_info_m0, ring_info_m1;
    ring_info_m0 = mol0.getRingInfo()->atomRings();
    ring_info_m1 = mol1.getRingInfo()->atomRings();

    std::vector<std::vector<int>> ring_comp_m0(l0.size(), std::vector<int>{-1}); // Initialize with -1

    for (const std::vector<int>& r0 : ring_info_m0) {
        std::string r0_label;
        for (int atomIdx : r0) {
            r0_label += l0[atomIdx];
        }

        std::string r0_label_rev = r0_label;
        std::reverse(r0_label_rev.begin(), r0_label_rev.end());

        for (const std::vector<int>& r1 : ring_info_m1) {
            if (r0.size() == r1.size()) {
                std::string r1_label;
                for (int atomIdx : r1) {
                    r1_label += l1[atomIdx];
                }

                std::vector<std::pair<std::string, int>> rotations = gen_rotations(r0_label);
                std::vector<std::pair<std::string, int>> inv_rotations = gen_rotations(r0_label_rev);

                std::vector<std::pair<std::string, int>> r0_rots;
                std::vector<std::pair<std::string, int>> r0_rev_rots;

                for (const auto& rot : rotations) {
                    if (r1_label == rot.first) {
                        r0_rots.push_back(rot);
                    }
                }
                for (const auto& rot : inv_rotations) {
                    if (r1_label == rot.first) {
                        r0_rev_rots.push_back(rot);
                    }
                }

                for (const auto& rot : r0_rots) {
                    for (size_t idx = 0; idx < r0.size(); ++idx) {
                        int targetIdx = (idx - rot.second + r1.size()) % r1.size();
                        if (ring_comp_m0[r0[idx]][0] == -1) {
                            ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                        } else {
                            ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                        }
                    }
                }

                for (const auto& rot : r0_rev_rots) {
                    for (size_t idx = 0; idx < r0.size(); ++idx) {
                        int targetIdx = (r1.size() - idx - rot.second + r1.size()) % r1.size();
                        if (ring_comp_m0[r0[idx]][0] == -1) {
                            ring_comp_m0[r0[idx]][0] = r1[targetIdx];
                        } else {
                            ring_comp_m0[r0[idx]].push_back(r1[targetIdx]);
                        }
                    }
                }
            }
        }
    }
    return ring_comp_m0;
}
