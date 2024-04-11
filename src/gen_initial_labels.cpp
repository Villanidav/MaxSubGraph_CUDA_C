#include "test.hpp"
#include <vector>
#include <set>
using namespace std;
#include <algorithm>
#include <unordered_set>

std::vector<std::string> find_common_strings(const std::vector<std::string>& l0, const std::vector<std::string>& l1) {
  // Use a set to efficiently store and check for unique common strings
  std::set<std::string> common_strings;

  // Add elements from l0 to the set
  for (const std::string& str : l0) {
    common_strings.insert(str);
  }

  // Find common elements (intersection) and store them in the result vector
  std::vector<std::string> result;
  std::set_intersection(l1.begin(), l1.end(), common_strings.begin(), common_strings.end(),
                        std::back_inserter(result));

  return result;
}

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string> l0, const std::vector<std::string> l1,     std::vector<std::vector<int> > ring_classes){
    std::vector<LabelClass> label_classes;
    const std::vector<std::string> common_labels = find_common_strings(l0,l1);
    for(std::string com : common_labels){ std::cout << "\n" << com; }

    for (const std::string& label : common_labels) {
        // Filter atoms and ring data based on label
        std::vector<int> g_elems;
        std::vector<std::vector<int> > g_ring_classes;
        for (size_t i = 0; i < l0.size(); ++i) {
            if (l0[i] == label) {
                g_elems.push_back(i);
                g_ring_classes.push_back(ring_classes[i]); // Assuming ring_classes access by index
            }
        }

        std::vector<int> h_elems;
        for (size_t j = 0; j < l1.size(); ++j) {
            if (l1[j] == label) {
                h_elems.push_back(j);
            }
        }

        LabelClass label_tmp(g_elems,h_elems,g_ring_classes,0,0);
        label_classes.push_back(label_tmp);

    }
    return label_classes;
}