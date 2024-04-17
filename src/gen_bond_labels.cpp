#include <vector>
#include <algorithm> // for std::sort, std::unique

std::vector<float> gen_bond_labels(const std::vector<std::vector<float> >& g0, const std::vector<std::vector<float> >& g1) {
  
  std::vector<float> all_labels;
  all_labels.push_back(0.0);
  all_labels.push_back(1.0);
  all_labels.push_back(1.5);
  all_labels.push_back(2.0);
  all_labels.push_back(3.0);
  all_labels.push_back(4.0);
  all_labels.push_back(5.0);
  all_labels.push_back(6.0);
  // Vector to store potential bond labels (extracted from both matrices)
  std::vector<float> intersection;

  // Iterate over rows and columns of g0
  for (size_t i = 0; i < g0.size(); ++i) {
    for (size_t j = 0; j < g0[i].size(); ++j) {
      float current_label = g0[i][j];

      // Check if the label exists in each row of g1 (avoid nested loops)
      bool found_in_g1 = false;
      for (const std::vector<float>& row : g1) {
        if (std::find(row.begin(), row.end(), current_label) != row.end()) {
          found_in_g1 = true;
          break;
        }
      }

      // If found in g1 and not already in the intersection, add it
      if (found_in_g1 && (intersection.empty() || intersection.back() != current_label)) {
        intersection.push_back(current_label);
      }
    }
  }


  // Sort the intersection vector for desired output
  std::sort(intersection.begin(), intersection.end());
  // Use unique to remove consecutive duplicates (may leave gaps)
  intersection.erase(std::unique(intersection.begin(), intersection.end()), intersection.end());

  // Resize the vector to remove empty space from erasing (optional)
  intersection.resize(intersection.size());
  
  // Return the intersection vector containing common bond labels (sorted and unique)
  return intersection;
}
