#include <vector>
#include <algorithm> // for std::sort, std::unique

std::vector<int> gen_bond_labels(const std::vector<std::vector<int> >& g0, const std::vector<std::vector<int> >& g1) {
  // Vector to store potential bond labels (extracted from both matrices)
  std::vector<int> all_labels;

  // Iterate over rows and columns of both matrices to extract edges
  for (size_t i = 0; i < g0.size(); ++i) {
    for (size_t j = 0; j < g0[i].size(); ++j) {
      all_labels.push_back(g0[i][j]);
    }
  }
  for (size_t i = 0; i < g1.size(); ++i) {
    for (size_t j = 0; j < g1[i].size(); ++j) {
      all_labels.push_back(g1[i][j]);
    }
  }

  // Sort the vector to group potential common edges
  std::sort(all_labels.begin(), all_labels.end());

  // Use unique to remove consecutive duplicates within all_labels
  std::vector<int>::iterator it = std::unique(all_labels.begin(), all_labels.end());

  // Resize the vector to remove elements after the unique iterator
  all_labels.resize(std::distance(all_labels.begin(), it));

  // Intersection vector to store common bond labels
  std::vector<int> intersection;

  // Iterate through the unique elements in all_labels
  for (int label : all_labels) {
    // Check if the label exists in each row of g1 (avoiding nested loops)
    bool found_in_g1 = false;
    for (const std::vector<int>& row : g1) {
      if (std::find(row.begin(), row.end(), label) != row.end()) {
        found_in_g1 = true;
        break;
      }
    }
    // If found in g1, add it to the intersection vector
    if (found_in_g1) {
      intersection.push_back(label);
    }
  }

  // Return the intersection vector containing common bond labels
  return intersection;
}