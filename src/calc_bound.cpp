#include <vector>
#include <algorithm> // for std::min
#include "test.hpp"

int calc_bound(const std::vector<LabelClass*>& label_classes) {
  int bound = 0;
  for (const LabelClass* label_class : label_classes) {
    // Consider only labels present in both graphs (i.e., adj == 1)
    if (label_class->adj == 1) {
      bound += std::min(label_class->g.size(), label_class->h.size());
    }
  }
  return bound;
}
