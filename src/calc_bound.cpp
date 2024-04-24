#include <vector>
#include <algorithm> // for std::min
#include "test.hpp"

int calc_bound( std::vector<LabelClass>& label_classes) {
  int bound = 0;
  for ( LabelClass& label_class : label_classes) {
      bound += std::min(label_class.g.size(), label_class.h.size());
  }
  return bound;
}
