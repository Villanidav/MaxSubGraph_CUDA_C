#include <vector>
#include <algorithm> // for std::min
#include "test.hpp"

int calc_bound( std::vector<LabelClass> label_classes) {
  int bound = 0;
  for ( const LabelClass& label_class : label_classes) {
      if( label_class.g.size() > label_class.h.size() ) bound = bound + label_class.h.size();
      else bound = bound + label_class.g.size();
  }
  return bound;
}
