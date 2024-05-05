#include <vector>
#include "test.hpp"

LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size) {

    
    int min_size = 999;
    LabelClass* label = nullptr;

    for ( LabelClass* c_label : label_classes) {
        if (c_label->adj == 1 || map_size == 0) {
            int c_max_size = std::max(c_label->g.size(), c_label->h.size());
            if (c_max_size < min_size) {
                min_size = c_max_size;
                label = c_label;
            }
        }
    }

    return label;
}