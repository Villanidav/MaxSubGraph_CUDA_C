#include "test.hpp"
#include <vector>
#include <set>
using namespace std;
#include <algorithm>
#include <unordered_set>


std::vector<std::string> find_common_strings(const std::vector<std::string>& l0, const std::vector<std::string>& l1) {
    // Utilizzare un set per memorizzare ed effettuare velocemente la ricerca di stringhe comuni uniche
    std::unordered_set<std::string> common_strings(l0.begin(), l0.end());

    // Vettore per memorizzare le stringhe comuni trovate
    std::vector<std::string> result;

    // Trovare le intersezioni tra le stringhe della seconda lista e le stringhe nel set
    for (const std::string& str : l1) {
        // Se la stringa è presente nel set delle stringhe comuni
        if (common_strings.find(str) != common_strings.end()) {
            // Aggiungila al risultato
            result.push_back(str);
            // Rimuovi la stringa dal set per evitare duplicati
            common_strings.erase(str);
        }
    }
       return result;
}

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string>& l0, const std::vector<std::string>& l1,     std::vector<std::vector<int> >& ring_classes){
    std::vector<LabelClass> label_classes;
    const std::vector<string> common_labels = find_common_strings(l0,l1);


    for (const std::string& label : common_labels) {
        // Filter atoms and ring data based on label
        std::vector<int> g_elems;
        std::vector<std::vector<int> > g_ring_classes;
        for (size_t i = 0; i < l0.size(); ++i) {
            if (l0[i] == label) {
                g_elems.push_back(i);
                if( !ring_classes.empty() ) g_ring_classes.push_back(ring_classes[i]); // Assuming ring_classes access by index

            }
        }

        std::vector<int> h_elems;
        for (size_t j = 0; j < l1.size(); ++j) {
            if (l1[j] == label) {
                h_elems.push_back(j);
            }
        }

        LabelClass label_tmp(g_elems,h_elems,g_ring_classes,0, label);
        label_classes.push_back(label_tmp);

    }
    return label_classes;
}