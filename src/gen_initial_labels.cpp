#include "test.hpp"
#include <vector>
#include <set>
using namespace std;
#include <algorithm>
#include <unordered_set>

std::vector<std::string> find_common_strings(const std::vector<std::string>& l0, const std::vector<std::string>& l1) {
    // Utilizzare un set per memorizzare ed effettuare velocemente la ricerca di stringhe comuni uniche
    std::set<std::string> common_strings;

    // Aggiungere le singole lettere dalla prima molecola al set
    for (const std::string& molecule : l0) {
        for (char letter : molecule) {
            std::string letter_str(1, letter); // Convertire il carattere in una stringa di lunghezza 1
            common_strings.insert(letter_str);
        }
    }

    // Trovare le intersezioni tra le singole lettere della seconda molecola e le stringhe nel set
    std::vector<std::string> result;
    for (const std::string& molecule : l1) {
        for (char letter : molecule) {
            std::string letter_str(1, letter); // Convertire il carattere in una stringa di lunghezza 1
            if (common_strings.find(letter_str) != common_strings.end()) {
                // Se la lettera è comune, aggiungila al risultato
                result.push_back(letter_str);
                // Rimuovi la lettera dal set per evitare duplicati
                common_strings.erase(letter_str);
            }
        }
    }

    return result;
}

std::vector<LabelClass> gen_initial_labels(const std::vector<std::string>& l0, const std::vector<std::string>& l1,     std::vector<std::vector<int> >& ring_classes){
    std::vector<LabelClass> label_classes;
    const std::vector<string> common_labels = find_common_strings(l0,l1);
    cout << "common labelssssaSAassaSasSAAS";
    for ( string r : common_labels )
        cout << " " << r;
    //for(std::string com : common_labels){ std::cout << "\n" << com; }

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