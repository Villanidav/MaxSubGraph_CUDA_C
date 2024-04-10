#include <iostream>
#include "test.hpp"


int main() {
    
    std::vector<int> g;
    std::vector<int> h;
    /*int adj;
    int label;*/
    std::vector<std::vector<int> > rings_g;

    g.push_back(8);
    g.push_back(10);
    g.push_back(2);
    
    h.push_back(45);
    h.push_back(40);
    h.push_back(25);
    
    rings_g.at(0).push_back(45);
    rings_g.at(0).push_back(40);

    rings_g.at(1).push_back(25);
    ;

    LabelClass label(g,h,rings_g,0,0);
    
    std::cout << "Stampo gli elementi di g ";
    for (int i : label.g)
    {
        std::cout << i << " ";
    }

    std::cout << "Stampo i match tra l'elemento con l'indice 8 all'interno della molecola ( il primo all'interno di g) ";
    std::vector<int> prova;
    prova.push_back(8);
    
    //label.get_ring_match_data(prova);
    
    

    
    return 0;
}