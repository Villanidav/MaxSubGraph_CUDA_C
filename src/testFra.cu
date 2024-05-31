//
// Created by davide on 4/19/24.
//

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
#include <ctime>
using namespace std;
using namespace RDKit;

void printMol(ROMol result) {
    std::vector<std::string> result_string;

    for (const auto &atom : result.atoms()) {
        result_string.push_back(atom->getSymbol());
    }

    cout << "[";
    for ( int idx = 0; idx < result_string.size(); idx++ ){
        if(idx == result_string.size()-1 ){
            cout <<"'"<<result_string.at(idx)<<"']"<<endl;
        }
        else cout <<"'"<<result_string.at(idx)<<"', ";
    }
}

int main(){

    ROMol result;
    std::string smile0 = "CN([C@@H]1C[C@@H](C1)C(=O)N2CCC2)c3ccc(cc3C(=O)n4cc(C(=O)O)c5ccccc45)[N+](=O)[O-]";
    std::string smile1 = "O=C([C@@H]1C[C@@H](C1)Nc2ccccc2C(=O)n3ccc4ccccc34)N5CCC5";
    cudaFree(0);
    clock_t start = clock();
    result = smiles_mcs(smile0, smile1, 1,1);
    clock_t end = clock();
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;
    cout<<"\nElapsed time: [ " << elapsed_seconds << " ]"<<endl;
    printMol(result);

    return 0;
}
