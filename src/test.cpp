//
// Created by davide on 4/19/24.
//
#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
#include <ctime>
using namespace std;
using namespace RDKit;



int main()
{
    string smile0 = "CC1CCc2c(C1)sc(NC(=O)c3ccc(cc3)C(C)(C)C)c2c4nnn(CC(=O)Nc5ccc6nc(oc6c5)c7ccccc7Cl)n4";
    string smile1 = "O=C(Cn1nnc(n1)c2c3CCCCc3sc2NC(=O)c4ccccc4)Nc5ccc6nc(oc6c5)c7ccccc7";
    clock_t start = clock();
    ROMol result;
    //cout<<"PRE FUNCTION" ;
    for ( int i = 0 ; i<1 ; ++i) {
        result = smiles_mcs(smile0, smile1, 1,1);
    }

    clock_t end = clock();

    // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    // Print the elapsed time in seconds
    std::cout << "\nElapsed time: " << elapsed_seconds << " seconds" << std::endl;


    std::vector<std::string> result_string;
    for (const auto &atom : result.atoms()) {
        result_string.push_back(atom->getSymbol());
    }

    cout << "\n MAX COMMON STRUCTURE ATOMS: \n[" ;
    for ( string idx : result_string )
        cout <<"'"<<idx<<"', ";
    cout << "]" ;

    return 0;
}
