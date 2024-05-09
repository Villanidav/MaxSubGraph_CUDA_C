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
    string smile0 = "Nc1c2C(=O)Oc3ccccc3c2cc(c4ccc(F)cc4)c1c5nnn(CC(=O)Nc6ccc7nc(oc7c6)c8ccccc8Cl)n5";
    string smile1 = "O=C(Cn1nnc(n1)c2cc3C(=O)Oc4ccccc4c3cc2c5ccccc5)Nc6ccc7nc(oc7c6)c8ccccc8";
    //string smile0 = "Cc1ccccc1c2nc(SCc3ccccc3)c(c4CCCc24)c5nnn(CC(=O)Nc6ccc7nc(oc7c6)c8ccccc8Cl)n5";
    //string smile1 = "O=C(Cn1nnc(n1)c2c3CCCc3c(nc2SCc4ccccc4)c5ccccc5)Nc6ccc7nc(oc7c6)c8ccccc8";
    clock_t start = clock();
    ROMol result;
    //cout<<"PRE FUNCTION" ;
    result = smiles_mcs(smile0, smile1, 1,1);


    clock_t end = clock();

    // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    // Print the elapsed time in seconds
    std::cout << "Elapsed time: " << elapsed_seconds << " seconds" << std::endl;

/*
    std::vector<std::string> result_string;
    for (const auto &atom : result.atoms()) {
        result_string.push_back(atom->getSymbol());
    }

    cout << "\n MAX COMMON STRUCTURE ATOMS: \n[" ;
    for ( string idx : result_string )
        cout <<"'"<<idx<<"', ";
    cout << "]" ;
*/
    return 0;
}
