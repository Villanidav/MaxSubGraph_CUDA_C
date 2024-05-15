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



int main()
{

ROMol result;
std::string smile0 = "CC1CN(CC2CCCCC2COCC(=C)Cc3cc(F)c(F)cc3F)C1";
std::string smile1 = "C=C(COCC1CCCCC1CN2CCC2)Cc3ccccc3";


result = smiles_mcs(smile0, smile1, 1,1);

            
        
    



    return 0;

}
