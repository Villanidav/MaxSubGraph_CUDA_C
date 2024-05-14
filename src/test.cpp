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

    /*std::string filename = "smiles.txt";
    std::vector<std::string> smiles;
    std::ifstream file(filename);

    if (file.is_open()) {
    std::string line;
        while (std::getline(file, line)) {
            smiles.push_back(line);
        }

    } else {
    std::cerr << "Error opening file: " << filename << std::endl;
  }

  

    // Close the file
    file.close();

  std::ofstream outfile("output.txt");
  std::streambuf* original_cout_buffer = std::cout.rdbuf();  // Save original buffer
  std::cout.rdbuf(outfile.rdbuf());
    

    

    clock_t start = clock();
    ROMol result;
    //cout<<"PRE FUNCTION" ;
    for ( int i = 0 ; i<smiles.size() -1; ++i) {
        for(int j = i+1; j < smiles.size(); j++){
            result = smiles_mcs(smiles.at(i), smiles.at(j), 1,1);

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
    }

    clock_t end = clock();

    // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    // Print the elapsed time in seconds
    std::cout << "\nElapsed time: " << elapsed_seconds << " seconds" << std::endl;


std::cout.rdbuf(original_cout_buffer);*/
    string s0 = "COCCCOc1cc(C[C@@H](C[C@H](NC(=O)OC(C)(C)C)C(O)CCCC(=O)N2CC3CCC(C3)C2)C(C)C)ccc1OC";
    string s1 = "O=C(CCCCCCCCc1ccccc1)N2CC3CCC(C3)C2";
    ROMol result = smiles_mcs(s0, s1 );
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
    cout<<"done\n\n";
    return 0;
}
