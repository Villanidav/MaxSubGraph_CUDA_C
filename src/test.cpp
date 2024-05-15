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

 

     std::string filename = "smiles.txt";
    std::vector<std::pair<std::string, std::string>> smiles;
    std::ifstream file(filename);
    int skip;
    std::string first, second;
    if (file.is_open()) {
        std::string line;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            
            // Skip the first int
            iss >> skip;
            // Read the first string
            iss >> first;
            // Skip the next four ints
            for (int i = 0; i < 4; ++i) {
                iss >> skip;
            }
            // Read the second string
            iss >> second;
            // Store the pair of strings
            smiles.push_back(std::make_pair(first, second));
        }
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
    }

    std::cout << "Read " << smiles.size() << " pairs" << std::endl;

    // Print the pairs for verification
    for (const auto& pair : smiles) {
        std::cout << pair.first << " " << pair.second << std::endl;
    }





    // Close the file
    file.close();

  std::ofstream outfile("output.txt");
  std::streambuf* original_cout_buffer = std::cout.rdbuf();  // Save original buffer
  std::cout.rdbuf(outfile.rdbuf());
 
  


    clock_t start = clock();
    ROMol result;
    std::vector<std::string> result_string;
    
    for ( int i = 0 ; i<smiles.size(); ++i) {
     
            result = smiles_mcs(smiles.at(i).first, smiles.at(i).second, 1,1);

            result_string.clear();
           
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

    clock_t end = clock();

    // Calculate elapsed time in seconds
    double elapsed_seconds = (double)(end - start) / CLOCKS_PER_SEC;

    

    std::cout.rdbuf(original_cout_buffer);


    return 0;

}
