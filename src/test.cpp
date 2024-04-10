#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"


int main() {
    std::string input = "example";
    std::vector<std::pair<std::string, int> >  result = gen_rotations(input);
 
    // Printing the result
    for (const auto& pair : result) {
        std::cout << pair.first << ", " << pair.second << std::endl;
    }

    return 0;
}