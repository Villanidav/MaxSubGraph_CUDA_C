#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
using namespace std;

int main() 
{
    //TEST GEN ROTATIONS
    /*
    std::string input = "example";
    std::vector<std::pair<std::string, int>>  result = gen_rotations(input);
 
    // Printing the result
    for (const auto& pair : result) {
        std::cout << pair.first << ", " << pair.second << std::endl;
    }*/

    //TEST SELECT VERTEX
    /* 
    std::vector<int> vtx_set = {0, 1, 2}; // Example vertex set
    std::vector<std::vector<int>> g = {{0, 1, 0}, {1, 1, 0}, {1, 0, 0}, }; // Example graph adjacency matrix

    int selected_vertex = select_vertex(vtx_set, g);
    std::cout << "Selected vertex: " << selected_vertex << std::endl;*/

    //TEST HOOD
        // Example graph represented as an adjacency matrix
    std::vector<std::vector<int>> g = {{0, 1, 0, 1},
                                        {1, 0, 1, 0},
                                        {0, 1, 0, 1},
                                        {1, 0, 1, 0}};
    
    // Example vertex and bond type
    int vtx = 1;  // Example vertex
    int edge = 1; // Example bond type

    // Print the hood of the given vertex
    std::cout << "Hood of vertex " << vtx << " with bond type " << edge << ": ";
    std::vector<int> neighbors = hood(vtx, g, edge);
    for (int neighbor : neighbors) {
        std::cout << neighbor << " ";
    }
    std::cout << std::endl;


    return 0;
}