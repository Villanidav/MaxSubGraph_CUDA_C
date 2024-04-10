#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"
#include "testFra.hpp"

using namespace std;

void print_label_info(const LabelClass* label) {
  if (label == nullptr) {
    std::cout << "No label selected." << std::endl;
    return;
  }

  std::cout << "Selected Label:" << std::endl;
  std::cout << "  g: ";
  for (int elem : label->g) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
  std::cout << "  h: ";
  for (int elem : label->h) {
    std::cout << elem << " ";
  }
  std::cout << std::endl;
  std::cout << "  adj: " << label->adj << std::endl;
  std::cout << "  label: " << label->label << std::endl;
}

int main() 
{
    //TEST GEN ROTATIONS
    /*
    std::string input = "example";
    std::vector<std::pair<std::string, int> >  result = gen_rotations(input);
 
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
    /*
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
    std::cout << std::endl;*/

    //TEST SELECT LABEL
      // Create some sample LabelClass objects
    std::vector<int> g1 = {1, 2, 3}, h1 = {4, 5};
    std::vector<int> g2 = {6, 7}, h2 = {8, 9, 10};
    std::vector<std::vector<int>> rings_g1, rings_g2; // Assuming rings_g is empty for this test

    LabelClass label_class1(g1, h1, rings_g1, 1, 0); // adj = 1, label = 0
    LabelClass label_class2(g2, h2, rings_g2, 0, 1); // adj = 0, label = 1

    // Create a list of LabelClass objects
    std::vector<LabelClass*> label_classes = {&label_class1, &label_class2};

    // Test with different map sizes
    int map_size1 = 0;
    int map_size2 = 2;

    // Call the function and print the selected label
    LabelClass* selected_label1 = select_label(label_classes, map_size1);
    std::cout << "Test with map_size = " << map_size1 << std::endl;
    print_label_info(selected_label1);

    LabelClass* selected_label2 = select_label(label_classes, map_size2);
    std::cout << "Test with map_size = " << map_size2 << std::endl;
    print_label_info(selected_label2);

    return 0;
}