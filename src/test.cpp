//
// Created by davide on 4/19/24.
//
#include <iostream>
#include <vector>
#include <string>
#include "test.hpp"

using namespace std;
using namespace RDKit;

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
    /*
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
    print_label_info(selected_label2);*/

    //TEST CALC BOUND
    /*
      // Create some sample LabelClass objects
    std::vector<int> g1 = {1, 2, 3}, h1 = {4, 5};
    std::vector<int> g2 = {6, 7}, h2 = {8, 9, 10};
    std::vector<std::vector<int>> rings_g1, rings_g2; // Assuming rings_g is empty for this test

    LabelClass label_class1(g1, h1, rings_g1, 1, 0); // adj = 1, label = 0, g size = 3, h size = 2
    LabelClass label_class2(g2, h2, rings_g2, 0, 1); // adj = 0, label = 1, g size = 2, h size = 3 (not considered)
    LabelClass label_class3(g1, h2, rings_g2, 1, 2); // adj = 1, label = 2, g size = 3, h size = 3

    // Create a list of LabelClass objects
    std::vector<LabelClass*> label_classes = {&label_class1, &label_class2, &label_class3};

    // Test calc_bound
    int bound = calc_bound(label_classes);
    std::cout << "Calculated bound: " << bound << std::endl;

    // Expected bound: min(g1, h1) + min(g3, h3) = 2 + 3 = 5 (considering only labels with adj=1)
    if (bound == 5) {
        std::cout << "Test passed!" << std::endl;
    } else {
        std::cout << "Test failed! Expected bound: 5" << std::endl;
    }*/

    //TEST GEN BOND LABELS
    /*

    // Test case 1: Empty intersection (assuming float values)
    std::vector<std::vector<float>> graph0 = {{3.8, 0.2, 2.1}, {0.8, 1.3, 1.7}, {2.4, 1.1, 0.9}};
    std::vector<std::vector<float>> graph1 = {{2.3, 1.8, 0.0}, {1.2, 2.7, 0.5}, {0.6, 0.1, 2.9}};
    std::vector<float> expected_labels = {}; // No common edges (adapt if expected is different)

    std::vector<float> test_labels = gen_bond_labels(graph0, graph1);

    if (test_labels == expected_labels) {
      std::cout << "Test passed! Expected labels (no common edges):";
    } else {
      std::cout << "Test failed! Got:";
    }
    for (float label : test_labels) {
      std::cout << " " << label;
    }
    std::cout << std::endl;

    // Test case 2: Common edges with duplicates (adapt values if needed)
    graph0 = {{1.0, 0.0}, {0.0, 1.0}};
    graph1 = {{1.0, 1.0}, {1.0, 0.0}};
    expected_labels = { 0.0 , 1.0 };

    test_labels = gen_bond_labels(graph0, graph1);

    if (test_labels == expected_labels) {
      std::cout << "Test passed! Expected labels:";
    } else {
      std::cout << "Test failed! Got:";
    }
    for (int label : test_labels) {
        std::cout << " " << label;
    }
    std::cout << std::endl;*/

    //TEST SMILES_MCS
    /*std::string s0 = "C[C@H](F)N";
    std::string s1 = "CN";
    std::pair<RWMol, RWMol> molecules = smiles_mcs(s0, s1, 1, 1);
    RWMol molA = molecules.first;  // Access the first molecule
    RWMol molB = molecules.second; // Access the second molecule

    for ( Atom *a : molA.atoms() ) {
        cout<<"\natomo prima molecola: "<< a->getSymbol();
    }
    for ( Atom *b : molB.atoms() ) {
        cout<<"\natomo seconda molecola: "<< b->getSymbol();
    }*/

    string smile0 = "C[C@H](F)N";
    string smile1 = "CN";
    ROMol result = smiles_mcs(smile0, smile1, 1,1);

    cout << "\n" << result.getNumAtoms();
    std::vector<std::string> result_string;
    for (const auto &atom : result.atoms()) {
        result_string.push_back(atom->getSymbol());
    }

    cout << "\n molecule: " <<endl;
    for ( string idx : result_string )
        cout << idx;


    return 0;
}