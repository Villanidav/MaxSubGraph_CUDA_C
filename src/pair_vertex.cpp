//
// Created by davide on 5/3/24.
//
#include "test.hpp"
#include <vector>
#include <string>
#include <queue>
#include <algorithm>

std::pair<int, int> pair_vertex(LabelClass lb,  std::vector<std::vector<float>> g0) {
    std::pair<int, int> pair;
    int v = select_vertex(lb.g, g0);
    pair.first = v;
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = {};
    std::vector<std::vector<int> > returned_data = {};
    returned_data = lb.get_ring_match_data(vector);

    if ( !returned_data.empty()) {
        v_ring_atoms = lb.get_ring_match_data(vector).at(0);
    }

    for ( int w : lb.h) {
        if( !v_ring_atoms.empty() ) {
            for(int x : v_ring_atoms){

                if( x == -1 ) {
                    break;
                }
                if( x == w ){
                    pair.second = w;
                    return pair;
                }
            }
        }
        pair.second = w;
        return pair;
    }
    pair.second = -1;
    return pair;
}