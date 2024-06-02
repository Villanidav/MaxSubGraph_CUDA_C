//
// Created by davide on 5/3/24.
//
#include "test.hpp"
#include <vector>
#include <string.h>
#include <string>
#include <queue>
#include <algorithm>
#include <unordered_set>

using namespace std;
const int DIM_POOL = 8;
std::vector<std::vector<float>> g0;
std::vector<std::vector<float>> g1;
std::vector<float> edge_labels;
//vector<pair<int,int>> m_best;

__device__ float *gpu_edge_labels;
__device__ int size_edge_labels;


__device__ float **gpu_g0;
__device__ int size_gpu_g0_row;
__device__ int size_gpu_g0_col;

__device__ float **gpu_g1;
__device__ int size_gpu_g1_row;
__device__ int size_gpu_g1_col;

typedef struct{
    int g_size;
    int h_size;
    int row_ring_size;
    int *col_ring_size;
    int *g;
    int *h;
    int adj;
    char label[4];
    int **rings_g;
}GpuLabelClass;


typedef struct{
    int first;
    int second;
}Pair;


typedef struct {
    int labels_size;
    int m_size;
    GpuLabelClass *labels;
    GpuLabelClass single_label;
    Pair *m_local;
}ThreadVar;


__shared__ Pair *m_best;
__shared__ int m_best_size;

LabelClass *select_label_DC(std::vector<LabelClass*>& label_classes, int map_size) {

<<<<<<< HEAD


void copyIntArray(int *a, int *b, int sizeb){
    for ( int i = 0 ; i < sizeb ; i++){
        a[i] = b[i];
    }
}


void copyIntMatrix(int **a, int **b, int rowsize, int *colsize )
{
    for( int i = 0 ; i < rowsize ; i++){
        for( int j = 0 ; j < colsize[i] ; j++){
            a[i][j] = b[i][j];
        }
    }
=======
    
    int min_size = 999;
    LabelClass* label = nullptr;

    for ( LabelClass* c_label : label_classes) {
        if (c_label->adj == 1 || map_size == 0) {
            int c_max_size = std::max(c_label->g.size(), c_label->h.size());
            if (c_max_size < min_size) {
                min_size = c_max_size;
                label = c_label;
            }
        }
    }

    return label;
>>>>>>> origin/main
}
// vtx_set: selected label class
// g: selected graph
int select_vertex(std::vector<int>& vtx_set, std::vector<std::vector<float> >& g) {
    // selects node from graph given a label, choosing an adjacent node with the maximum degree

    int max_deg = -1;
    int vtx = 0;
    for (int c_vtx : vtx_set) {
        int deg = 0;
        for (float i : g[c_vtx]) {
            if (i != 0) {
                deg += 1;
            }
        }

        if (deg > max_deg) {
            max_deg = deg;
            vtx = c_vtx;
        }
    }
    return vtx;
}
std::vector<std::string> find_common_strings(const std::vector<std::string>& l0, const std::vector<std::string>& l1) {
    // Utilizzare un set per memorizzare ed effettuare velocemente la ricerca di stringhe comuni uniche
    std::unordered_set<std::string> common_strings(l0.begin(), l0.end());

    // Vettore per memorizzare le stringhe comuni trovate
    std::vector<std::string> result;

    // Trovare le intersezioni tra le stringhe della seconda lista e le stringhe nel set
    for (const std::string& str : l1) {
        // Se la stringa è presente nel set delle stringhe comuni
        if (common_strings.find(str) != common_strings.end()) {
            // Aggiungila al risultato
            result.push_back(str);
            // Rimuovi la stringa dal set per evitare duplicati
            common_strings.erase(str);
        }
    }
    return result;
}
std::vector<LabelClass> cpu_gen_initial_labels(const std::vector<std::string>& l0, const std::vector<std::string>& l1,     std::vector<std::vector<int> >& ring_classes){
    std::vector<LabelClass> label_classes;
    const std::vector<string> common_labels = find_common_strings(l0,l1);


    for (const std::string& label : common_labels) {
        // Filter atoms and ring data based on label
        std::vector<int> g_elems;
        std::vector<std::vector<int> > g_ring_classes;
        for (size_t i = 0; i < l0.size(); ++i) {
            if (l0[i] == label) {
                g_elems.push_back(i);
                if( !ring_classes.empty() ) g_ring_classes.push_back(ring_classes[i]); // Assuming ring_classes access by index

            }
        }

        std::vector<int> h_elems;
        for (size_t j = 0; j < l1.size(); ++j) {
            if (l1[j] == label) {
                h_elems.push_back(j);
            }
        }

        LabelClass label_tmp(g_elems,h_elems,g_ring_classes,0, label);
        label_classes.push_back(label_tmp);

    }
    return label_classes;
}
std::vector<float> cpu_gen_bond_labels(const std::vector<std::vector<float> >& g0, const std::vector<std::vector<float> >& g1) {

    std::vector<float> all_labels;
    all_labels.push_back(0.0);
    all_labels.push_back(1.0);
    all_labels.push_back(1.5);
    all_labels.push_back(2.0);
    all_labels.push_back(3.0);
    all_labels.push_back(4.0);
    all_labels.push_back(5.0);
    all_labels.push_back(6.0);
    // Vector to store potential bond labels (extracted from both matrices)
    std::vector<float> intersection;

    // Iterate over rows and columns of g0
    for (size_t i = 0; i < g0.size(); ++i) {
        for (size_t j = 0; j < g0[i].size(); ++j) {
            float current_label = g0[i][j];

            // Check if the label exists in each row of g1 (avoid nested loops)
            bool found_in_g1 = false;
            for (const std::vector<float>& row : g1) {
                if (std::find(row.begin(), row.end(), current_label) != row.end()) {
                    found_in_g1 = true;
                    break;
                }
            }

            // If found in g1 and not already in the intersection, add it
            if (found_in_g1 && (intersection.empty() || intersection.back() != current_label)) {
                intersection.push_back(current_label);
            }
        }
    }


    // Sort the intersection vector for desired output
    std::sort(intersection.begin(), intersection.end());
    // Use unique to remove consecutive duplicates (may leave gaps)
    intersection.erase(std::unique(intersection.begin(), intersection.end()), intersection.end());

    // Resize the vector to remove empty space from erasing (optional)
    intersection.resize(intersection.size());

    // Return the intersection vector containing common bond labels (sorted and unique)
    return intersection;
}
std::vector<int> hood(int vtx, const std::vector<std::vector<float>>& g, float edge) {
    // Return the neighbors of a specified node, with the specified bond type.
    std::vector<int> friends;
    for (std::size_t i = 0; i < g.size(); ++i) {
        if (g[i][vtx] == edge && static_cast<std::size_t>(vtx) != i) {
            friends.push_back(i);
        }
    }

    return friends;
}




__device__
bool contains(int value, int *arr, int size) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return true;
        }
    }
    return false;
}


bool host_contains(int value, int *arr, int size) {
    for (int i = 0; i < size; ++i) {
        if (arr[i] == value) {
            return true;
        }
    }
    return false;
}

//      puts into a 2D array the data regarding indexes of rings related to the array of elements
//      2D-array that will contain the result that will be modified
//      1D-array containing idxList
//      1D-array of elements
//      int size of elements
__device__
int kernel_get_ring_match_data(int *dim_col, int **result, int *idxList, int *elems, int elem_size, GpuLabelClass *lc){
    
    int index;
    int k=0;
    for( int i = 0; i < elem_size ; ++i){
        index = 0;
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) idxList[k] = index;
            index++;
            k++;
        }
    }

    for ( int k = 0 , i = 0 ; k < index ; ++k , ++i ){
        dim_col[i] = lc->col_ring_size[idxList[k]];
        for ( int j = 0 ; j < dim_col[i] ; ++j ){
             result[i][j] = lc->rings_g[idxList[k]][j];
        }
    }
    return index;
}

int host_get_ring_match_data(int *dim_col, int **result, int *idxList, int *elems, int elem_size, GpuLabelClass *lc){
    
    int index;
    int k=0;
    for( int i = 0; i < elem_size ; ++i){
        index = 0;
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) idxList[k] = index;
            index++;
            k++;
        }
    }

    for ( int k = 0 , i = 0 ; k < index ; ++k , ++i ){
        dim_col[i] = lc->col_ring_size[idxList[k]];
        for ( int j = 0 ; j < dim_col[i] ; ++j ){
             result[i][j] = lc->rings_g[idxList[k]][j];
        }
    }
    return index;
}
// function used in matchable, returns a 1D array containing the ring data related to vertex v
// returns the size of the 1D array
__device__
int kernel_matchable_ring_data(int **result, int v, GpuLabelClass *lc, int *idxList, int *dim_col){
    int index;
    int k=0;
    index = 0;
    for( int j = 0 ; j < lc->g_size ; ++j ){
        if( lc->g[j] == v ) idxList[k] = index;
        index++;
        k++;
    }


    for ( int k = 0 , i = 0 ; k < index ; ++k , ++i ){
        dim_col[i] = lc->col_ring_size[idxList[k]];
        for ( int j = 0 ; j < dim_col[i] ; ++j ){
             result[i][j] = lc->rings_g[idxList[k]][j];
        }
    }
    return index;
}


int host_matchable_ring_data(int *result, int v, GpuLabelClass *lc){
    int index=0;
    for ( int i = 0 ; i < lc->g_size ; ++i){
        if ( lc->g[i] == v ) result[index] = index;
        index++;
    }
    return index;

}

// return the best select label given an array of labels
__device__
void kernel_select_label(GpuLabelClass *label , GpuLabelClass *lcs, int map_size, int lcs_size){
    int min = 999;
    int max = 0;
    for( int i = 0 ; i < lcs_size ; ++i ){
        if( lcs[i].adj == 1 || map_size == 0 ){
            if( lcs[i].g_size > lcs[i].h_size ) max = lcs[i].g_size;
            else max = lcs[i].h_size;
            if( max < min ){
                min = max;
                *label = lcs[i];
            }
        }
    }
    return;
}


void host_select_label(GpuLabelClass *label , GpuLabelClass *lcs, int map_size, int lcs_size){
    int min = 999;
    int max = 0;
    for( int i = 0 ; i < lcs_size ; ++i ){
        if( lcs[i].adj == 1 || map_size == 0 ){
            if( lcs[i].g_size > lcs[i].h_size ) max = lcs[i].g_size;
            else max = lcs[i].h_size;
            if( max < min ){
                min = max;
                *label = lcs[i];
            }
        }
    }
    return;
}

// compute the bound given a 1D array of struct GpuLabelClass and its size
__device__
int kernel_calc_bound(GpuLabelClass *lcs, int lc_size) {
    int bound = 0;
    for( int i = 0 ; i < lc_size ; ++i){
        if ( lcs[i].g_size > lcs[i].h_size ) bound = bound + lcs[i].h_size;
        else bound = bound + lcs[i].g_size;
    }
    return bound;
}


//return = size of the friends
//friend is the OUTPUT
__device__
int hoodG(int *friends,int vtx, float edge, float **g0) {
    int size = 0;
    for (int i = 0; i < size_gpu_g0_row; i++) {
        if (g0[i][vtx] == edge && vtx != i) {
            friends[size] = i;
            size++;
        }
    }
    return size;
}


int host_hoodG(int *friends,int vtx, float edge, float **g0) {
    int size = 0;
    for (int i = 0; i < size_gpu_g0_row; i++) {
        if (g0[i][vtx] == edge && vtx != i) {
            friends[size] = i;
            size++;
        }
    }
    return size;
}

// result == size of generated label
// output is l_draft
// input : v
__device__
void kernel_resize(int *array, int size_arr, int place_availabel){
    int count = 0;
    bool flag = true;
    for(int i = 0; i < size_arr && flag ; ++i){
        if(array[i] == -1){continue;}
        if( i == count){ count ++; continue; }
        array[count] = array[i]; count++;
        if(count == place_availabel){flag = false;}
    }
}

// result == size of generated label
// output is l_draft
// input : v

void host_resize(int *array, int size_arr, int place_availabel){
    int count = 0;
    bool flag = true;
    for(int i = 0; i < size_arr && flag ; ++i){
        if(array[i] == -1){continue;}
        if( i == count){ count ++; continue; }
        array[count] = array[i]; count++;
        if(count == place_availabel){flag = false;}
    }
}


__device__
int kernel_gen_new_Labels(GpuLabelClass *l_draft ,  int v, int w, GpuLabelClass *lcs, int lcs_size, int *idxList) {
    int vs,ws, draft_size = 0;
    int dim_row;
    int count = 0;
    for ( int i = 0 ; i < lcs_size ; ++i ){
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize = hoodG(l_draft[draft_size].g, v, gpu_edge_labels[j], gpu_g0 );
            for ( int k = 0 , vs = 0 ; k < friendsize ; ++k ){
                if( contains(l_draft[draft_size].g[k], lcs[i].g, lcs[i].g_size) ){  vs++;  }
                else{ l_draft[draft_size].g[k] = -1;}
            }
            kernel_resize(l_draft[draft_size].g, friendsize, vs );

            dim_row = kernel_get_ring_match_data(l_draft[draft_size].col_ring_size, l_draft[draft_size].rings_g, idxList ,l_draft[draft_size].g, vs+1, &lcs[i] );

            friendsize = hoodG(l_draft[draft_size].h, w, gpu_edge_labels[j], gpu_g1 );
            for ( int k = 0 , ws = 0 ; k < friendsize ; ++k ){
                if( contains(l_draft[draft_size].h[k], lcs[i].h, lcs[i].h_size) ){  ws++; }
                else {
                    l_draft[draft_size].h[k] = -1;
                }
            }
            kernel_resize(l_draft[draft_size].h, friendsize, ws );

            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) adj = 1;
                else adj = 0;
                l_draft[draft_size].g_size = vs;
                l_draft[draft_size].h_size = ws;
                l_draft[draft_size].row_ring_size = dim_row;
                l_draft[draft_size].adj = adj;
                for( int c = 0 ; c < 4 ; c++){
                    l_draft[draft_size].label[c] = lcs[i].label[c];
                }
                
                draft_size++;
            }
        }
    }
    return draft_size;
}


int host_gen_new_Labels(GpuLabelClass *l_draft ,  int v, int w, GpuLabelClass *lcs, int lcs_size, int *idxList) {
    int vs,ws, draft_size = 0;
    int dim_row;
    int count = 0;
    for ( int i = 0 ; i < lcs_size ; ++i ){
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize = host_hoodG(l_draft[draft_size].g, v, gpu_edge_labels[j], gpu_g0 );
            for ( int k = 0 , vs = 0 ; k < friendsize ; ++k ){
                if( host_contains(l_draft[draft_size].g[k], lcs[i].g, lcs[i].g_size) ){  vs++;  }
                else{ l_draft[draft_size].g[k] = -1;}
            }
            host_resize(l_draft[draft_size].g, friendsize, vs );

            dim_row = host_get_ring_match_data(l_draft[draft_size].col_ring_size, l_draft[draft_size].rings_g, idxList ,l_draft[draft_size].g, vs+1, &lcs[i] );

            friendsize = host_hoodG(l_draft[draft_size].h, w, gpu_edge_labels[j], gpu_g1 );
            for ( int k = 0 , ws = 0 ; k < friendsize ; ++k ){
                if( host_contains(l_draft[draft_size].h[k], lcs[i].h, lcs[i].h_size) ){  ws++; }
                else {
                    l_draft[draft_size].h[k] = -1;
                }
            }
            host_resize(l_draft[draft_size].h, friendsize, ws );

            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) adj = 1;
                else adj = 0;
                l_draft[draft_size].g_size = vs;
                l_draft[draft_size].h_size = ws;
                l_draft[draft_size].row_ring_size = dim_row;
                l_draft[draft_size].adj = adj;
                for( int c = 0 ; c < 4 ; c++){
                    l_draft[draft_size].label[c] = lcs[i].label[c];
                }
                
                draft_size++;
            }
        }
    }
    return draft_size;
}

// ------ FRA functions //


//helper function
void vectorToPointerEdge(float *gpu_edge_labels){
    if(edge_labels.size() == 0){
        gpu_edge_labels = nullptr;
        return;
    }
    int size = 0;
    for(float edg : edge_labels){
        gpu_edge_labels[size] = edg;
        size++;
    }

    return;
}


void vectorToPointerMatrix(const std::vector<std::vector<float>>& g,
                           float** gpu_g
) {
    // Get dimensions of the vector
    int numRows = g.size();
    if (numRows == 0) {
        // Empty vector, set pointers to nullptr
        gpu_g = nullptr;
        return;
    }

    int numCol = g[0].size();

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCol; ++j) {
            gpu_g[i][j] = g[i][j];
        }
    }
}


void LabelFromCpuToGpu(GpuLabelClass *new_label, const vector<LabelClass>& old_label ){

    for (int idx = 0 ; idx < old_label.size() ; ++idx ){
        new_label[idx].g_size = old_label.at(idx).g.size();
        new_label[idx].h_size = old_label.at(idx).h.size();
        new_label[idx].row_ring_size = old_label.at(idx).rings_g.size();
        new_label[idx].adj = old_label.at(idx).adj;
        strcpy(new_label[idx].label, old_label.at(idx).label.c_str() );

        for ( int j = 0 ; j < old_label.at(idx).g.size() ; ++j ){
            new_label[idx].g[j] = old_label.at(idx).g.at(j);
        }
        for ( int j = 0 ; j < old_label.at(idx).h.size() ; ++j ){
            new_label[idx].h[j] = old_label.at(idx).h.at(j);}

        for ( int row = 0 ; row < old_label.at(idx).rings_g.size() ; ++row ){
            new_label[idx].col_ring_size[row] = old_label.at(idx).rings_g.at(row).size();
            for ( int col = 0 ; col < old_label.at(idx).rings_g.at(row).size() ; ++col ){
                new_label[idx].rings_g[row][col] = old_label.at(idx).rings_g.at(row).at(col);
            }
        }

    }
}


__device__ bool kernel_matchable(int**  v_ring_atoms, int v, int w, GpuLabelClass *lc) {
    /*int size  = kernel_matchable_ring_data(v_ring_atoms, v, lc );
    if( size > 0  ) {
        for(int i = 0; i< size; i++){
            if( v_ring_atoms[i] == -1 )return false;
            if( v_ring_atoms[i] == w ) return true;
        }
        return false;
    }*/
    return true;
}

bool host_matchable(int*  v_ring_atoms, int v, int w, GpuLabelClass *lc) {
    int size  = host_matchable_ring_data(v_ring_atoms, v, lc );
    if( size > 0  ) {
        for(int i = 0; i< size; i++){
            if( v_ring_atoms[i] == -1 )return false;
            if( v_ring_atoms[i] == w ) return true;
        }
        return false;
    }
    return true;
}



// vtx_set: selected label class
// g: selected graph
__device__ void kernel_select_vertex(int *result, int *result_pos, int *vtx_set, int vtx_size, float **g, int num_row, int num_column) {
    int max_deg = -1;
    int vtx = 0;

    for(int i = 0; i < vtx_size; i++){
        int deg = 0;
        for(int j = 0; j < num_column; j++){
            int consider = g[vtx_set[i]][j];
            if(consider != 0){
                deg++;
            }
        }

        if(deg>max_deg){
            max_deg = deg;
            *result = vtx_set[i];
            *result_pos = i;
        }
    }

}


// vtx_set: selected label class
// g: selected graph
void host_select_vertex(int *result, int *result_pos, int *vtx_set, int vtx_size, float **g, int num_row, int num_column) {
    int max_deg = -1;
    int vtx = 0;

    for(int i = 0; i < vtx_size; i++){
        int deg = 0;
        for(int j = 0; j < num_column; j++){
            int consider = g[vtx_set[i]][j];
            if(consider != 0){
                deg++;
            }
        }

        if(deg>max_deg){
            max_deg = deg;
            *result = vtx_set[i];
            *result_pos = i;
        }
    }
    return ;
}







size_t calcSize(const vector<LabelClass>& lcs) {
    size_t result=0;
    for( const LabelClass& lc : lcs) {
        result = result + (lc.g.size()+lc.h.size());
        for( vector<int> i : lc.rings_g ) {
            result += i.size();
        }
        result += 3;
    }
    return result;
}

void printLabelClass(const LabelClass& lb) {
        cout<< lb.label << " [ ";
        cout<< " G("<< lb.g.size() << "): ";
        if(!lb.g.empty()) {for ( int i : lb.g ) cout<<"["<<i<<"]";}
        cout<< " H("<< lb.h.size() << "): ";
        if(!lb.h.empty()) {for ( int i : lb.h ) cout<<"["<<i<<"]";}
        cout<< " RINGS("<< lb.rings_g.size() << "): [";
        for( vector<int> i : lb.rings_g ){cout<<"("<<i.size()<<")"<<"["; for( int j: i ) cout<<j<<", ";  cout<<" ]";}
        cout<<"]";
        cout<< " edge : " <<lb.adj<<" " ;
        cout<< lb.label << " ] "<<endl;
}




LabelClass *select_label(std::vector<LabelClass*>& label_classes, int map_size);

bool matchable(const int v,const int w, LabelClass lc ) {
    std::vector<int> vector;
    vector.push_back(v);
    std::vector<int>  v_ring_atoms = {};
    v_ring_atoms = lc.get_ring_match_data(vector).at(0);

    if( !v_ring_atoms.empty() ) {
        for(const int x : v_ring_atoms){
            if( x == -1 )return false;
            if( x == w ) return true;}
        return false;
    }
    return true;
}

vector<LabelClass> genNewLabels(int v, int w, const vector<LabelClass>& lcs) {
    vector<LabelClass> l_draft;

    for(LabelClass label : lcs){
        for(float edge_l : edge_labels){
            std::vector<int> v_conn;
            std::vector<int> w_conn;
            std::vector<std::vector<int> > v_c_rings;

            for(int vtx : hood(v,g0,edge_l)){if( std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() )   v_conn.push_back(vtx);}

            v_c_rings = label.get_ring_match_data(v_conn);

            for(int vtx : hood(w,g1,edge_l)){if(std::find(label.h.begin(),label.h.end(),vtx) != label.h.end() )  w_conn.push_back(vtx);}

            int adj;
            if(!v_conn.empty() && !w_conn.empty()){
                if(edge_l != 0.0 || label.adj == 1) adj = 1;
                else adj = 0;
                LabelClass tmp(v_conn,w_conn,v_c_rings,adj, label.label);
                l_draft.push_back(tmp);
            }
        }
    }
    return l_draft;
}







__global__ 
void kernel_function( ThreadVar *thread_pool_read, ThreadVar *thread_pool_write ){
    int index = threadIdx.x + blockIdx.x * blockDim.x;
    int space = 4;
    GpuLabelClass *label = nullptr;

    kernel_select_label(label, thread_pool_read[index].labels, 
                        thread_pool_read[index].m_size, thread_pool_read[index].labels_size);

    if ( (thread_pool_read[index].m_size + kernel_calc_bound(thread_pool_read[index].labels, thread_pool_read[index].labels_size  )  
        < m_best_size) || !label  ) return;

    int jump = 0;
    for( int v_idx = 0 ; v_idx < label->g_size ; ++v_idx){
        for( int w_idx = 0 ; w_idx < label->h_size ; ++w_idx){
            if( !kernel_matchable(label->rings_g, label->g[v_idx], label->h[w_idx], label )) continue;

            for( int z = 0 ; z < thread_pool_write[4*index + jump].m_size ; ++z ){
                thread_pool_write[4*index + jump].m_local[z].first = thread_pool_read[index].m_local[z].first;
                thread_pool_write[4*index + jump].m_local[z].second = thread_pool_read[index].m_local[z].second;
            }

            /*kernel_gen_new_labels( thread_pool_write[4*index + jump].labels, label.g[v_idx], label.h[w_idx] ,
                                   thread_pool_read[index].labels , thread_pool_read[index].labels_size ,
                                   );*/

            if( thread_pool_write[4*index + jump].m_size > m_best_size ){
                for( int z = 0 ; z < thread_pool_write[4*index + jump].m_size ; ++z ){
                    m_best[z].first = thread_pool_write[4*index + jump].m_local[z].first;
                    m_best[z].second = thread_pool_write[4*index + jump].m_local[z].second;
                }
            }
        }
    }
}






void cpyThreadPool( ThreadVar *thread_pool_read, ThreadVar *thread_pool_write ){
    int r_idx = 0;
    for( int w_idx = 0; w_idx < DIM_POOL*DIM_POOL ; w_idx++){ 
        thread_pool_read[r_idx].labels_size = thread_pool_write[w_idx].labels_size;
        //da modificare l'if
        if ( thread_pool_write[w_idx].labels_size > 0 ){
            thread_pool_read[r_idx].m_size = thread_pool_write[w_idx].m_size;
            thread_pool_read[r_idx].m_local->first = thread_pool_write[w_idx].m_local->first;
            thread_pool_read[r_idx].m_local->second = thread_pool_write[w_idx].m_local->second;
            for ( int l_idx = 0 ; l_idx < thread_pool_write[w_idx].labels_size ; l_idx++ ){
                thread_pool_read[r_idx].labels[l_idx].row_ring_size = thread_pool_write[w_idx].labels[l_idx].row_ring_size;
                thread_pool_read[r_idx].labels[l_idx].g_size = thread_pool_write[w_idx].labels[l_idx].g_size;
                thread_pool_read[r_idx].labels[l_idx].h_size = thread_pool_write[w_idx].labels[l_idx].h_size;
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].g , thread_pool_write[w_idx].labels[l_idx].g , thread_pool_write[w_idx].labels[l_idx].g_size );
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].h , thread_pool_write[w_idx].labels[l_idx].h , thread_pool_write[w_idx].labels[l_idx].h_size );
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].col_ring_size , thread_pool_write[w_idx].labels[l_idx].col_ring_size , thread_pool_write[w_idx].labels[l_idx].row_ring_size );
                strcpy(thread_pool_write[w_idx].labels[l_idx].label, thread_pool_write[w_idx].labels[l_idx].label );
                copyIntMatrix( thread_pool_write[w_idx].labels[l_idx].rings_g , thread_pool_write[w_idx].labels[l_idx].rings_g, 
                                thread_pool_write[w_idx].labels[l_idx].row_ring_size, thread_pool_write[w_idx].labels[l_idx].col_ring_size );
            }
            r_idx++;
        }
        thread_pool_write[w_idx].labels_size = 0;
    }
}






vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                   const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                   std::vector<std::vector<int> >& ring_classes) {
    //vars
    g0 = g00;
    g1 = g11;
    size_t N = DIM_POOL * DIM_POOL;
    ThreadVar *thread_pool_read;
    ThreadVar *thread_pool_write;
    edge_labels =cpu_gen_bond_labels(g0, g1);
    int min_mol_size = std::min(l0.size(), l1.size());
    std::vector <LabelClass> initial_label_classes = cpu_gen_initial_labels(l0, l1, ring_classes);
    GpuLabelClass **gpu_initial_label_classes;
    int gpu_initial_label_classes_size = initial_label_classes.size();
    Pair m_local;
    vector <LabelClass> lcs;
    //cuda Mallocs
    //cuda malloc edge labels
    cudaMallocManaged(&gpu_edge_labels, sizeof(float) * edge_labels.size());
    //cuda malloc adj matrix mol 0
    cudaMallocManaged((void **) &gpu_g0, l0.size() * sizeof(float *));
    for (int i = 0; i < l0.size(); ++i) { cudaMallocManaged((void **) &(gpu_g0[i]), l0.size() * sizeof(float)); }
    //cuda malloc adj matrix mol 1
    cudaMallocManaged((void **) &gpu_g1, l1.size() * sizeof(float *));
    for (int i = 0; i < l1.size(); ++i) { cudaMallocManaged((void **) &(gpu_g1[i]), l1.size() * sizeof(float)); }
    //cuda malloc GpuLabelClass array
    int size = initial_label_classes.size() * 2;
    cudaMallocManaged(&gpu_initial_label_classes, initial_label_classes.size() * sizeof(GpuLabelClass *));
    for (int k = 0; k < initial_label_classes.size(); ++k) {
        cudaMallocManaged(&gpu_initial_label_classes[k], sizeof(GpuLabelClass) * size);
        for (int i = 0; i < size; i++) {
            cudaMallocManaged(&gpu_initial_label_classes[k][i].col_ring_size, sizeof(int *) * min_mol_size);
            cudaMallocManaged(&gpu_initial_label_classes[k][i].g, sizeof(int) * l0.size());
            cudaMallocManaged(&gpu_initial_label_classes[k][i].h, sizeof(int) * l1.size());
            cudaMallocManaged(&gpu_initial_label_classes[k][i].rings_g, sizeof(int *) * min_mol_size);
            for (int j = 0; j < l0.size(); ++j) {
                cudaMallocManaged(&(gpu_initial_label_classes[k][i].rings_g[j]), l0.size() * sizeof(int));}}}
    //cudamalloc / initialize pool
    cudaMallocManaged(&thread_pool_read, sizeof(ThreadVar) * N);
    for (int j = 0; j < N; ++j) {
        thread_pool_read->labels_size = 0;
        thread_pool_read->m_size = 0;
        cudaMallocManaged(&thread_pool_read[j].labels, sizeof(gpu_initial_label_classes));
        cudaMallocManaged(&thread_pool_read[j].m_local, sizeof(Pair) * min_mol_size);
    }
    cudaMallocManaged(&thread_pool_write, sizeof(ThreadVar) * N);
    for (int j = 0; j < N; ++j) {
        thread_pool_write->labels_size = 0;
        thread_pool_write->m_size = 0;
        cudaMallocManaged(&thread_pool_write[j].labels, sizeof(gpu_initial_label_classes));
        cudaMallocManaged(&thread_pool_write[j].m_local, sizeof(Pair) * min_mol_size);
    }

    //initialize
    //init edge labels
    vectorToPointerEdge(gpu_edge_labels);
    size_edge_labels = edge_labels.size();
    //init adj matrix mol0
    vectorToPointerMatrix(g0, gpu_g0);
    size_gpu_g0_row = g0.size();
    size_gpu_g0_col = g0[0].size();
    //init adj matrix mol 1
    vectorToPointerMatrix(g1, gpu_g1);
    size_gpu_g1_row = g1.size();
    size_gpu_g1_col = g1[0].size();

    cout << "CPU: Initializing thread pool" << endl;

    //init n_thread at 1
    int v,w,n_threads=1;
    for( LabelClass lc : initial_label_classes ) {
        v = select_vertex(lc.g,g0);
        w = select_vertex(lc.h,g1);
        if( !matchable(v,w,lc ) ) continue;
        m_local.first = v;
        m_local.second = w;
        lcs = genNewLabels(v,w,initial_label_classes);
        LabelFromCpuToGpu(gpu_initial_label_classes[n_threads],lcs);
        thread_pool_read[n_threads].labels = gpu_initial_label_classes[n_threads];
        thread_pool_read[n_threads].labels_size = lcs.size();
        thread_pool_read[n_threads].m_size = 1;
        thread_pool_read[n_threads].m_local[0] = m_local;
        n_threads++;
    }

    LabelFromCpuToGpu(gpu_initial_label_classes[0],initial_label_classes);
    host_select_label(&thread_pool_write[0].single_label, gpu_initial_label_classes[0],0, initial_label_classes.size() );
    
    std::vector<LabelClass*> label_class_pointers;
    if( !initial_label_classes.empty() ) label_class_pointers.reserve(initial_label_classes.size()+1); // Reserve space for the pointers

    for (LabelClass& item : initial_label_classes) {
        label_class_pointers.push_back(&item); // Add the address of each element to the new vector
    }

    LabelClass* single_label_class_pointer;
    single_label_class_pointer = select_label_DC(label_class_pointers, 0);
    
    cout << "GPU: label selected : " << thread_pool_write[0].single_label.label << endl;
    cout << "CPU : label selected : " << single_label_class_pointer->label << endl;


    int v_tmp = select_vertex(single_label_class_pointer->g,g0);
    int w_tmp = select_vertex(single_label_class_pointer->h,g1);
    int *result_v;
    int *result_w;
    int *result_pos;
    cudaMallocManaged(&result_v, sizeof(int));
    cudaMallocManaged(&result_w, sizeof(int));
    cudaMallocManaged(&result_pos, sizeof(int));
    host_select_vertex(result_v,result_pos, thread_pool_write[0].single_label.g,  thread_pool_write[0].single_label.g_size, gpu_g0, size_gpu_g0_row, size_gpu_g0_col );
    host_select_vertex(result_w,result_pos, thread_pool_write[0].single_label.h,  thread_pool_write[0].single_label.h_size, gpu_g1, size_gpu_g1_row, size_gpu_g1_col );
    

    cout << "GPU: vertex V selected : " << *result_v << endl;
    cout << "CPU: vertex V selected : " << v_tmp << endl;

    cout << "GPU: vertex W selected : " << *result_w << endl;
    cout << "CPU: vertex W selected : " << w_tmp << endl;


    lcs = genNewLabels(v_tmp,w_tmp,initial_label_classes);
    int *idxList;
    cudaMallocManaged(&idxList, sizeof(int) * min_mol_size);
    int new_labels_size = host_gen_new_Labels(thread_pool_write[1].labels, v_tmp, w_tmp, gpu_initial_label_classes[0], initial_label_classes.size(), idxList);

    
    //stampa thread pool read
    /*for ( int j = 0 ; j < n_threads ; ++j ){
        for ( int k = 0 ; k < thread_pool_read[j].labels_size ; ++k ){
            cout<< "\nLABEL : "<<thread_pool_read[j].labels[k].label;
            cout<< "\nADJ : "<<thread_pool_read[j].labels[k].adj;
            cout<<"\n G :  ";
            for ( int row = 0 ; row < thread_pool_read[j].labels[k].g_size ; row++ ){
                cout<<"[ "<< thread_pool_read[j].labels[k].g[row]<<" ]";
            }
            cout<<"\n H :  ";
            for ( int row = 0 ; row < thread_pool_read[j].labels[k].h_size ; row++ ){
                cout<<"[ "<< thread_pool_read[j].labels[k].h[row]<<" ]";
            }
            cout<<"\n ALL RINGS :  ";
            for ( int row = 0 ; row < thread_pool_read[j].labels[k].row_ring_size ; row++ ){
                cout<<"\n{{ ";
                for ( int col = 0 ; col < thread_pool_read[j].labels[k].col_ring_size[row] ; col++ ){
                    cout<<"[ "<< thread_pool_read[j].labels[k].rings_g[row][col] <<" ]";
                }
                cout<<" }} ";
            }
        }
    }*/

    size_t threadsPerBlock;

    size_t numberOfBlocks;

    threadsPerBlock = 8;
    numberOfBlocks = 8;
    int h = 0;
/*

    do{
        printf("\nnew kernel call\n");
    kernel_function<<< threadsPerBlock , numberOfBlocks >>>( thread_pool_read , thread_pool_write);
    cudaDeviceSynchronize();
    
    cout<<"\n\n\nSTAMPA DI THPOOL WRITE POST CHIAMATA A FUNZIONE (stampo lo zero )\n";
    for( int i = 0; i < N ; i++){
        cout<<"[ "<< i <<"] "<<thread_pool_write[i].m_local->first<<endl;
    }
    cpyThreadPool( thread_pool_read , thread_pool_write );

    cout<<"\n\n\nSTAMPA DI THPOOL READ POST COPY (non stampo lo zero)\n";
    for( int i = 0; i < N ; i++){
        if( thread_pool_read[i].m_local->first > 0 )
            cout<<"[ "<<i<<"] "<<thread_pool_read[i].m_local->first<<endl;
    }
    /*int j = 0;
    for( int i = 0; i < N ; i++){ 
        if ( thread_pool_write[i].m_local->first != 0 ){
            thread_pool_read[j].m_local->first = thread_pool_write[i].m_local->first;
            j++;
        }
    }
   


    h++;
    }while( h < 2);


*/











    //cudaFree
    cudaFree( gpu_edge_labels );
    for (int i = 0; i < l0.size(); ++i) {cudaFree(gpu_g0[i]);}
    cudaFree(gpu_g0);
    for (int i = 0; i < l1.size(); ++i) {cudaFree(gpu_g1[i]);}
    cudaFree(gpu_g1);
    for ( int k = 0; k < initial_label_classes.size() ; ++k ) {
        for (int i = 0; i < size; i++) {
            cudaFree(gpu_initial_label_classes[k][i].col_ring_size);
            cudaFree(gpu_initial_label_classes[k][i].g);
            cudaFree(gpu_initial_label_classes[k][i].h);
            for (int j = 0; j < l0.size(); ++j) {
                cudaFree(gpu_initial_label_classes[k][i].rings_g[j]);
            }
            cudaFree(gpu_initial_label_classes[k][i].rings_g);
        }
        cudaFree(gpu_initial_label_classes[k]);
    }
    cudaFree( gpu_initial_label_classes);

    vector<pair<int,int>> m;

    return m;
}