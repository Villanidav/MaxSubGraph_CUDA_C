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
const int DIM_POOL = 16;
std::vector<std::vector<float>> g0;
std::vector<std::vector<float>> g1;
std::vector<float> edge_labels;
vector<pair<int,int>> m_best;

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
    Pair *m_local;
}ThreadVar;

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

    for( int i = 0; i < elem_size ; ++i){
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) idxList[index] = index;
            index++;
        }
    }
    for ( int k = 0 , i = 0 ; k < index+1 ; ++k , ++i ){
        result[i] = lc->rings_g[idxList[k]];
        dim_col[i] = lc->col_ring_size[idxList[k]];
    }
    return index+1;
}

int host_get_ring_match_data(int *dim_col, int **result, int *idxList, int *elems, int elem_size, GpuLabelClass *lc){
    int index;

    for( int i = 0; i < elem_size ; ++i){
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) idxList[index] = index;
            index++;
        }
    }
    for ( int k = 0 , i = 0 ; k < index+1 ; ++k , ++i ){
        result[i] = lc->rings_g[idxList[k]];
        dim_col[i] = lc->col_ring_size[idxList[k]];
    }
    return index+1;
}
// function used in matchable, returns a 1D array containing the ring data related to vertex v
// returns the size of the 1D array
__device__
int kernel_matchable_ring_data(int *result, int v, GpuLabelClass *lc){
    int index=0;
    for ( int i = 0 ; i < lc->g_size ; ++i){
        if ( lc->g[i] == v ) result[index] = index;
        index++;
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

// compute the bound given a 1D array of struct GpuLabelClass and its size
__device__
int calc_bound(GpuLabelClass *lcs, int lc_size) {
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
int kernel_gen_new_Labels(GpuLabelClass *l_draft ,  int v, int w, GpuLabelClass *lcs, int lcs_size, int *v_conn, int *w_conn, int **v_c_rings, int *friends, int *idxList, int *dim_col) {
    int vs,ws, draft_size = 0;
    int dim_row;
    for ( int i = 0 ; i < lcs_size ; ++i ){
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize = hoodG(friends, v, gpu_edge_labels[j], gpu_g0 );
            for ( int k = 0 , vs = 0 ; k < friendsize ; ++k ){
                if( contains(friends[k], lcs[i].g, lcs[i].g_size) ){ v_conn[vs] = friends[k]; vs++; }
            }

            dim_row = kernel_get_ring_match_data(dim_col, v_c_rings, idxList ,v_conn, vs+1, &lcs[i] );

            friendsize = hoodG(friends, w, gpu_edge_labels[j], gpu_g1 );
            for ( int k = 0 , ws = 0 ; k < friendsize ; ++k ){
                if( contains(friends[k], lcs[i].h, lcs[i].h_size) ){ w_conn[ws] = friends[k]; ws++; }
            }

            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) adj = 1;
                else adj = 0;
                l_draft[draft_size].g_size = vs;
                l_draft[draft_size].h_size = ws;
                l_draft[draft_size].row_ring_size = dim_row;
                l_draft[draft_size].col_ring_size = dim_col;
                l_draft[draft_size].g = v_conn;
                l_draft[draft_size].h = w_conn;
                l_draft[draft_size].adj = adj;
                for( int c = 0 ; c < 4 ; c++){
                    l_draft[draft_size].label[c] = lcs[i].label[c];
                }
                l_draft[draft_size].rings_g = v_c_rings;
                draft_size++;
            }
        }
    }
    return draft_size;
}


int host_gen_new_Labels(GpuLabelClass *l_draft,int v, int w, GpuLabelClass *lcs, int lcs_size, int *v_conn, int *w_conn, int **v_c_rings, int *friends, int *idxList, int *dim_col) {
    int vs,ws, draft_size = 0;
    int dim_row;
    for ( int i = 0 ; i < lcs_size ; ++i ){
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize = host_hoodG(friends, v, gpu_edge_labels[j], gpu_g0 );
            for ( int k = 0 , vs = 0 ; k < friendsize ; ++k ){
                if( host_contains(friends[k], lcs[i].g, lcs[i].g_size) ){ v_conn[vs] = friends[k]; vs++; }
            }

            dim_row = host_get_ring_match_data(dim_col, v_c_rings, idxList ,v_conn, vs+1, &lcs[i] );

            friendsize = host_hoodG(friends, w, gpu_edge_labels[j], gpu_g1 );
            for ( int k = 0 , ws = 0 ; k < friendsize ; ++k ){
                if( host_contains(friends[k], lcs[i].h, lcs[i].h_size) ){ w_conn[ws] = friends[k]; ws++; }
            }

            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) adj = 1;
                else adj = 0;
                l_draft[draft_size].g_size = vs;
                l_draft[draft_size].h_size = ws;
                l_draft[draft_size].row_ring_size = dim_row;
                l_draft[draft_size].col_ring_size = dim_col;
                l_draft[draft_size].g = v_conn;
                l_draft[draft_size].h = w_conn;
                l_draft[draft_size].adj = adj;
                for( int c = 0 ; c < 4 ; c++){
                    l_draft[draft_size].label[c] = lcs[i].label[c];
                }
                l_draft[draft_size].rings_g = v_c_rings;
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


__device__ bool kernel_matchable(int*  v_ring_atoms, int v, int w, GpuLabelClass *lc) {
    int size  = kernel_matchable_ring_data(v_ring_atoms, v, lc );
    if( size > 0  ) {
        for(int i = 0; i< size; i++){
            if( v_ring_atoms[i] == -1 )return false;
            if( v_ring_atoms[i] == w ) return true;
        }
        return false;
    }
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
int host_select_vertex(int result, int *result_pos, int *vtx_set, int vtx_size, float **g, int num_row, int num_column) {
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
            result = vtx_set[i];
            *result_pos = i;
        }
    }
    return result;
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




bool solve_mcs() {

    /*queue_elem elem =  Q.back();

    Q.pop_back();

    vector<LabelClass> lcs = elem.labels;
    vector<pair<int,int >> m_local = elem.m_local;

    std::vector<LabelClass*> label_class_pointers;
    label_class_pointers.reserve(lcs.size());
    for (LabelClass& item : lcs) {label_class_pointers.push_back(&item);}


    LabelClass *lcc = select_label(label_class_pointers, m_local.size());
    if ( m_local.size() + calc_bound(lcs) <= m_best.size() || ( !lcc && !m_local.empty() )  ){ if( !Q.empty() ){ return true; } return false;}

    queue_elem qel;
    LabelClass lc = *lcc;
    pair<int,int> m_temp;

    for( int v : lc.g )  {
        for ( int w : lc.h ) {
            if ( !matchable(v,w,lc) ) continue;
            m_temp.first = v;
            m_temp.second = w;
            m_local.push_back(m_temp);
            qel.labels = genNewLabels(v,w,lcs);
            qel.m_local = m_local;
            Q.push_back(qel);
            if ( m_local.size() > m_best.size() ) m_best = m_local;
            m_local.pop_back();
        }
    }

    return true;*/
}






vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g00, const std::vector<std::vector<float>>& g11,
                                   const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                   std::vector<std::vector<int> >& ring_classes) {
    cout << "CPU: Initializing ..." << endl;
    //vars
    g0 = g00;
    g1 = g11;
    ThreadVar *thread_pool;
    size_t N = DIM_POOL * DIM_POOL;
    edge_labels = cpu_gen_bond_labels(g0, g1);
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
    cudaMallocManaged(&thread_pool, sizeof(ThreadVar) * N);
    for (int j = 0; j < N; ++j) {
        thread_pool->labels_size = 0;
        thread_pool->m_size = 0;
        cudaMallocManaged(&thread_pool[j].labels, sizeof(gpu_initial_label_classes));
        cudaMallocManaged(&thread_pool[j].m_local, sizeof(Pair) * min_mol_size);
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


    int v,w,n_threads=0;
    for( LabelClass lc : initial_label_classes ) {
        v = select_vertex(lc.g,g0);
        w = select_vertex(lc.h,g1);
        if( !matchable(v,w,lc ) ) continue;
        m_local.first = v;
        m_local.second = w;
        lcs = genNewLabels(v,w,initial_label_classes);
        LabelFromCpuToGpu(gpu_initial_label_classes[n_threads],lcs);
        thread_pool[n_threads].labels = gpu_initial_label_classes[n_threads];
        thread_pool[n_threads].labels_size = lcs.size();
        thread_pool[n_threads].m_size = 1;
        thread_pool[n_threads].m_local[0] = m_local;
        n_threads++;
    }



    for ( int j = 0 ; j < n_threads ; ++j ){
        for ( int k = 0 ; k < thread_pool[j].labels_size ; ++k ){
            cout<< "\nLABEL : "<<thread_pool[j].labels[k].label;
            cout<< "\nADJ : "<<thread_pool[j].labels[k].adj;
            cout<<"\n G :  ";
            for ( int row = 0 ; row < thread_pool[j].labels[k].g_size ; row++ ){
                cout<<"[ "<< thread_pool[j].labels[k].g[row]<<" ]";
            }
            cout<<"\n H :  ";
            for ( int row = 0 ; row < thread_pool[j].labels[k].h_size ; row++ ){
                cout<<"[ "<< thread_pool[j].labels[k].h[row]<<" ]";
            }
            cout<<"\n ALL RINGS :  ";
            for ( int row = 0 ; row < thread_pool[j].labels[k].row_ring_size ; row++ ){
                cout<<"\n{{ ";
                for ( int col = 0 ; col < thread_pool[j].labels[k].col_ring_size[row] ; col++ ){
                    cout<<"[ "<< thread_pool[j].labels[k].rings_g[row][col] <<" ]";
                }
                cout<<" }} ";
            }
        }
    }



    /*bool flag;
    do{flag = solve_mcs();}while(flag);*/




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


    return m_best;
}