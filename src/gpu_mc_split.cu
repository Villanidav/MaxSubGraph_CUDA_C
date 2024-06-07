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



//struct 
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
    int *idxList;
    Pair *m_local;
}ThreadVar;


//global variables
const int DIM_POOL = 50;
std::vector<std::vector<float>> g00;
std::vector<std::vector<float>> g11;
std::vector<float> edge_labels;

__shared__ float *gpu_edge_labels;
__shared__ int size_edge_labels;


__shared__ float **gpu_g0;
__shared__ int size_gpu_g0_row;
__shared__ int size_gpu_g0_col;

__shared__ float **gpu_g1;
__shared__ int size_gpu_g1_row;
__shared__ int size_gpu_g1_col;


__shared__ Pair *m_best;
__shared__ int m_best_size;



//util functions
void checkError(cudaError_t r) {
  if (r != cudaSuccess) {
    printf("CUDA error on line %d: %s\n", cudaGetErrorString(r));
    exit(0);
  }
}
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
void printLabelClass(GpuLabelClass lb) {
        cout<< lb.label << " [ ";
        cout<< " G("<< lb.g_size << "): ";
        if(!lb.g_size == 0) {for ( int i = 0 ; i < lb.g_size; i++ ) cout<<"["<<lb.g[i]<<"]";}
        cout<< " H("<< lb.h_size << "): ";
        if(!lb.h_size == 0 ) {for ( int i = 0; i <  lb.h_size; i++ ) cout<<"["<<lb.h[i]<<"]";}
        cout<< " RINGS("<< lb.row_ring_size << "): [";
        for( int i = 0; i< lb.row_ring_size; i++){cout<<"("<<lb.col_ring_size[i]<<")"<<"["; for( int j = 0; j <  lb.col_ring_size[i]; j++) cout<<lb.rings_g[i][j]<<", ";  cout<<" ]";}
        cout<<"]";
        cout<< " edge : " <<lb.adj<<" " ;
        cout<< lb.label << " ] "<<endl;
}   
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

}
void cpyGpuLabelClass(GpuLabelClass *l1, GpuLabelClass l2){
    l1->adj = l2.adj;
    l1->row_ring_size = l2.row_ring_size;
    l1->g_size = l2.g_size;
    l1->h_size = l2.h_size;
    strcpy(l1->label , l2.label);
    copyIntArray( l1->g , l2.g, l2.g_size);
    copyIntArray( l1->h, l2.h, l2.h_size);
    copyIntArray( l1->col_ring_size, l2.col_ring_size , l2.row_ring_size);
    copyIntMatrix( l1->rings_g, l2.rings_g, l2.row_ring_size, l2.col_ring_size );/**/
}
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
void vectorToPointerMatrix(const std::vector<std::vector<float>>& g,float** gpu_g) {
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
int cpyThreadPool( ThreadVar *thread_pool_read, ThreadVar *thread_pool_write ){
    int r_idx = 0;
    for( int w_idx = 0; w_idx < DIM_POOL*DIM_POOL ; w_idx++){ 
        thread_pool_read[r_idx].labels_size = thread_pool_write[w_idx].labels_size;
        //da modificare l'if
        if ( thread_pool_write[w_idx].labels_size > 0 ){
            thread_pool_read[r_idx].single_label.g_size = 0;
            thread_pool_read[r_idx].single_label.h_size = 0;
            thread_pool_read[r_idx].single_label.row_ring_size = 0;
            thread_pool_read[r_idx].m_size = thread_pool_write[w_idx].m_size;
            thread_pool_read[r_idx].m_local->first = thread_pool_write[w_idx].m_local->first;
            thread_pool_read[r_idx].m_local->second = thread_pool_write[w_idx].m_local->second;
            for ( int l_idx = 0 ; l_idx < thread_pool_write[w_idx].labels_size ; l_idx++ ){
                thread_pool_read[r_idx].labels[l_idx].adj = thread_pool_write[w_idx].labels[l_idx].adj;
                thread_pool_read[r_idx].labels[l_idx].row_ring_size = thread_pool_write[w_idx].labels[l_idx].row_ring_size;
                thread_pool_read[r_idx].labels[l_idx].g_size = thread_pool_write[w_idx].labels[l_idx].g_size;
                thread_pool_read[r_idx].labels[l_idx].h_size = thread_pool_write[w_idx].labels[l_idx].h_size;
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].g , thread_pool_write[w_idx].labels[l_idx].g , thread_pool_write[w_idx].labels[l_idx].g_size );
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].h , thread_pool_write[w_idx].labels[l_idx].h , thread_pool_write[w_idx].labels[l_idx].h_size );
                copyIntArray( thread_pool_read[r_idx].labels[l_idx].col_ring_size , thread_pool_write[w_idx].labels[l_idx].col_ring_size , thread_pool_write[w_idx].labels[l_idx].row_ring_size );
                for( int c = 0 ; c < 4 ; c++){
                    thread_pool_read[r_idx].labels[l_idx].label[c] = thread_pool_write[w_idx].labels[l_idx].label[c];
                }
                copyIntMatrix( thread_pool_read[r_idx].labels[l_idx].rings_g , thread_pool_write[w_idx].labels[l_idx].rings_g, 
                                thread_pool_write[w_idx].labels[l_idx].row_ring_size, thread_pool_write[w_idx].labels[l_idx].col_ring_size );
            }
            
            for( int ms = 0 ; ms < thread_pool_write[w_idx].m_size ; ms ++){
                thread_pool_read[r_idx].m_local[ms].first = thread_pool_write[w_idx].m_local[ms].first;
                 thread_pool_read[r_idx].m_local[ms].second = thread_pool_write[w_idx].m_local[ms].second;
            }
            thread_pool_read[r_idx].m_size = thread_pool_write[w_idx].m_size;
            r_idx++;
        }
        thread_pool_write[w_idx].labels_size = 0;
    }
    return r_idx;
}



// vtx_set: selected label class
// g: selected graph
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
int host_get_ring_match_data(int *dim_col, int **result, int *idxList, int *elems, int elem_size, GpuLabelClass *lc){
    int index;
    int idx_list_size= 0 ;


    for( int i = 0; i < elem_size ; ++i){
        index = 0;
        for( int j = 0 ; j < lc->g_size ; ++j ){
            if( lc->g[j] == elems[i] ) {idxList[idx_list_size] = index; idx_list_size++;}
            index++;
        }
    }


    for ( int i = 0  ; i < idx_list_size ; ++i  ){
        dim_col[i] = lc->col_ring_size[idxList[i]];
        for ( int j = 0 ; j < dim_col[i] ; ++j ){
            result[i][j] = lc->rings_g[idxList[i]][j];
        }
    }



    return idx_list_size;
}


// return the best select label given an array of labels
void host_select_label(GpuLabelClass *label , GpuLabelClass *lcs, int map_size, int lcs_size){
    int min = 999;
    int max ;
    for( int i = 0 ; i < lcs_size ; ++i ){
        //printf("LABEL CLASSES INTERNE[%d]\n", i);
        //printLabelClass(lcs[i]);
        if( lcs[i].adj == 1 || map_size == 0 ){
            if( lcs[i].g_size > lcs[i].h_size ) max = lcs[i].g_size;
            else max = lcs[i].h_size;
            //printf("\nMAX : %d\n", max);
            if( max < min ){
                
                min = max;
                cpyGpuLabelClass(label, lcs[i] );
            }
        }
    }
    return;
}

// compute the bound given a 1D array of struct GpuLabelClass and its size
int host_calc_bound(GpuLabelClass *lcs, int lc_size) {
    int bound = 0;
    for( int i = 0 ; i < lc_size ; ++i){
        if ( lcs[i].g_size > lcs[i].h_size ) bound = bound + lcs[i].h_size;
        else bound = bound + lcs[i].g_size;
    }
    return bound;
}


//return = size of the friends
//friend is the OUTPUT
int host_hoodG(int *friends,int vtx, float edge, float **g, int size_g) {
    int size = 0;
    
    for (int i = 0; i < size_g; i++) {
        if ( g[i][vtx] == edge && vtx != i) {
            friends[size] = i;
            size++;
        }
    }

    return size;
}


// result == size of generated label
// output is l_draft
// input : v
void host_resize(int *array, int size_arr, int place_availabel){
    int count = 0;
    /*bool flag = true;
    for(int i = 0; i < size_arr && flag ; ++i){
        if(array[i] == -1){continue;}
        if( i == count){ count ++; continue; }
        array[count] = array[i]; count++;
        if(count == place_availabel){flag = false;}
    }*/


    for( int i = 0 ; i < size_arr && place_availabel > 0 ; i++){
        if( array[i] != -1 ){
            place_availabel--;
            array[count] = array[i];
            count ++;
        }
    }
}


int host_gen_new_labels(GpuLabelClass *l_draft ,  int v, int w, GpuLabelClass *lcs, int lcs_size, int *idxList) {
    int vs,ws, draft_size = 0;
    int dim_row;
    int count = 0;
    for ( int i = 0 ; i < lcs_size ; ++i ){
       
        for ( int j = 0 ; j < size_edge_labels ; ++j ){
            int friendsize;
            friendsize = host_hoodG( l_draft[draft_size].g , v , gpu_edge_labels[j], gpu_g0 , size_gpu_g0_row);
        
            vs = 0;
            for ( int k = 0; k < friendsize ; ++k ){
                if( host_contains(l_draft[draft_size].g[k] , lcs[i].g, lcs[i].g_size) ){ vs++;  }
                else{ l_draft[draft_size].g[k] = -1;}
            }
        
            host_resize(l_draft[draft_size].g, friendsize, vs );

            dim_row = host_get_ring_match_data(l_draft[draft_size].col_ring_size, l_draft[draft_size].rings_g, idxList ,l_draft[draft_size].g, vs, &lcs[i] );
            

            friendsize = host_hoodG(l_draft[draft_size].h, w, gpu_edge_labels[j], gpu_g1, size_gpu_g1_row );
            //printf("\n esco da hood 2");
            ws = 0;
            for ( int k = 0 ; k < friendsize ; ++k ){
                
                if( host_contains(l_draft[draft_size].h[k], lcs[i].h, lcs[i].h_size) ){  ws++; }
                else {
                    l_draft[draft_size].h[k] = -1;
                }
            }
            host_resize(l_draft[draft_size].h, friendsize, ws );
    
            int adj;
            if ( ws > 0 && vs > 0 ){
                if( gpu_edge_labels[j] != 0.0 || lcs[i].adj == 1 ) {adj = 1;}
                else { adj = 0; }

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


//given two atoms from the same label, return true if they are matchable, false otherwise
// based on how their rings matches
bool host_matchable(int **v_ring_atoms, int v, int w, GpuLabelClass *lc, int *idxList) {

    host_get_ring_match_data(lc->col_ring_size, v_ring_atoms, idxList , &v, 1 ,lc);
    if( lc->col_ring_size[idxList[0]] > 0  ){
        for(int i = 0; i < lc->col_ring_size[idxList[0]]  ; i++){
            if( v_ring_atoms[0][i] == -1 )return false;
            if( v_ring_atoms[0][i] == w ) return true;
        }
        return false;
    }
    return true;
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



//function used in the main -- equivalent to the ones above but using LabelClasses objects, not struct
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

            for(int vtx : hood(v,g00,edge_l)){
                if( std::find(label.g.begin(),label.g.end(),vtx) != label.g.end() ){
                    v_conn.push_back(vtx);
                }
            }
            v_c_rings = label.get_ring_match_data(v_conn);

            for(int vtx : hood(w,g11,edge_l)){
                if(std::find(label.h.begin(),label.h.end(),vtx) != label.h.end() ){
                    w_conn.push_back(vtx);
                }
            }

            int adj;
            if(!v_conn.empty() && !w_conn.empty()){
                if(edge_l != 0.0 || label.adj == 1){
                    adj = 1;
                }else{
                    adj = 0;
                }
                LabelClass tmp(v_conn,w_conn,v_c_rings,adj, label.label);
                l_draft.push_back(tmp);
            }
        }
    }
    return l_draft;
}




void host_parallel_solve_mcs( ThreadVar *thread_pool_read, ThreadVar *thread_pool_write, int n_threads ){
    
    int index ;
    int z = 0;
    int space = 4;
    int flag = 0;
    
    for( index = 0 ; index < n_threads ; index++ ){
        GpuLabelClass *label = &thread_pool_read[index].single_label;
        host_select_label(label, thread_pool_read[index].labels, 
                            thread_pool_read[index].m_size, thread_pool_read[index].labels_size);
        
        flag = 0;
        if ( (thread_pool_read[index].m_size + host_calc_bound(thread_pool_read[index].labels, thread_pool_read[index].labels_size  )  
            < m_best_size) || !label  ) { flag = 1; }
        if( flag == 0 ){
            int jump = 0;
            for( int v_idx = 0 ; v_idx < label->g_size ; ++v_idx){
                for( int w_idx = 0 ; w_idx < label->h_size ; ++w_idx){

                    if( !host_matchable(label->rings_g, label->g[v_idx], label->h[w_idx], label, thread_pool_read[index].idxList ) ) continue;
                    
                    for( z = 0 ; z < thread_pool_read[index].m_size ; z++ ){
                        thread_pool_write[4*index + jump].m_local[z].first = thread_pool_read[index].m_local[z].first;
                        thread_pool_write[4*index + jump].m_local[z].second = thread_pool_read[index].m_local[z].second;
                    }
                    thread_pool_write[4*index + jump].m_size = thread_pool_read[index].m_size + 1;
                    thread_pool_write[4*index + jump].m_local[z].first = label->g[v_idx];
                    thread_pool_write[4*index + jump].m_local[z].second = label->h[v_idx];
                    

                    int l_s = host_gen_new_labels( thread_pool_write[(4*index + jump)].labels , label->g[v_idx], label->h[w_idx] ,
                                                                                            thread_pool_read[index].labels , thread_pool_read[index].labels_size ,
                                                                                            thread_pool_read[index].idxList);
                    
                    thread_pool_write[(4*index + jump)].labels_size = l_s;
                    if( thread_pool_write[4*index + jump].m_size > m_best_size ){
                        m_best_size = thread_pool_write[4*index + jump].m_size ;
                        for( int z = 0 ; z < thread_pool_write[4*index + jump].m_size ; z++ ){
                            m_best[z].first = thread_pool_write[4*index + jump].m_local[z].first;
                            m_best[z].second = thread_pool_write[4*index + jump].m_local[z].second;
                        }

                    }
                jump ++;
                }
            }
       }
    }
}



vector<pair<int,int>> gpu_mc_split(const std::vector<std::vector<float>>& g000, const std::vector<std::vector<float>>& g111,
                                   const std::vector<std::string>& l0, const std::vector<std::string>& l1,
                                   std::vector<std::vector<int> >& ring_classes) {
    g00 = g000;
    g11 = g111;
    size_t N = DIM_POOL * DIM_POOL;
    edge_labels = cpu_gen_bond_labels(g00, g11);
    std::vector <LabelClass> initial_label_classes = cpu_gen_initial_labels(l0, l1, ring_classes);

    
 
    //cuda Mallocs
    cudaMallocManaged(&m_best , sizeof(Pair)* l1.size());
    //cuda malloc edge labels
    checkError(cudaMallocManaged(&gpu_edge_labels, sizeof(float) * edge_labels.size()));
    //cuda malloc adj matrix mol 0
    checkError(cudaMallocManaged((void **) &gpu_g0, l0.size() * sizeof(float *)));
    for (int i = 0; i < l0.size(); ++i) { checkError(cudaMallocManaged((void **) &(gpu_g0[i]), l0.size() * sizeof(float))); }
    //cuda malloc adj matrix mol 1
    checkError(cudaMallocManaged((void **) &gpu_g1, l1.size() * sizeof(float *)));
    for (int i = 0; i < l1.size(); ++i) { checkError(cudaMallocManaged((void **) &(gpu_g1[i]), l1.size() * sizeof(float))); }
    //cuda malloc GpuLabelClass array
    int size = initial_label_classes.size() * 3;
    printf("INITIAL SIZE %d", size);
    int min_mol_size = std::min(l0.size(), l1.size());

    //cudamalloc / initialize pool
    ThreadVar *thread_pool_read;
    ThreadVar *thread_pool_write;

    checkError(cudaMallocManaged(&thread_pool_read, sizeof(ThreadVar) * N/5));
    for (int j = 0; j < N/5; ++j) {

        checkError(cudaMallocManaged(&thread_pool_read[j].single_label.col_ring_size, sizeof(int ) * min_mol_size));
        checkError(cudaMallocManaged(&thread_pool_read[j].single_label.g, sizeof(int) * l0.size()));
        checkError(cudaMallocManaged(&thread_pool_read[j].single_label.h, sizeof(int) * l1.size()));
        checkError(cudaMallocManaged(&thread_pool_read[j].single_label.rings_g, sizeof(int *) * l0.size()));
        for (int h = 0; h < l0.size(); ++h) {
            checkError(cudaMallocManaged(&(thread_pool_read[j].single_label.rings_g[h]), l0.size() * sizeof(int)));
            }
        checkError(cudaMallocManaged(&thread_pool_read[j].m_local, sizeof(Pair) * min_mol_size));
        checkError(cudaMallocManaged(&thread_pool_read[j].idxList, sizeof(int) * min_mol_size));

        checkError(cudaMallocManaged(&thread_pool_read[j].labels, size * sizeof(GpuLabelClass)));
        for (int k = 0; k < size; k++) {
                checkError(cudaMallocManaged(&thread_pool_read[j].labels[k].col_ring_size, sizeof(int ) * min_mol_size));
                checkError(cudaMallocManaged(&thread_pool_read[j].labels[k].g, sizeof(int) * l0.size()));
                checkError(cudaMallocManaged(&thread_pool_read[j].labels[k].h, sizeof(int) * l1.size()));
                checkError(cudaMallocManaged(&thread_pool_read[j].labels[k].rings_g, sizeof(int *) * l0.size()));
                for (int h = 0; h < l0.size(); ++h) {
                    checkError(cudaMallocManaged(&(thread_pool_read[j].labels[k].rings_g[h]), l0.size() * sizeof(int)));
                }
        }
    
    }
    checkError(cudaMallocManaged(&thread_pool_write, sizeof(ThreadVar) * N));
    for (int j = 0; j < N; ++j) {
        checkError(cudaMallocManaged(&thread_pool_write[j].labels, size * sizeof(GpuLabelClass)));
        
        for (int k = 0; k < size; k++) {
                checkError(cudaMallocManaged(&thread_pool_write[j].labels[k].col_ring_size, sizeof(int ) * min_mol_size));
                checkError(cudaMallocManaged(&thread_pool_write[j].labels[k].g, sizeof(int) * l0.size()));
                checkError(cudaMallocManaged(&thread_pool_write[j].labels[k].h, sizeof(int) * l1.size()));
                checkError(cudaMallocManaged(&thread_pool_write[j].labels[k].rings_g, sizeof(int *) * l0.size()));
                for (int h = 0; h < l0.size(); ++h) {
                    checkError(cudaMallocManaged(&(thread_pool_write[j].labels[k].rings_g[h]), l0.size() * sizeof(int)));
                }
        }
        checkError(cudaMallocManaged(&thread_pool_write[j].single_label.col_ring_size, sizeof(int ) * min_mol_size));
        checkError(cudaMallocManaged(&thread_pool_write[j].single_label.g, sizeof(int) * l0.size()));
        checkError(cudaMallocManaged(&thread_pool_write[j].single_label.h, sizeof(int) * l1.size()));
        checkError(cudaMallocManaged(&thread_pool_write[j].single_label.rings_g, sizeof(int *) * l0.size()));
        for (int h = 0; h < l0.size(); ++h) {
            checkError(cudaMallocManaged(&(thread_pool_write[j].single_label.rings_g[h]), l0.size() * sizeof(int)));
            }
        checkError(cudaMallocManaged(&thread_pool_write[j].m_local, sizeof(Pair) * min_mol_size));
        checkError(cudaMallocManaged(&thread_pool_write[j].idxList, sizeof(int) * min_mol_size));
    
    }

    //initialize
    //init edge labels
    vectorToPointerEdge(gpu_edge_labels);
    size_edge_labels = edge_labels.size();
    //init adj matrix mol0
    vectorToPointerMatrix(g00, gpu_g0);
    size_gpu_g0_row = g00.size();
    size_gpu_g0_col = g00[0].size();
    //init adj matrix mol 1
    vectorToPointerMatrix(g11, gpu_g1);
    size_gpu_g1_row = g11.size();
    size_gpu_g1_col = g11[0].size();
    Pair m_local;
    vector <LabelClass> lcs;
    cout << "CPU: Initializing thread pool" << endl;

    //init n_thread at 1
    int size_lcs = 0;
    int v,w,n_threads=0;
    for( LabelClass lc : initial_label_classes ) {
        v = select_vertex(lc.g,g00);
        w = select_vertex(lc.h,g11);
        if( !matchable(v,w,lc ) ) continue;
        m_local.first = v;
        m_local.second = w;
        lcs = genNewLabels(v,w,initial_label_classes);
        size_lcs = lcs.size();
        LabelFromCpuToGpu(thread_pool_read[n_threads].labels,lcs);
        thread_pool_read[n_threads].labels_size = lcs.size();
        thread_pool_read[n_threads].m_size = 1;
        thread_pool_read[n_threads].m_local[0] = m_local;
        n_threads++;
    }

    int max = n_threads;
    int level = 0;
    do{
        host_parallel_solve_mcs(thread_pool_read, thread_pool_write , n_threads);
        n_threads = cpyThreadPool( thread_pool_read , thread_pool_write );
        level++;
    }while( n_threads > 0 );

    vector<pair<int,int>> m;
    pair<int, int> tmp;
    for( int best = 0 ; best < m_best_size ; best++ ){
        tmp.first = m_best[best].first;
        tmp.second = m_best[best].second;
        m.push_back(tmp);
    }

    //cudaFree
    cudaFree( gpu_edge_labels );
    for (int i = 0; i < l0.size(); ++i) {cudaFree(gpu_g0[i]);}
    cudaFree(gpu_g0);
    for (int i = 0; i < l1.size(); ++i) {cudaFree(gpu_g1[i]);}
    cudaFree(gpu_g1);
    
    for (int j = 0; j < N; ++j) {
        for (int k = 0; k < size; k++) {
                checkError(cudaFree(thread_pool_write[j].labels[k].col_ring_size));
                checkError(cudaFree(thread_pool_write[j].labels[k].g));
                checkError(cudaFree(thread_pool_write[j].labels[k].h));
                for (int h = 0; h < l0.size(); ++h) {
                    checkError(cudaFree((thread_pool_write[j].labels[k].rings_g[h])));
                }
                checkError(cudaFree(thread_pool_write[j].labels[k].rings_g));
        }
        checkError(cudaFree(thread_pool_write[j].labels));

        checkError(cudaFree(thread_pool_write[j].single_label.col_ring_size));
        checkError(cudaFree(thread_pool_write[j].single_label.g));
        checkError(cudaFree(thread_pool_write[j].single_label.h));
        for (int h = 0; h < l0.size(); ++h) {
            checkError(cudaFree((thread_pool_write[j].single_label.rings_g[h])));
            }
        checkError(cudaFree(thread_pool_write[j].single_label.rings_g));
        checkError(cudaFree(thread_pool_write[j].m_local));
        checkError(cudaFree(thread_pool_write[j].idxList));
    }
    checkError(cudaFree(thread_pool_write));

    for (int j = 0; j < N/5; ++j) {
        for (int k = 0; k < size; k++) {
                checkError(cudaFree(thread_pool_read[j].labels[k].col_ring_size));
                checkError(cudaFree(thread_pool_read[j].labels[k].g));
                checkError(cudaFree(thread_pool_read[j].labels[k].h));
                for (int h = 0; h < l0.size(); ++h) {
                    checkError(cudaFree((thread_pool_read[j].labels[k].rings_g[h])));
                }
                checkError(cudaFree(thread_pool_read[j].labels[k].rings_g));
        }
        checkError(cudaFree(thread_pool_read[j].labels));

        checkError(cudaFree(thread_pool_read[j].single_label.col_ring_size));
        checkError(cudaFree(thread_pool_read[j].single_label.g));
        checkError(cudaFree(thread_pool_read[j].single_label.h));
        for (int h = 0; h < l0.size(); ++h) {
            checkError(cudaFree((thread_pool_read[j].single_label.rings_g[h])));
            }
        checkError(cudaFree(thread_pool_read[j].single_label.rings_g));
        checkError(cudaFree(thread_pool_read[j].m_local));
        checkError(cudaFree(thread_pool_read[j].idxList));
    }
    checkError(cudaFree(thread_pool_read));

    return m;
}

