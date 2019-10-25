#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <list>
#include <ctime>
#include <sstream>
#include <vector>
#include <sys/stat.h>
#include <algorithm> 

#include "bwt.h"
#include "SuffixTree.h"
#include "read_genes.h"

using namespace std; 

std::map<char, int> mapOfGenes = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}};
std::map<int, char> mapOfNumbers = {{0, '$'}, {1, 'A'}, {2, 'C'}, {3, 'G'}, {4, 'T'}};

std::map<char, char> map_complement = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}};

// dummy variables used by the comparator
/* string temp_str; */
/* int string_length; */

int* counts_table(const string& str);
int** occurrences_matrix(const string& str); 

// use when we want to create a new index
BWT::BWT(string& rawString, vector< pair<int, int> > holes) {
    this->reference = rawString;
    this->encode();

    this->holes = holes;
    this->C = counts_table(this->encoded);
    this->Occ = occurrences_matrix(this->encoded);
}

// used to create BWT by reading the index from a directory
BWT::BWT(const string& directory){
    recover_index(directory);
}

// make this faster using Ukkonen's algorithm
void BWT::encode() {
    SuffixTree* tree = new SuffixTree(this->reference);
    tree->build();

    this->location_array = new int[this->reference.length()];
    this->encoded = tree->bwt(location_array);
}

// recreates the encoded bwt string based on the occurences matrix
void BWT::recover_encoded_string(int string_length){
    this->encoded.resize(string_length);

    // find first character
    for (int j = 0; j < NUM_BASES; j ++){
        if (this->Occ[0][j] == 1){
            this->encoded[0] = mapOfNumbers[j];
        }
    }

    // fill in rest of the characters
    for(int i = 0; i < string_length; i ++){
        for (int j = 0; j < NUM_BASES; j ++){
            if (this->Occ[i][j] > this->Occ[i][j]){
                this->encoded[0] = mapOfNumbers[j];
            }
        }
    }

    // also recostruct the reference string
    this->reference.resize(string_length);
    for (int i = 0; i < string_length; i ++){
        this->reference[(this->location_array[i] +string_length - 1) % string_length] = this->encoded[i];
    }
}

// This function assumes that location_array is also stored, else use the one below
// void BWT::decode() {
//     unsigned long string_length = this->encoded.length();
//     cout << string_length << endl;
//     this->reference[string_length];
//     cout << "1" << endl;

//     for (unsigned long i = 0; i < string_length; ++i)
//         this->reference[(this->location_array[i] + string_length - 1) % string_length] = this->encoded[i];
// }

void serialize_1d(ostream& outfile, int* arr, int elems) {
    outfile << elems << " ";
    for (int i = 0; i < elems; i++)
        outfile << arr[i] << " ";
}

int* deserialize_1d(istream& file) {
    int elems;
    file >> elems;
    int* arr = new int[elems];
    for (int i = 0; i < elems; i++) {
        file >> arr[i];
    }
    return arr;
}

void serialize_2d(ostream& outfile, int** arr, int rows, int cols) {
    outfile << rows << " ";
    outfile << cols << " ";
    for (int i = 0; i < rows; i++)
        for(int j = 0; j < cols; j++)
            outfile << arr[i][j] << " ";
}

int** deserialize_2d(istream& file, int* string_length) {
    int rows, cols;
    file >> rows;
    file >> cols;
    int** arr = new int*[rows];
    for (int i = 0; i < rows; i++) {
        arr[i] = new int[cols];
        for(int j = 0; j < cols; j++){
            file >> arr[i][j];
        }
    }

    // also change the length of the encoded string
    *string_length = rows;
    return arr;
}

void serialize_vector(ostream& outfile, vector< pair<int, int> > holes){
    outfile << holes.size() << " ";

    for(vector< pair<int, int> >::size_type i = 0; i != holes.size(); i++) {
        outfile << holes[i].first << " ";
        outfile << holes[i].second << " ";
    }
}

vector< pair<int, int> > deserialize_vector(istream& file){
    int size; 
    file >> size;

    vector< pair<int, int> > holes;
    for (int i = 0; i < size; i++){
        int first, second;
        file >> first;
        file >> second;
        holes.push_back(make_pair(first, second));
    }
    
    return holes;
}

// we have to store the location_array, the counts table and the occurences matrix, as well as the holes
// TODO should we also store the encoded? probably no ...
void BWT::store_index(const string& directory){
    
    struct stat info;

    if( stat( directory.c_str(), &info ) != 0 ){
        cout << "Creating directory " << directory << endl;

        const int dir_err = mkdir(directory.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);    
        
        if (dir_err < 0 && dir_err < EEXIST)
        {
            cout << "Error creating directory " << directory << endl;
            exit(1);
        }
    }

    ofstream file_location_array (directory + "/" + LOCATION_ARRAY_FILENAME + FILENAME_TYPE, ios::out);
    serialize_1d(file_location_array, this->location_array, this->reference.length());
    file_location_array.close();

    ofstream file_counts_array (directory + "/" + COUNTS_ARRAY_FILENAME + FILENAME_TYPE, ios::out);
    serialize_1d(file_counts_array, this->C, NUM_BASES);
    file_counts_array.close();

    ofstream file_occurrences_matrix (directory + "/" + OCCURRENCES_MATRIX_FILENAME + FILENAME_TYPE, ios::out);
    serialize_2d(file_occurrences_matrix, this->Occ, this->reference.length(), NUM_BASES);
    file_occurrences_matrix.close();

    ofstream file_holes (directory + "/" + HOLES_FILENAME + FILENAME_TYPE, ios::out);
    serialize_vector(file_holes, this->holes);
    file_holes.close();

    return;
}

void BWT::recover_index(const string& directory){
    struct stat info;

    if( stat( directory.c_str(), &info ) != 0 ){
        cout << "Directory specified: " << directory << " does not exist" << endl;
        exit(1);
    }   

    ifstream file_location_array;
    getFile(directory + "/" + LOCATION_ARRAY_FILENAME + FILENAME_TYPE, file_location_array);
    this->location_array = deserialize_1d(file_location_array);
    file_location_array.close();

    ifstream file_counts_array;
    getFile(directory + "/" + COUNTS_ARRAY_FILENAME + FILENAME_TYPE, file_counts_array);
    this->C = deserialize_1d(file_counts_array);
    file_counts_array.close();

    int* string_length = new int[1];
    ifstream file_occurrences_matrix;
    getFile(directory + "/" + OCCURRENCES_MATRIX_FILENAME + FILENAME_TYPE, file_occurrences_matrix);
    this->Occ = deserialize_2d(file_occurrences_matrix, string_length);
    file_occurrences_matrix.close();

    ifstream file_holes;
    getFile(directory + "/" + HOLES_FILENAME + FILENAME_TYPE, file_holes);
    this->holes = deserialize_vector(file_holes);
    file_holes.close();

    recover_encoded_string(*string_length);
}

// TODO make decode faster by using the location_array !!
/* void BWT::decode(){ */
/*     temp_str = this->encoded; */

/*     int *indices = new int[this->encoded.length()]; */

/*     for (unsigned long i = 0; i < this->encoded.length(); ++i) { */
/*       indices[i] = i; */
/*     } */
    
/*     stable_sort(indices, indices + this->encoded.length(), comparator); */ 

/*     int startIndex = 0; */
/*     for (; this->encoded[startIndex] != '$'; ++startIndex); */

/*     char *builder = new char[this->encoded.length() + 1]; */
/*     for (unsigned long i = 0; i < this->encoded.length(); ++i) { */
/*         startIndex = indices[startIndex]; */
/*         char c = this->encoded[startIndex]; */
/*         builder[i] = c; */
/*     } */
/*     builder[this->encoded.length()] = 0; */

/*     this->reference = string(builder); */
/* } */

string compress(const string& str) {
    vector<char> compressed;

    // previous_char == '*' means that no previous chars need to be stored
    char previous_char = '*';
    unsigned char count = 0;

    for (unsigned long i = 0; i < str.length(); ++i) {
        if (previous_char == '*') {
            previous_char = str[i];
            count = 1;
        } else if (previous_char == str[i]) {
            count ++;
            if (count == 255) {
                compressed.push_back(count);
                compressed.push_back(previous_char);
                previous_char = '*';
                count = 0;
            }
        } else {
            compressed.push_back(count);
            compressed.push_back(previous_char);
            previous_char = str[i];
            count = 1;
        }
    }

    if (count > 0) {
        compressed.push_back(count);
        compressed.push_back(previous_char);
    }

    return string(compressed.begin(), compressed.end());
}

string uncompress(const string str){
    vector<char> uncompressed;

    for (unsigned long int i = 0; i < str.length(); i += 2) {
        for (int j = 0; j < str[i]; j ++)
            uncompressed.push_back(str[i + 1]);
    }

    return string(uncompressed.begin(), uncompressed.end());
}

int* counts_table(const string& str) {
    int* C = new int[NUM_BASES];
    fill_n(C, NUM_BASES, 0);

    for (unsigned long i = 0; i < str.length(); ++i) {
        if (str[i] == 'T')
            continue;

        C[mapOfGenes[str[i]] + 1] ++;
    }
   
    C[2] += C[1];
    C[3] += C[2];
    C[4] += C[3];   

    return C;
}

int** occurrences_matrix(const string& str) {
    int **Occ = new int*[str.length()];
    int current_values[NUM_BASES] = {0};

    for(unsigned long i = 0; i < str.length(); ++i) {
        Occ[i] = new int[NUM_BASES];
        
        current_values[mapOfGenes[str[i]]] += 1;
        for (int j = 0; j < NUM_BASES; j ++) {
            Occ[i][j] = current_values[j];
        }
    }  

    return Occ;
}

block* BWT::get_matches(const string& test) {
    // int* C = counts_table(this->encoded);
    // int** Occ = occurrences_matrix(this->encoded);

    unsigned long i = test.length() - 1;
    int start = C[mapOfGenes[test[i]]];
    int end;

    // for test[i] + 1 assume that all nucleic ashes types have been observed at least once!
    if (test[i] == 'T')
        end = this->encoded.length() - 1;
    else
        end = this->C[mapOfGenes[test[i]] + 1] - 1;

    while (start <= end and i > 0) {
        --i;
        start = this->C[mapOfGenes[test[i]]] + this->Occ[start - 1][mapOfGenes[test[i]]];
        end = this->C[mapOfGenes[test[i]]] + this->Occ[end][mapOfGenes[test[i]]] - 1;
    }

    block* bl = new block;
    bl->start = start;
    bl->end = end;

    //cout << start << ' ' << end << endl;
    return bl;
}

int BWT::get_location(int pos) {
    return this->location_array[pos];
}
