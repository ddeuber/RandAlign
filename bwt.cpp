#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <list>
#include <ctime>
#include <sstream>
#include <vector>
#include <algorithm> 

#include "bwt.h"
#include "SuffixTree.h"

using namespace std; 

std::map<char, int> mapOfGenes = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}};

std::map<char, char> map_complement = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}};

// dummy variables used by the comparator
/* string temp_str; */
/* int string_length; */

int* counts_table(const string& str);
int** occurrences_matrix(const string& str); 

BWT::BWT(string& rawString, bool is_encoded) {
    if (is_encoded){
        this->encoded = rawString;
        this->decode();
    } else {
        this->reference = rawString;
        this->encode();

        this->C = counts_table(this->encoded);
        this->Occ = occurrences_matrix(this->encoded);
    }   
}

// make this faster using Ukkonen's algorithm
void BWT::encode() {
    SuffixTree* tree = new SuffixTree(this->reference);
    tree->build();

    this->location_array = new int[this->reference.length()];
    this->encoded = tree->bwt(location_array);
}

/* bool comparator(int a, int b) { */
/*     return temp_str[a] < temp_str[b]; */
/* } */

// This function assumes that location_array is also stored, else use the one below
void BWT::decode() {
    unsigned long string_length = this->encoded.length();
    cout << string_length << endl;
    this->reference[string_length];
    cout << "1" << endl;

    for (unsigned long i = 0; i < string_length; ++i)
        this->reference[(this->location_array[i] + string_length - 1) % string_length] = this->encoded[i];
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

    cout << start << ' ' << end << endl;
    return bl;
}

int BWT::get_location(int pos) {
    return this->location_array[pos];
}