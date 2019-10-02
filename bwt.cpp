#include <iostream>
#include <fstream>
#include <string>
#include <assert.h>
#include <list>
#include <ctime>
#include <bits/stdc++.h>
#include <map>
#include <algorithm> 

#include "bwt.hpp"

#define SEED_LENGTH 12 
#define NUM_BASES 5
#define MISMATCHES_THRESHOLD 5

using namespace std; 

// dummy variables used by the comparator
string temp_str;
int string_length;

// map from symbols to ints
map<char, int> mapOfGenes = {{'$', 0}, {'A', 1}, {'C', 2}, {'G', 3}, {'T', 4}};
map<char, char> map_complement = {{'A', 'T'}, {'C', 'G'}, {'G', 'C'}, {'T', 'A'}};

// used to return the indices for the matching lines found
typedef struct {
    int start;
    int end;
} block;

typedef struct {
    string id;
    string forward_read;
    string forward_read_quality;
    string backward_read;
    string backward_read_quality;
} read_block;

bool indices_comparator(int a, int b) {
    int pos = 0;

    while (temp_str[(a + pos) % string_length] == temp_str[(b + pos) % string_length]){
        pos ++;
    }

    return temp_str[(a + pos) % string_length] < temp_str[(b + pos) % string_length];
}

// TODO make this faster using Ukkonen's algorithm
void encode(string rawString, string& encoded, int* location_array) {
    temp_str = rawString;

    string_length = rawString.length();

    for (int i = 0; i < string_length; ++i) {
            location_array[i] = i;
    }

    sort(location_array, location_array + string_length, indices_comparator); 
    
    char *encoded_ptr = new char[string_length + 1];
    for (unsigned long i = 0; i < rawString.length(); ++i) {
        encoded_ptr[i] = rawString[(location_array[i] + string_length - 1) % string_length];
    }    
    encoded_ptr[string_length] = 0;

    encoded = string(encoded_ptr);
}

bool comparator(int a, int b) {
    return temp_str[a] < temp_str[b];
}

// TODO make decode faster by using the location_array !!
string decode(string encoded){
    temp_str = encoded;

    int *indices = new int[encoded.length()];

    for (unsigned long i = 0; i < encoded.length(); ++i) {
      indices[i] = i;
    }
    
    stable_sort(indices, indices + encoded.length(), comparator); 

    int startIndex = 0;
    for (; encoded[startIndex] != '$'; ++startIndex);

    char *builder = new char[encoded.length() + 1];
    for (unsigned long i = 0; i < encoded.length(); ++i) {
        startIndex = indices[startIndex];
        char c = encoded[startIndex];
        builder[i] = c;
    }
    builder[encoded.length()] = 0;

    return string(builder);
}

string compress(string str) {
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

string uncompress(string str){
    vector<char> uncompressed;

    for (unsigned long int i = 0; i < str.length(); i += 2) {
        for (int j = 0; j < str[i]; j ++)
            uncompressed.push_back(str[i + 1]);
    }

    return string(uncompressed.begin(), uncompressed.end());
}

string join(const vector<string>& vec, const char* delim)
{
    stringstream res;
    copy(vec.begin(), vec.end(), ostream_iterator<string>(res, delim));
    return res.str();
}

int* counts_table(string str) {
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

int** occurrences_matrix(string str) {
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

block* get_matches(string ref, string test) {
    int* C = counts_table(ref);
    int** Occ = occurrences_matrix(ref);

    unsigned long i = test.length() - 1;
    int start = C[mapOfGenes[test[i]]];
    int end;

    // for test[i] + 1 assume that all nucleic ashes types have been observed at least once!
    if (test[i] == 'T')
        end = ref.length() - 1;
    else
        end = C[mapOfGenes[test[i]] + 1] - 1;

    while (start <= end and i > 0) {
        --i;
        start = C[mapOfGenes[test[i]]] + Occ[start - 1][mapOfGenes[test[i]]];
        end = C[mapOfGenes[test[i]]] + Occ[end][mapOfGenes[test[i]]] - 1;
    }

    block* bl = new block;
    bl->start = start;
    bl->end = end;

    return bl;
    
}

string read_reference_gene(string file){
    string gene;
    vector<string> vector_genes;

    ifstream read_file(file);

    if (read_file){
        getline (read_file, gene);

        while(getline(read_file, gene)) {
            vector_genes.push_back(gene);
        }
        vector_genes.push_back("$");

        return join(vector_genes, "");
    } else {
        cout << "File: " << file << "cannot be opened" << endl; 
        exit(1);
    }
}

string extract_id(string s) {
    return s.substr(1, s.find("/") - 1);
}

read_block* get_one_read(ifstream& file_1, ifstream& file_2){
    read_block* rb = new read_block;

    try {
        getline(file_1, rb->forward_read);
        rb->id = extract_id(rb->forward_read);

        getline(file_1, rb->forward_read);
        getline(file_1, rb->forward_read_quality);
        getline(file_1, rb->forward_read_quality);

        getline(file_2, rb->backward_read);
        assert(!extract_id(rb->backward_read).compare(rb->id));
        getline(file_2, rb->backward_read);
        getline(file_2, rb->backward_read_quality);
        getline(file_2, rb->backward_read_quality);
    } catch(std::out_of_range e){
        return NULL;
    }

    return rb;
}

bool checkFileExistence(const string& filename)
{
    ifstream f(filename.c_str());
    return f.is_open();
}
    
void getFile(string filename, /*out*/ ifstream& file)
{
    const bool file_exists = checkFileExistence(filename);
    if (!file_exists) {
        cout << "File " << filename << " not found." << endl;
        exit(1);
    }
    file.open(filename.c_str());
}

int number_of_mismatches(string s1, string s2) {

    // assume len(s1) == len(s2) 
    int mismatches = 0;
    for (unsigned long i = 0; i < s1.length(); ++i)
        if (s1[i] != s2[i])
            mismatches ++;

    return mismatches;
}

// TODO modify this method to add to file using SAM format!
void print_matches(block* bl, string id, string test, int* location_array, string reference) {
    if (bl->start > bl->end)
        cout << "No matches found" << endl;
    else {
        // assume that from all matches the one corresponding to the least amount of mismatches is the allignment 

        // foreach match, here only take first
        int pos = location_array[bl->start];
        int mismatches;
        
        // chack that the position of seed found is not too close to the end of the string
        if (pos + test.length() > reference.length())
            mismatches = MISMATCHES_THRESHOLD + 1;
        else
            mismatches = number_of_mismatches(test, reference.substr(pos, test.length()));       
        for (int i = bl->start + 1; i <= bl->end; ++i) {
            int new_pos = location_array[i];
            if (new_pos + test.length() > reference.length())
                continue;

            int new_mismatches = number_of_mismatches(test, reference.substr(new_pos, test.length()));       
            if (new_mismatches < mismatches) {
                mismatches = new_mismatches;
                pos = new_pos;
            }
        }

        if (mismatches > MISMATCHES_THRESHOLD) {
            cout << "No matches found" << endl;
            return;
        }


        // calculate score of match 
        // for start do not use minimum edit distance but assume only mismatches can be made and no insertions or deletions
        cout << "Number of mismatches: " << mismatches << " for read with id: " << id << " at location: " << pos << endl;
    }
}

// Function to reverse a string and take its complement
// TODO make this faster!
string reverse_complement(string const &s)
{
    char *rev = new char[s.length()];
    for (unsigned long i = 0; i < s.length(); ++i)
        rev[i] = map_complement[s[s.length() - 1 - i]];

    
	/* string rev(s1.rbegin(), s1.rend()); */
	return string(rev);
}

int main(int argc, char** argv) {

    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <reference_gene_file> <reads_gene_file_1> <reads_gene_file2>" << endl;
        exit(1);
    }

    string reference = read_reference_gene(argv[1]);
    int *location_array = new int[reference.length()];
    
    string encoded;

    encode(reference, encoded, location_array);

    cout << encoded << endl;
    string decoded = decode(encoded);

    ifstream file_1;
    getFile(argv[2], file_1);
    ifstream file_2;
    getFile(argv[3], file_2);


    // for both forward and abckward extract a seed of size 12 and match both forward and backward seed to reference
    // based on these matches -> see if the whole read matches by employing edit distance for oen read as well as both
    // forward an backward read!

    /* string seed = rb->forward_read.substr(0, SEED_LENGTH); */

    // assume forward direction first for both (ny iverting the second) and then the opposite!
    
    // TODO also think abou how to nake parallel
    // TODO also think about how to efficiently store the genome
    // TODO suffix arrays can be built in O(n) ... can we use the same to construct bwt in O(n) ? 
    // TODO founding a match with one orientation means that we do not need to check the other direction!

    string seed, reverse_complement_str;
    block * bl;
    while (true) {
        read_block* rb = get_one_read(file_1, file_2);
        if (rb == NULL)
            break;

        cout << "---------------------------------------------" << endl;
        cout << rb->id << endl;

        seed = rb->forward_read.substr(0, SEED_LENGTH);
        bl = get_matches(encoded, seed);
        print_matches(bl, rb->id, rb->forward_read, location_array, reference);

        seed = rb->backward_read.substr(0, SEED_LENGTH);
        bl = get_matches(encoded, seed);   
        print_matches(bl, rb->id, rb->backward_read, location_array, reference);

        reverse_complement_str = reverse_complement(rb->forward_read);
        seed = reverse_complement_str.substr(0, SEED_LENGTH);
        bl = get_matches(encoded, seed);
        print_matches(bl, rb->id, reverse_complement_str, location_array, reference);

        reverse_complement_str = reverse_complement(rb->backward_read);
        seed = reverse_complement_str.substr(0, SEED_LENGTH);
        bl = get_matches(encoded, seed);   
        print_matches(bl, rb->id, reverse_complement_str, location_array, reference);
    }

    return 0;
        



    
    /* clock_t begin = clock(); */
    /* clock_t end = clock(); */
    /* double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC; */
    /* cout << "first" << elapsed_secs << endl; */

    /* string compressed = compress(encoded); */
    /* cout << "Reference size: " << reference.length() << " and compressed size: " << compressed.length() << endl; */   
    /* string uncompressed = uncompress(compressed); */
    /* assert(!encoded.compare(uncompressed)); */

    /* string decoded = decode(encoded); */
    /* assert(!reference.compare(decoded)); */

    /* string test = get_one_read(argv[2]); */

    /* block* bl = get_matches(encoded, test); */   
    /* if (bl->start > bl->end) */
    /*     cout << "No matches found" << endl; */
    /* else */
    /*     cout << "String searching found for: " << bl->start << ' ' << bl->end << endl; */
}