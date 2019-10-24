#ifndef BWT_H
#define BWT_H 

#include <map>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>

// files where the index will be stored
#define LOCATION_ARRAY_FILENAME "location_array"
#define COUNTS_ARRAY_FILENAME "counts_array"
#define OCCURRENCES_MATRIX_FILENAME "occurrences_matrix"
#define HOLES_FILENAME "holes"
#define FILENAME_TYPE ".save"

#define SEED_LENGTH 12
#define NUM_BASES 5
#define MISMATCHES_THRESHOLD 5

// used to return the indices for the matching lines found
typedef struct { 
    int start;
    int end;
} block;


extern std::map<char, int> mapOfGenes;

extern std::map<char, char> map_complement;

class BWT {
	private:
        int* C; 
        int** Occ; 
        
        int* location_array;
        std::string encoded;

        void encode();
        void decode();

        void recover_encoded_string(int);
        void recover_index(const std::string&);

	public:
        std::vector< std::pair<int, int> > holes;

        // create BWT and index
        BWT(std::string&, std::vector< std::pair<int, int> >);

        // create BWT by reading index from filesystem
        BWT(const std::string&);

        void store_index(const std::string&);

        std::string reference;
        
        block* get_matches(const std::string&);

        // returns location in original string based on the location array 
        int get_location(int);
};

#endif
