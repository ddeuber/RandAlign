#ifndef BWT_H
#define BWT_H 

#include <map>
#include <string>
#include <sstream>
#include <iostream>

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
	public:
        // create BWT by providing either the reference (is_encoded = false) or the encoded string (is_encoded = true)
        BWT(std::string& rawString, bool is_encoded);

        std::string reference;
        
        block* get_matches(const std::string& test);

        // returns location in original string based on the location array 
        int get_location(int pos);
};

#endif
