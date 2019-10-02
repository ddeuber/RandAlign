#ifndef BWT_H
#define BWT_H 

#include <map>

#define SEED_LENGTH 12
#define NUM_BASES 5
#define MISMATCHES_THRESHOLD 5

// used to return the indices for the matching lines found
typedef struct { 
    int start;
    int end;
} block;


// map from symbols to ints
extern std::map<char, int> mapOfGenes;
// map of complementary genes
extern std::map<char, char> map_complement;

class BWT {
	private:
        int* location_array;
        std::string encoded;

        void encode();
        // TODO decode should also get the location array matrix (or recover is threw a lot of computation..)
        void decode();
	public:
        // create BWT by providing either the reference (is_encoded = false) or the encoded string (is_encoded = true)
        BWT(std::string rawString, bool is_encoded);

        std::string reference;
        
        block* get_matches(std::string test);

        // returns location in original string based on the location array 
        int get_location(int pos);
};

#endif
