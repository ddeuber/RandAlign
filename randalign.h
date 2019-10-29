#ifndef RAND_ALIGN_H 
#define RAND_ALIGN_H 

#include "bwt.h"
#include <string>
#include <vector>
#include "samfile.h"
#include "read_genes.h"

typedef struct{
	bool failed; 
	// Read 1
	std::string readSeq1; 
	std::string qualSeq1; 
	int pos1; 
	std::string cigar1; 
	// Read 2
	std::string readSeq2; 
	std::string qualSeq2; 
	int pos2; 
	std::string cigar2; 
	bool read1Reversed; 
} results_block; 

/*
 * Score for local alignment.
 */
int score(char a, char b);

/*
 * The string s1 should be longer than s2.
 * The following algorithm is a slight adaptation of the Smith–Waterman algorithm.
 * The difference is that the whole sequence s2 is aligned to s1, not only a substring.
 * Also there will not be any insertions at the beginning of s2.
 * The function returns (negative score + readlength), so this can be used as an edit distance.
 */
int quasi_local_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string &s2Aligned); 


class RandomizedAligner {
	private:
		BWT* bwt;
		SAMFile* samFile;
	
	public:
		RandomizedAligner(BWT* bwt, SAMFile* samfile);

		// returns position of possible match, -1 if nothing is found. The according cigar string is stored in cigarOutput, and the edit distance in editDistance
		// Note that the actual mean of the seeds will be meanSeedLength + 10 
		int get_alignment_candidate(std::string const& read, std::string const& qualString, int meanSeedLength, int maxDist, std::string& cigarOutput, int &editDistance, bool mismatchOnly);
		void get_alignment_candidates(std::string const& read, std::string const& qualString, int meanSeedLength, int maxShift, int maxDistance, std::vector<int> &refPositions, std::vector<std::string> &cigarStrings, std::vector<int> &editDistances, bool mismatchOnly);

		// align and print into SAMFile
		results_block* align_and_print(read_block* rb, int maxIter=100);
};

#endif
