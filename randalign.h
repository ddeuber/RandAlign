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
    /*
     * RandomizedAligner finds matching positions for reads by taking two random seeds and matching them to the reference genome. 
     * In particular, get_alignment_candidate is a function which applies this randomized match-finding for one read and returns one of the matches,
     * where a match with a lower edit distance has a higher probability of being returned.
     * This is then used in align_and_print for both read1 and read2 to get possible matching positions, and the position is accepted
     * if the forward and the reverse read have a distance of no more than 450+readlength. 
     */
	private:
		BWT* bwt;
		SAMFile* samFile;

		/*  
         *  Returns position of possible match, -1 if nothing is found. The according cigar string is stored in cigarOutput, and the edit distance in editDistance.
         *  If mismatchOnly is true, insertions and deletions will not be considered for the alignment (this makes it much faster since we do not need to call quasi_local_alignment).
		 *  Note that the actual mean of the seeds will be meanSeedLength + 10.
         */
		int get_alignment_candidate(std::string const& read, std::string const& qualString, int meanSeedLength, int maxDist, std::string& cigarOutput, int &editDistance, bool mismatchOnly);
	
	public:
		RandomizedAligner(BWT* bwt, SAMFile* samfile);

		/* 
         * Aligns and prints into SAMFile one read pair.
         */
		results_block* align_and_print(read_block* rb, int maxIter=100);
};

#endif
