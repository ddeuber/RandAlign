#ifndef RAND_ALIGN_H 
#define RAND_ALIGN_H 

#include "bwt.h"
#include <string>
#include <vector>
#include "samfile.h"
#include "read_genes.h"

int cost(char a, char b);
int score(char a, char b);

// returns edit distance and aligned strings are saved in s1aligned and a2aligned
int global_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string &s2Aligned); 

int quasi_local_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string &s2Aligned); 
			 
class RandomizedAligner {
	private:
		BWT* bwt;
		SAMFile* samFile;
	
	public:
		RandomizedAligner(BWT* bwt, SAMFile* samfile);

		// returns position of possible match, -1 if nothing is found. The according cigar string is stored in cigarOutput, and the edit distance in editDistance
		// Note that the actual mean of the seeds will be meanSeedLength + 5
		int get_alignment_candidate(std::string const& read, int meanSeedLength, std::string& cigarOutput, int &editDistance);

		// align and print into SAMFile
		void align_and_print(read_block* rb, int maxIter=100);
};

#endif
