#include "randalign.h"
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "samfile.h"
#include "read_genes.h"


int score(char a, char b){
	if (a==b)
		return 1;
	else if (a=='-' || b=='-')
		return -1;
	else
		return 0;
}	

int quasi_local_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string& s2Aligned) {
	int table[s1.length()+1][s2.length()+1];

	int IMIN = std::numeric_limits<int>::min();

	// initialize border
	table[0][0] = 0;
	for (unsigned int i=1; i<=s1.length(); ++i)
		table[i][0] = 0;
	for (unsigned int i=1; i<=s2.length(); ++i)
		table[0][i] = IMIN/2; //initialize to a very small number, so that there is no insertion at the begging of s2

	int v1, v2, v3;

	// fill table as in Smith-Waterman
	for (unsigned int i=1; i<=s1.length(); ++i){
		for (unsigned int j=1; j<=s2.length(); ++j){
			v1 = table[i-1][j-1] + score(s1[i-1], s2[j-1]);
			v2 = table[i-1][j] + score(s1[i-1], '-');
			v3 = table[i][j-1] + score('-', s2[j-1]);
			table[i][j] = std::max(std::max(v1,v2), v3);
		}
	}

	int maxscore = IMIN;
	int value; 
	int maxindexi = -1;
	int maxindexj = -1;

    // find maximal score
	for (int i=0; i<=s1.length(); ++i){
		for (int j=0; j<=s2.length(); ++j) {
			value = table[i][j];
			if (value > maxscore) {
				maxscore = value;
				maxindexi = i;
				maxindexj = j;
			}
		}
	}
	
	// find matchoptimal alignment by backtracking
	int v;
	int i = maxindexi;
	int j = maxindexj;
	std::string a1 = "";
	std::string a2 = "";

	while(j > 0 && i>=0) {
		v = table[i][j];
		
		v1 = (i>0 && j>0) ? table[i-1][j-1] + score(s1[i-1], s2[j-1]) : IMIN;
		v2 = (i>0) ? table[i-1][j] + score(s1[i-1], '-') : IMIN;
		v3 = (j>0) ? table[i][j-1] + score('-', s2[j-1]) : IMIN;

		if (v == v1) {
			a1 += s1[i-1];
			a2 += s2[j-1];
			--i;
			--j;
		} else if (v == v2) {
			a1 += s1[i-1];
			a2 += '-';
			--i;
		} else {
			a1 += '-';
			a2 += s2[j-1];
			--j;
		}
	}
	
    // construct alignment strings
	std::reverse(a1.begin(), a1.end());
	std::reverse(a2.begin(), a2.end());

	s1Aligned = s1.substr(0, i) + a1 + s1.substr(maxindexi);
	s2Aligned = std::string(i, '-') + a2 + s2.substr(maxindexj);
	s2Aligned.resize(s1Aligned.length(), '-');

	return -maxscore + (int)s2.length(); 
}

RandomizedAligner::RandomizedAligner(BWT* bwt, SAMFile* samFile){
	this->bwt = bwt;
	this->samFile = samFile;
}

int RandomizedAligner::get_alignment_candidate(std::string const& read, std::string const &qualString, int meanSeedLength, int maxDist, std::string &cigarOutput, int &editDistance, bool mismatchOnly){
	int n = read.length();
	double p = 1.0 * meanSeedLength / n;
	int N_SEEDS = 2;
	
	// get binomially distributed seed length
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::binomial_distribution<int> binom(n, p);
	int seedLengths[N_SEEDS];
	for (unsigned int i=0; i<N_SEEDS; ++i)
		seedLengths[i] = binom(generator) + 10; // seeds must not get too small

	// get random seeds, where the probability for a seed is higher if it has higher minimal phred score on its positions
	std::discrete_distribution<int> distribution;
	std::vector<int> weights;
	int seedStarts[N_SEEDS];
	
	int minval = 100;

	do {
		for (unsigned int i=0; i<N_SEEDS; ++i) {
			weights = std::vector<int>(n-seedLengths[i]-1);
			for (int j=0; j<weights.size(); ++j){
				for (int k=0; k<seedLengths[i]; ++k){
					if (qualString[j+k] < minval)
						minval = qualString[j+k];
				weights[j] = minval - 33; //minimal phred score
				weights[j] = weights[j] * weights[j];
				}
			}
			distribution = std::discrete_distribution<int>(weights.begin(), weights.end());
			seedStarts[i] = distribution(generator);
		}
	} while( std::abs(seedStarts[0] - seedStarts[1]) < 20);


	// exchange seeds such that first one is the first in the read
	if (seedStarts[0] > seedStarts[1]) {
		std::swap(seedStarts[0], seedStarts[1]);
		std::swap(seedLengths[0], seedLengths[1]);
	}

	// get matches from BWT index
	block* bl; 
	std::vector<int> positions[N_SEEDS];
	for (unsigned int i=0; i<N_SEEDS; ++i){
		bl = bwt->get_matches(read.substr(seedStarts[i], seedLengths[i]));
        for (int j = bl->start; j <= bl->end; ++j){
			positions[i].push_back(bwt->get_location(j));
		}
		std::sort(positions[i].begin(), positions[i].end());
	}

	// take matches where the distance betweend seed1 and seed2 is more or less preserved (up to MAX_SHIFT)
	// Note that this works because the positions are sorted.
	int MAX_SHIFT = mismatchOnly ? 0 : 3;
	std::vector<int> pairPositions;
	int i=0;
	int j=0;
	int p0, p1;
	
	while(i<positions[0].size() && j<positions[1].size()){
		p0 = positions[0][i];
		p1 = positions[1][j];
		
		if (p1-p0 < seedStarts[1] - seedStarts[0] - MAX_SHIFT) // p1 before p0 or distance is to small
			++j;
		else if (p1-p0 <= seedStarts[1] - seedStarts[0] + MAX_SHIFT) { //distance is right
			pairPositions.push_back(p0);
			++i;
		} else // distance is to big
			++i;
	}

	if (pairPositions.empty())
		return -1;
	
	// compare edit distance of found matches
	int npos = pairPositions.size();
	int editDistances[npos];
	weights = std::vector<int>(npos); // weights which are then used to sample one of the matches
	int refPos[npos];
	std::string dnaAligned[npos];
	std::string readAligned[npos];
	std::string refSeq;

	for (unsigned int i=0; i<npos; ++i){
		if (mismatchOnly){
			refPos[i] = std::max(pairPositions[i]-seedStarts[0], 0); 
			refSeq = bwt->reference.substr(refPos[i], read.length());
			editDistances[i] = number_of_mismatches(refSeq, read); 
			weights[i] = std::max((int)read.length()-(int)read.length()/maxDist*(editDistances[i]), 0);
			dnaAligned[i] = refSeq;
			readAligned[i] = read;
		} else {
			refPos[i] = std::max(pairPositions[i]-seedStarts[0]-5, 0); 
			refSeq = bwt->reference.substr(refPos[i], read.length()+10);
			editDistances[i] = quasi_local_alignment(refSeq, read, dnaAligned[i], readAligned[i]);
			weights[i] = std::max((int)read.length()-(int)read.length()/maxDist*(editDistances[i]), 0);
		}
	}

	// check that there is a strictly positive weight, otherwise return
	bool allzero = true;
	for (auto i : weights){
		if (i!=0) {
			allzero = false;
			break;
		}
	}
	
	if (allzero)
		return -1;

    // One of the possible read matches is sampled, where matches with lower edit distance are more likely,
    // and if editDistance >= maxDistance, the probability for the match is zero
	std::discrete_distribution<int> categorical(weights.begin(), weights.end());
	int match = categorical(generator);
	editDistance = editDistances[match];	
	
	// find position where read starts in alignment
	int readStart = readAligned[match].find_first_not_of('-');
	int readEnd = readAligned[match].find_last_not_of('-') + 1;
	int readSize = readEnd - readStart;
	cigarOutput = SAMFile::alignment_to_CIGAR(dnaAligned[match].substr(readStart, readSize), readAligned[match].substr(readStart, readSize));
	return bwt->convert_index_to_original_index(refPos[match] + readStart);
}

results_block* RandomizedAligner::align_and_print(read_block* rb, int maxIter){
	// Return variable 
	results_block* res_block = new results_block; 

	int pos1, pos1rev;
	int pos2, pos2rev;
	int dist1, dist1rev;
	int dist2, dist2rev;
	int distsum1, distsum2;
	std::string cigar1, cigar1rev;
	std::string cigar2, cigar2rev;
	
	std::string &read1 = rb->forward_read;
	std::string &read2 = rb->backward_read;
	std::string &qualSeq1 = rb->forward_read_quality;
	std::string &qualSeq2 = rb->backward_read_quality;

	std::string rev1 = reverse_complement(read1);
	std::string rev2 = reverse_complement(read2);
	std::string qualSeq1rev = qualSeq1; 
	std::reverse(qualSeq1rev.begin(), qualSeq1rev.end());
	std::string qualSeq2rev = qualSeq2;
	std::reverse(qualSeq2rev.begin(), qualSeq2rev.end());

	bool mismatchOnly;
	int meanSeedLength = 15;
	int maxDist;

	for(int i=0; i<maxIter; ++i){
		// try out both reads as reverse reads
		mismatchOnly = i < maxIter/2;

		if (mismatchOnly)
			maxDist = std::min(i/5 + 1, 10);
		else
			maxDist = std::min((i-maxIter)/5 + 1, 10);

		pos1 = get_alignment_candidate(read1, qualSeq1, meanSeedLength, maxDist, cigar1, dist1, mismatchOnly);
		pos2 = get_alignment_candidate(read2, qualSeq2, meanSeedLength, maxDist, cigar2, dist2, mismatchOnly);
		
		pos1rev = get_alignment_candidate(rev1, qualSeq1rev, meanSeedLength, maxDist, cigar1rev, dist1rev, mismatchOnly);
		pos2rev = get_alignment_candidate(rev2, qualSeq2rev, meanSeedLength, maxDist, cigar2rev, dist2rev, mismatchOnly);

        // In case a match was found for one read but not its reverse read, try again
		for (int j=0; j<5 && pos1 >= 0 && pos2rev == -1; ++j)
			pos2rev = get_alignment_candidate(rev2, qualSeq2rev, meanSeedLength, maxDist, cigar2rev, dist2rev, mismatchOnly);
		for (int j=0; j<5 && pos2rev >= 0 && pos1 == -1; ++j)
			pos1 = get_alignment_candidate(read1, qualSeq1, meanSeedLength, maxDist, cigar1, dist1, mismatchOnly);
		for (int j=0; j<5 && pos2 >= 0 && pos1rev == -1; ++j)
			pos1rev = get_alignment_candidate(rev1, qualSeq1rev, meanSeedLength, maxDist, cigar1rev, dist1rev, mismatchOnly);
		for (int j=0; j<5 && pos1rev >= 0 && pos2 == -1; ++j)
			pos2 = get_alignment_candidate(read2, qualSeq2, meanSeedLength, maxDist, cigar2, dist2, mismatchOnly);

		bool read1rev2 = (pos1 > 0) && (pos2rev>0) && std::abs(pos1 - pos2rev) < 450 + read1.length(); // if combination read1 rev2 is possible
		bool read2rev1 = (pos2 > 0) && (pos1rev>0) && std::abs(pos2 - pos1rev) < 450 + read1.length(); // if combination read2 rev1 is possible
	
		if (!read1rev2 && !read2rev1)
			continue;
		else {
			if (read1rev2 && read2rev1) { //if both combinations are possible take one with lower sum of edit distance
				distsum1 = dist1 + dist2rev;
				distsum2 = dist2 + dist1rev;
			} else if (read1rev2) {
				distsum1 = 0;
				distsum2 = 1;
			} else {
				distsum1 = 1;
				distsum2 = 0;
			}
			
			if(distsum1 <= distsum2){
				std::reverse(qualSeq2.begin(), qualSeq2.end());
				// Write to output 
				res_block->failed = false; 
				// Read 1
				res_block->readSeq1 = read1; 
				res_block->qualSeq1 = qualSeq1; 
				res_block->pos1 = pos1+1;
				res_block->cigar1 = cigar1; 

				//Read 2
				res_block->readSeq2 = rev2; 
				res_block->qualSeq2 = qualSeq2; 
				res_block->pos2 = pos2rev+1;
				res_block->cigar2 = cigar2rev; 

				res_block->read1Reversed = false; 

				return res_block;
			} else {
				std::reverse(qualSeq1.begin(), qualSeq1.end());
				// Write to output 
				res_block->failed = false; 
				// Read 1
				res_block->readSeq1 = rev1; 
				res_block->qualSeq1 = qualSeq1; 
				res_block->pos1 = pos1rev+1;
				res_block->cigar1 = cigar1rev; 

				//Read 2
				res_block->readSeq2 = read2; 
				res_block->qualSeq2 = qualSeq2; 
				res_block->pos2 = pos2+1;
				res_block->cigar2 = cigar2; 

				res_block->read1Reversed = true; 

				return res_block; 
			}
		}
	}

	// If failed

	// Write to output 
	res_block->failed = true; 
	// Read 1
	res_block->readSeq1 = read1; 
	res_block->qualSeq1 = qualSeq1; 
	res_block->pos1 = 0;
	res_block->cigar1 = "*"; 

	//Read 2
	res_block->readSeq2 = rev2; 
	res_block->qualSeq2 = qualSeq2; 
	res_block->pos2 = 0;
	res_block->cigar2 = "*"; 

	res_block->read1Reversed = false; 
	
	// std::cout << "Failed for read " << rb->id << std::endl;
	// std::cout << pos1 << " " << pos2rev << "\t" << cigar1 << " " << cigar2rev << "\t" << dist1 << " " << dist2rev << std::endl;
	// std::cout << pos2 << " " << pos1rev << "\t" << cigar2 << " " << cigar1rev << "\t" << dist2 << " " << dist1rev << std::endl;
	// std::cout << std::endl;
	// samFile->add_paired_read_entry(rb->id, read1, qualSeq1, 0, "*", rev2, qualSeq2, 0, "*", false);

	return res_block;
}

