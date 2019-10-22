#include "randalign.h"
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include "samfile.h"
#include "read_genes.h"

int cost(char a, char b) {
	if (a==b)
		return 0;
	else if (a=='-' || b=='-')
		return 2;
	else 
		return 1;
}

int score(char a, char b){
	if (a==b)
		return 1;
	else if (a=='-' || b=='-')
		return -1;
	else
		return 0;
}	

// returns negative score to be compatible with editDistance
// s1 should be longer than s2
int quasi_local_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string& s2Aligned) {
	int table[s1.length()+1][s2.length()+1];

	int IMIN = std::numeric_limits<int>::min();

	// initialize border
	table[0][0] = 0;
	for (unsigned int i=1; i<=s1.length(); ++i)
		table[i][0] = 0;
	for (unsigned int i=1; i<=s2.length(); ++i)
		table[0][i] = IMIN/2;

	int v1, v2, v3;

	// fill table
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
	
	// find optimal alignment
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
	
	std::reverse(a1.begin(), a1.end());
	std::reverse(a2.begin(), a2.end());

	s1Aligned = s1.substr(0, i) + a1 + s1.substr(maxindexi);
	s2Aligned = std::string(i, '-') + a2 + s2.substr(maxindexj);
	s2Aligned.resize(s1Aligned.length(), '-');

	return maxscore; 
}

int global_alignment(std::string const &s1, std::string const &s2, std::string& s1Aligned, std::string& s2Aligned) {
	// For simplicity, the Needleman-Wunsch algorithm is used
	int table[s1.length()+1][s2.length()+1];

	int IMAX = std::numeric_limits<int>::max();

	// initialize border
	table[0][0] = 0;
	for (unsigned int i=1; i<=s1.length(); ++i)
		table[i][0] = table[i-1][0] + cost('A', '-');
	for (unsigned int i=1; i<=s2.length(); ++i)
		table[0][i] = table[0][i-1] + cost('-', 'A');

	int v1, v2, v3;

	// fill table
	for (unsigned int i=1; i<=s1.length(); ++i){
		for (unsigned int j=1; j<=s2.length(); ++j){
			v1 = table[i-1][j-1] + cost(s1[i-1], s2[j-1]);
			v2 = table[i-1][j] + cost(s1[i-1], '-');
			v3 = table[i][j-1] + cost('-', s2[j-1]);
			table[i][j] = std::min(std::min(v1,v2), v3);
		}
	}
	
	// find optimal alignment
	int v;
	int i = s1.length();
	int j = s2.length();
	std::string a1 = "";
	std::string a2 = "";

	while(i > 0  || j > 0) {
		v = table[i][j];
		
		v1 = (i>0 && j>0) ? table[i-1][j-1] + cost(s1[i-1], s2[j-1]) : IMAX;
		v2 = (i>0) ? table[i-1][j] + cost(s1[i-1], '-') : IMAX;
		v3 = (j>0) ? table[i][j-1] + cost('-', s2[j-1]) : IMAX;

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
	
	std::reverse(a1.begin(), a1.end());
	std::reverse(a2.begin(), a2.end());

	s1Aligned = a1;
	s2Aligned = a2;
	return table[s1.length()][s2.length()];
}


RandomizedAligner::RandomizedAligner(BWT* bwt, SAMFile* samFile){
	this->bwt = bwt;
	this->samFile = samFile;
}

int RandomizedAligner::get_alignment_candidate(std::string const& read, int meanSeedLength, std::string &cigarOutput, int &editDistance, bool mismatchOnly){
	int n = read.length();
	double p = 1.0 * meanSeedLength / n;
	int N_SEEDS = 2;
	int MAX_OVERLAP = 10;
	
	// get binomially distributed seed length
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
	std::default_random_engine generator(seed);
	std::binomial_distribution<int> binom(n, p);
	int seedLengths[N_SEEDS];
	for (unsigned int i=0; i<N_SEEDS; ++i)
		seedLengths[i] = binom(generator) + 5; //such that seeds to not get too small

	// get random seeds
	std::uniform_int_distribution<int> unif;
	int seedStarts[N_SEEDS];
	do {
		for (unsigned int i=0; i<N_SEEDS; ++i) {
			unif = std::uniform_int_distribution<int>(0, n-seedLengths[i]-1);
			seedStarts[i] = unif(generator);
		}
	} while( std::abs(seedStarts[0] - seedStarts[1]) < 20);

	// exchange seeds such that first one is the first in the read
	if (seedStarts[0] > seedStarts[1]) {
		std::swap(seedStarts[0], seedStarts[1]);
		std::swap(seedLengths[0], seedLengths[1]);
	}

	// get matches
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
	int MAX_SHIFT = mismatchOnly ? 0 : 15;
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
	
	// compare edit distance of found matches, find optimal one
	int npos = pairPositions.size();
	int editDistances[npos];
	int refPos[npos];
	std::string dnaAligned[npos];
	std::string readAligned[npos];
	std::string refSeq;

	for (unsigned int i=0; i<npos; ++i){
		if (mismatchOnly){
			refPos[i] = std::max(pairPositions[i]-seedStarts[0], 0); 
			refSeq = bwt->reference.substr(refPos[i], read.length());
			editDistances[i] = number_of_mismatches(refSeq, read); 
			dnaAligned[i] = refSeq;
			readAligned[i] = read;
		} else {
			refPos[i] = std::max(pairPositions[i]-seedStarts[0]-15, 0); 
			refSeq = bwt->reference.substr(refPos[i], read.length()+30);
			// editDistances[i] = global_alignment(refSeq, read, dnaAligned[i], readAligned[i]);
			editDistances[i] = quasi_local_alignment(refSeq, read, dnaAligned[i], readAligned[i]);
		}
	}

	int* minDist = std::min_element(editDistances, editDistances + npos);
	editDistance = *minDist;
	int optimalMatch = minDist - editDistances; 
	
	// if mismatchesOnly, return nothing found if there are more than 5 mismatches
	if (mismatchOnly) {
		if (editDistance > 5)
			return -1;
	}

	// find position where read starts in alignment
	int readStart = readAligned[optimalMatch].find_first_not_of('-');
	int readEnd = readAligned[optimalMatch].find_last_not_of('-') + 1;
	int readSize = readEnd - readStart;
	cigarOutput = SAMFile::alignment_to_CIGAR(dnaAligned[optimalMatch].substr(readStart, readSize), readAligned[optimalMatch].substr(readStart, readSize));
	return refPos[optimalMatch] + readStart;
}

void RandomizedAligner::align_and_print(read_block* rb, int maxIter){
	
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

	bool mismatchOnly;
	int meanSeedLength = 15;

	for(int i=0; i<maxIter; ++i){
		// try out both reads as reverse reads
		mismatchOnly = i < maxIter/2;

		pos1 = get_alignment_candidate(read1, meanSeedLength, cigar1, dist1, mismatchOnly);
		pos2 = get_alignment_candidate(read2, meanSeedLength, cigar2, dist2, mismatchOnly);
		
		pos1rev = get_alignment_candidate(rev1, meanSeedLength, cigar1rev, dist1rev, mismatchOnly);
		pos2rev = get_alignment_candidate(rev2, meanSeedLength, cigar2rev, dist2rev, mismatchOnly);

		for (int j=0; j<5 && pos1 >= 0 && pos2rev == -1; ++j)
			pos2rev = get_alignment_candidate(rev2, meanSeedLength, cigar2rev, dist2rev, mismatchOnly);
		for (int j=0; j<5 && pos2rev >= 0 && pos1 == -1; ++j)
			pos1 = get_alignment_candidate(read1, meanSeedLength, cigar1, dist1, mismatchOnly);
		for (int j=0; j<5 && pos2 >= 0 && pos1rev == -1; ++j)
			pos1rev = get_alignment_candidate(rev1, meanSeedLength, cigar1rev, dist1rev, mismatchOnly);
		for (int j=0; j<5 && pos1rev >= 0 && pos2 == -1; ++j)
			pos2 = get_alignment_candidate(read2, meanSeedLength, cigar2, dist2, mismatchOnly);

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
				samFile->add_paired_read_entry(rb->id, read1, qualSeq1, pos1+1, cigar1, rev2, qualSeq2, pos2rev+1, cigar2rev, false);
				return;
			} else {
				std::reverse(qualSeq1.begin(), qualSeq1.end());
				samFile->add_paired_read_entry(rb->id, rev1, qualSeq1, pos1rev+1, cigar1rev, read2, qualSeq2, pos2+1, cigar2, true);
				return;
			}
		}
	}
	
	// if failed.
	std::cout << "Failed for read " << rb->id << std::endl;
	std::cout << pos1 << " " << pos2rev << "\t" << cigar1 << " " << cigar2rev << std::endl;
	std::cout << pos2 << " " << pos1rev << "\t" << cigar2 << " " << cigar1rev << std::endl;
	samFile->add_paired_read_entry(rb->id, read1, qualSeq1, 0, "*", rev2, qualSeq2, 0, "*", false);
}
 
