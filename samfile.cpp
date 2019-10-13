#include "samfile.h"
#include <iostream>

SAMFile::SAMFile(std::string const &fileName, std::string const &refSeqName, int refSeqLength){
	file.open(fileName);
	file << "@HD\tN:1.4\tSO:unsorted";
    for (int i=0; i<8; ++i)
        file << tab;
    file << endline;

	file << "@SQ\tSN:" << refSeqName << "\tLN:" << refSeqLength;
    for (int i=0; i<8; ++i)
        file << tab;
    file << endline;

	this->refSeqName = refSeqName;
}

void SAMFile::close(){
	file.close();
}

void SAMFile::add_paired_read_entry(std::string const &readName, std::string const &readSeq1, std::string const &qualSeq1, int pos1, std::string const &cigar1,
										std::string const &readSeq2, std::string const &qualSeq2, int pos2, std::string const &cigar2,
										bool read1Reversed, int mapQuality1, int mapQuality2){
	int tempLen;
	if (read1Reversed) {
		tempLen = pos1 - pos2 + readSeq1.length();
		add_single_entry(readName, true, readSeq1, qualSeq1, true, pos1, cigar1, pos2, -tempLen, mapQuality1);
		add_single_entry(readName, false, readSeq2, qualSeq2, false, pos2, cigar2, pos1, tempLen, mapQuality2);
	} else {
		tempLen = pos2 - pos1 + readSeq2.length();
		add_single_entry(readName, true, readSeq1, qualSeq1, false, pos1, cigar1, pos2, tempLen, mapQuality1);
		add_single_entry(readName, false, readSeq2, qualSeq2, true, pos2, cigar2, pos1, -tempLen, mapQuality2);	
	}
}

void SAMFile::add_single_entry(std::string const &readName, bool isRead1, std::string const &readSeq, std::string const &qualSeq, bool isReverse, int pos, std::string const &cigar, int posNext, int tempLen, int mapQuality){
	file << readName << tab;
	file << generate_flag(pos, isReverse, isRead1, posNext) << tab; 
	file << refSeqName << tab;
	file << pos << tab;
	file << mapQuality << tab;
	file << cigar << tab;
	file << "=" << tab;
	file << posNext << tab;
	file << tempLen << tab;
	file << readSeq << tab;
	file << qualSeq << endline;
}

int  SAMFile::generate_flag(int pos, bool isReverse, bool isRead1, int posOtherRead){
	int  flag = 1;
	
	if (pos > 0 && posOtherRead > 0)
		flag |= 2;
	if (pos == 0)
		flag |= 4;
	if (posOtherRead == 0)
		flag |= 8;

	if (isReverse)
		flag |= 16;
	else
		flag |= 32;

	if (isRead1)
		flag |= 64;
	else
		flag |= 128;
		
	return flag;
}

std::string SAMFile::alignment_to_CIGAR(std::string const &dnaAligned, std::string const &readAligned){
	std::string cifar = "";
	char lastmode = '!'; //one of 'M' (match), 'I' (insert in reference), 'D' (deletion in reference), 'X' mismatch
	char mode; 
	int count; 
	char c1, c2;
		
	for (unsigned int i=0; i < dnaAligned.length(); ++i){
		c1 = dnaAligned[i];
		c2 = readAligned[i];
		
		if (c1 == c2)
			mode = 'M';
		else if (c1 == '-')
			mode = 'I';
		else if (c2 == '-')
			mode = 'D';
		else
			mode = 'X';
		
		if (mode == lastmode) {
			++count;
		} else if (lastmode == '!'){
			lastmode = mode;
			count = 1;
		} else {
			cifar += std::to_string(count);
			cifar += lastmode;
			lastmode = mode;
			count = 1;
		}
	}

	cifar += std::to_string(count);
	cifar += lastmode;
	return cifar;
}


