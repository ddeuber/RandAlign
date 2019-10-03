#ifndef SAMFILE_H
#define SAMFILE_H

#include <string>
#include <fstream>

class SAMFile {
	private:
		std::ofstream file;
		std::string refSeqName;
		std::string tab = "\t";
		std::string endline = "\n";
		
		/*
		 * The Following function adds the entry for one read to the file.
		 * 
		 * readName: name of read
		 * isRead1: true if the read is read1, false for read2
		 * readSeq: nucleotide sequence of read. IMPORTANT: if the read is a reverse read, readSeq needs to be already reversed and inverted
		 * qualSeq: quality sequence of phred scores. IMPORTANT: if the read is a reverse read, qualSeq needs to be already reversed
		 * isReverse: true if the read is a reverse read
		 * pos: position of read in reference sequence. Note that the first nucleotide has position 1. If pos=0, this means that the read is unmapped.
		 * cigar: CIGAR string 
		 * posNext: pos of paired read
		 * tempLen: length of segment from the first nucleotide of the forward read to the last nucleotide of the reverse read. 
		 * 			The absolute value of tempLen is the same for both reads, but the reverse read tempLen has a negative sign. 
		 * mapQuality = −10 log10 Pr{mapping position is wrong}, mapQuality=255 means that mapping quality is not available. default is 99
		 *
		 */
		void add_single_entry(std::string const &readName, bool isRead1, std::string const &readSeq, std::string const &qualSeq, bool isReverse, int pos, std::string const &cigar, int posNext, int tempLen, int mapQuality=99); 

		/*
		 * Generates FLAG for a read.
		 * 
		 * pos: position of read in reference sequence. Note that the first nucleotide has position 1. If pos=0, this means that the read is unmapped.
		 * isReverse: true if read is reverse read
		 * isRead1: true if the read is read1, false for read2
		 * posOtherRead: pos of the other read of the pair
		*/	
		int generate_flag(int pos, bool isReverse, bool isRead1, int posOtherRead);


	public:
		/*
		 * Initializes SAMFile and prints header. 
		 * 
		 * fileName: desired name of the SAM file. Filetype (e.g. ".sam") is not added automatically but needs to be written in fileName.
		 * refSeqName: name of the reference sequence
		 * refSeqLength: legnth of the reference sequence
		 */
		SAMFile(std::string const &fileName, std::string const &refSeqName, int refSeqLength);

		/*
		 * File needs to be closed in order to flush the changes to disk.
		 */
		void close();

		/*
		 * The Following function adds two lines corresponding to a read pair to the SAM file.
		 * 
		 * readName: name of read
		 * readSeq1: nucleotide sequence of read1. IMPORTANT: if read1 is a reverse read, readSeq1 needs to be already reversed and inverted
		 * qualSeq1: quality sequence of phred scores of read1. IMPORTANT: if read1  is a reverse read, qualSeq1 needs to be already reversed
		 * pos1: position of read1 in reference sequence. Note that the first nucleotide has position 1. If pos=0, this means that the read is unmapped.
		 * cigar1: CIGAR string of read1
		 * readSeq2: nucleotide sequence of read2. IMPORTANT: if read2 is a reverse read, readSeq2 needs to be already reversed and inverted
		 * qualSeq2: quality sequence of phred scores of read2. IMPORTANT: if read2  is a reverse read, qualSeq2 needs to be already reversed
		 * pos2: position of read2 in reference sequence. Note that the first nucleotide has position 1. If pos=0, this means that the read is unmapped.
		 * cigar2: CIGAR string of read2	
		 * read1Reversed: true if read1 is the reverse read, false if read2 is the reverse read 
		 * mapQuality1/2 = −10 log10 Pr{mapping position is wrong} of read1/2, mapQuality=255 means that mapping quality is not available. default is 99
		 *
		 */
		void add_paired_read_entry(std::string const &readName, std::string const &readSeq1, std::string const &quallSeq1, int pos1, std::string const &cigar1,
													std::string const &readSeq2, std::string const &qualSeq2, int pos2, std::string const &cigar2,
													bool read1Reversed, int mapQuality1=99, int mapQuality2=99);

		/*
		 * Converts a global alignment between a part of the reference sequence (dnaAligned) and the read (readAligned) into a CIGAR string. 
		 */	
		static std::string alignment_to_CIGAR(std::string const &dnaAligned, std::string const &readAligned);
};

#endif
