#include <iostream>
#include <fstream>
#include <string>
#include <string.h>
#include <assert.h>
#include <list>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <algorithm> 

#include "read_genes.h"
#include "bwt.h"

using namespace std;

void getFile(const string& filename, /*out*/ ifstream& file);

ReadGenes::ReadGenes(const string& filename_1, const string& filename_2) {
    getFile(filename_1, this->file_1);
    getFile(filename_2, this->file_2);
}


string join(const vector<string>& vec, const char* delim)
{
    stringstream res;
    copy(vec.begin(), vec.end(), ostream_iterator<string>(res, delim));
    return res.str();
}

// returns the newly read genen and
// a vector of pairs containing the starting position of each sequence of N's and the number of N's in this position
string read_reference_gene(const string& filename, string& refName, int &refLength, vector< pair<int, int> >& holes){
    string gene;
    vector<string> vector_genes;

    ifstream read_file(filename);

    if (read_file){
        getline (read_file, gene);
		refName = gene.substr(1);

        while(getline(read_file, gene)) {
            vector_genes.push_back(gene);
        }
        vector_genes.push_back("$");

        string str = join(vector_genes, "");
		refLength = str.length() - 1;

        // find potential N's in the string
        // all position should be relative to the new index
        int i = 0;
        int new_index = 0;

        while (i < str.length()){
            // check consistency of reference gene
            if (!strchr("NACGT$", str[i])){
                cout << "Invalid chacater " << str[i] << " encountered at input reference string. Exiting ..." << endl;
                exit(1);
            }

            if (str[i] == 'N'){
                int j = i + 1;
                while(j < str.length() && str[j] == 'N'){
                    j ++;
                }

                holes.push_back(make_pair(new_index, j - i));
                i = j;
            }
            else {
                ++i;
                ++new_index;
            }
        }

        str.erase(remove(str.begin(), str.end(), 'N'), str.end());
        read_file.close();
        return str;
    } else {
        cout << "File: " << filename << "cannot be opened" << endl; 
        exit(ERROR_EXIT_CODE);
    }
}

string extract_id(const string& s) {
    return s.substr(1, s.find("/") - 1);
}

// check consistency for read
void check_consistent_read(string str){
    for (int i = 0; i < str.length(); ++i){
        if (!strchr("ACGT", str[i])){
            cout << "Invalid chacater " << str[i] << " encountered at read string. Exiting ..." << endl;
            exit(1);
        }
    }
}

read_block* ReadGenes::get_one_read(){
    read_block* rb = new read_block;

    try {
        getline(file_1, rb->forward_read);
        rb->id = extract_id(rb->forward_read);

        getline(file_1, rb->forward_read);
        check_consistent_read(rb->forward_read);
        
        getline(file_1, rb->forward_read_quality);
        getline(file_1, rb->forward_read_quality);

        getline(file_2, rb->backward_read);
        assert(!extract_id(rb->backward_read).compare(rb->id));
        getline(file_2, rb->backward_read);
        check_consistent_read(rb->backward_read);

        getline(file_2, rb->backward_read_quality);
        getline(file_2, rb->backward_read_quality);
    } catch(std::out_of_range e){
        return NULL;
    }

    return rb;
}

bool checkFileExistence(const string& filename)
{
    ifstream f(filename.c_str());
    return f.is_open();
}
    
void getFile(const string& filename, /*out*/ ifstream& file)
{
    const bool file_exists = checkFileExistence(filename);
    if (!file_exists) {
        cout << "File " << filename << " not found." << endl;
        exit(ERROR_EXIT_CODE);
    }
    file.open(filename.c_str());
}

int number_of_mismatches(const string &s1, const string& s2) {

    // assume len(s1) == len(s2) 
    int mismatches = 0;
    for (unsigned long i = 0; i < s1.length(); ++i)
        if (s1[i] != s2[i])
            mismatches ++;   

    return mismatches;
}

// TODO modify this method to add to file using SAM format!
void print_matches(block* bl, const string& id, const string& test, BWT* bwt) {	
    if (bl->start > bl->end)
        cout << "No matches found" << endl;
    else {
        // assume that from all matches the one corresponding to the least amount of mismatches is the allignment 
        int pos = bwt->get_location(bl->start);
        int mismatches;
        
        // check that the position of seed found is not too close to the end of the string
        if (pos + test.length() > bwt->reference.length())
            mismatches = MISMATCHES_THRESHOLD + 1;
        else
            mismatches = number_of_mismatches(test, bwt->reference.substr(pos, test.length()));       

        for (int i = bl->start + 1; i <= bl->end; ++i) {
            int new_pos = bwt->get_location(i);
            if (new_pos + test.length() > bwt->reference.length())
                continue;

            int new_mismatches = number_of_mismatches(test, bwt->reference.substr(new_pos, test.length()));       
            if (new_mismatches < mismatches) {
                mismatches = new_mismatches;
                pos = new_pos;
            }
        }

        if (mismatches > MISMATCHES_THRESHOLD) {
            cout << "No matches found" << endl;
            return;
        }

        // calculate score of match 
        // for start do not use minimum edit distance but assume only mismatches can be made and no insertions or deletions
        cout << "Number of mismatches: " << mismatches << " for read with id: " << id << " at location: " << pos << endl;
    }
}

// Function to reverse a string and take its complement
// TODO make this faster!
string reverse_complement(const string& s)
{
    char *rev = new char[s.length() + 1];
    for (unsigned long i = 0; i < s.length(); ++i)
        rev[i] = map_complement[s[s.length() - 1 - i]];

    rev[s.length()] = '\0';
	/* string rev(s1.rbegin(), s1.rend()); */
	return string(rev);
}
