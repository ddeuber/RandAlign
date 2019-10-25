#include <iostream>
#include <vector>
#include <string.h>
#include <algorithm>
#include "bwt.h"
#include "read_genes.h"
#include "samfile.h"
#include "randalign.h"

using namespace std;

// takes as input a number that corresponds to a position to the current reference gene (after preprocessing the N's)
// and returns the relative position on the original string
// TODO This function could be done much faster using a binary search scheme.
// Is this really necessary? In the large genome file only around 10 regions of N's can be spotted
int convert_index_to_original_index(int index, vector< pair<int, int> > holes) {
    
    int new_index = index;

    for(vector< pair<int, int> >::size_type i = 0; i != holes.size(); i++) {
        if (index >= holes[i].first)
            new_index += holes[i].second;
        else
            break;
    }

    return new_index;
}


// char* getCmdOption(char ** begin, char ** end, const std::string & option)
// {
//     char ** itr = std::find(begin, end, option);
//     if (itr != end && ++itr != end)
//     {
//         return *itr;
//     }
//     return 0;
// }

// bool cmdOptionExists(char** begin, char** end, const std::string& option)
// {
//     return std::find(begin, end, option) != end;
// }


int main(int argc, char** argv) {

    if (argc < 2){
        cout << "No option provided" << endl;
    }

    string refName;
	int refLength;
    // vector that contains pairs of positions and lengths on sequences of N's
    vector< pair<int, int> > holes;
    string reference;
    BWT* bwt;

    /* available options
       index: create index, store and exit
            ./<appname> index <genome_file> 
       index_run: create index, store and run
            ./<appname> index_run <genome_file> <read_file_1> <read_file_2>
       recover: recover the index and run
            ./<appname> recover <read_file_1> <read_file_2>
       run: just run
            ./<appname> run <genome_file> <read_file_1> <read_file_2>
    */
    if ((strcmp(argv[1], "index") == 0 && argc != 3) || (strcmp(argv[1], "index_run") == 0&& argc != 5) \
     || (strcmp(argv[1], "recover") == 0 && argc != 4) || (strcmp(argv[1], "run") == 0&& argc != 5)){
         cout << "Invalid number of arguments" << endl;
         exit(1);
     }


    int allignment_file_index;

    if (strcmp(argv[1], "index") == 0 || strcmp(argv[1], "index_run") == 0){
        
        cout << "Creating index in directory: \"index\"" << endl;

        reference = read_reference_gene(argv[2], refName, refLength, holes);
        bwt = new BWT(reference, holes);

        bwt->store_index("index");

        if (strcmp(argv[1], "index") == 0){
            cout << "Created index and exiting" << endl;
            return 0;
        }

        // to get files with allignments skip first 2 words corresponding to the words "index" and <file_with_genome>
        allignment_file_index = 3;

    } else if (strcmp(argv[1], "recover") == 0){
        string directory;
        cout << "Recovering index from directory: \"index\"" << endl;

        bwt = new BWT("index");

        // to get files with allignments skip first 1 word corresponding to the words "recover"
        allignment_file_index = 2;
    } else if (strcmp(argv[1], "run") == 0){
        // in this case create the index and run
        reference = read_reference_gene(argv[2], refName, refLength, holes);
        bwt = new BWT(reference, holes);

        // to get files with allignments skip first 1 word corresponding to the words "recover"
        allignment_file_index = 3;
    } else {
        cout << "Unkownn option: " << argv[1] << ". Exiting.." << endl;
        exit(1);
    }

    cout << "Running allignment" << endl;


	string input1(argv[allignment_file_index]);
	string output = input1.substr(0, input1.length()-4);
	output += ".generated.mod.sam";
    
	SAMFile samFile(output, refName, refLength);
	RandomizedAligner randAlign(bwt, &samFile);
   
    ReadGenes* rg = new ReadGenes(argv[allignment_file_index], argv[allignment_file_index + 1]);

    string seed, reverse_complement_str;
    //block * bl;
    while (true) {
        read_block* rb = rg->get_one_read();
        if (rb == NULL)
            break;
            
        //cout << "---------------------------------------------" << endl;
        //cout << rb->id << endl;

		randAlign.align_and_print(rb);

		/*
        seed = rb->forward_read.substr(0, SEED_LENGTH);
        bl = bwt->get_matches(seed);
        print_matches(bl, rb->id, rb->forward_read, bwt);

        seed = rb->backward_read.substr(0, SEED_LENGTH);
        bl = bwt->get_matches(seed);
        print_matches(bl, rb->id, rb->backward_read, bwt);

        reverse_complement_str = reverse_complement(rb->forward_read);
        seed = reverse_complement_str.substr(0, SEED_LENGTH);
        bl = bwt->get_matches(seed);
        print_matches(bl, rb->id, reverse_complement_str, bwt);

        reverse_complement_str = reverse_complement(rb->backward_read);
        seed = reverse_complement_str.substr(0, SEED_LENGTH);
        bl = bwt->get_matches(seed);
        print_matches(bl, rb->id, reverse_complement_str, bwt);
		*/
    }

	samFile.close();
    return 0;
}
