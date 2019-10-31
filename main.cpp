#include <iostream>
#include <vector>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include "bwt.h"
#include "read_genes.h"
#include "samfile.h"
#include "randalign.h"

#if PARALLEL_RUN
#include "tbb/parallel_for.h"
#endif

using namespace std;

// #define PARALLEL_READS 4096

int main(int argc, char** argv) {

    if (argc < 2){
        cout << "No option provided" << endl;
    }

    string refName;
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

        reference = read_reference_gene(argv[2], refName, holes);
        bwt = new BWT(reference, refName, holes);

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
        reference = read_reference_gene(argv[2], refName, holes);
        bwt = new BWT(reference, refName, holes);

        // to get files with allignments skip first 1 word corresponding to the words "recover"
        allignment_file_index = 3;
    } else {
        // should never reach this
        cout << "Unkownn option: " << argv[1] << ". Exiting.." << endl;
        exit(1);
    }

    cout << "Running allignment" << endl;


	string input1(argv[allignment_file_index]);
	string output = input1.substr(0, input1.length()-4);
	output += ".generated.mod.sam";
    cout << "Storing results at " << output << endl; 

    int refLength = bwt->reference.length() - 1;
	for (auto h : bwt->holes)
		refLength += h.second;
		
	SAMFile samFile(output, bwt->refName, refLength);
	RandomizedAligner randAlign(bwt, &samFile);
   
    ReadGenes* rg = new ReadGenes(argv[allignment_file_index], argv[allignment_file_index + 1]);

    unsigned real_count; 
    bool finished = false; 
    string seed, reverse_complement_str;
 
    //Read enviromental variable 
    int par_reads;
    char *v = getenv("PARALLEL_READS");
    
    if(v == NULL){
        par_reads=4096; //set default chunk size
    } else{
        par_reads=atoi(v); //convert string to integer
    }

    cout << "Total reads performed in parallel: " << par_reads << endl; 

    read_block * read_array[par_reads]; 
    results_block * res_block[par_reads];

    while (!finished) {
        // Initialize
        real_count = 0u; 
        for(int i=0; i<par_reads; i++){
            read_block* rb = rg->get_one_read();
            if (rb == NULL){
                finished = true; 
                break; 
            }
            else {
                read_array[i] = rb; 
                ++real_count; 
            }
        }

#if PARALLEL_RUN
        /* Parallel Execution */
        tbb::parallel_for(0u, (unsigned)real_count, [&](unsigned i){
            res_block[i] = randAlign.align_and_print(read_array[i]);
        });
#else 
        for(int i = 0; i < real_count; i++){
            res_block[i] = randAlign.align_and_print(read_array[i]);
        }
#endif
		
        /* Print to SAMFle */
        for(int i=0; i<real_count; i++){
            if (res_block[i]->failed){
                std::cout << "Failed for read " << read_array[i]->id << std::endl;
                std::cout << std::endl;
            }
            samFile.add_paired_read_entry(
                        read_array[i]->id, 
                        res_block[i]->readSeq1, 
                        res_block[i]->qualSeq1, 
                        res_block[i]->pos1, 
                        res_block[i]->cigar1, 
                        res_block[i]->readSeq2, 
                        res_block[i]->qualSeq2, 
                        res_block[i]->pos2, 
                        res_block[i]->cigar2, 
                        res_block[i]->read1Reversed);
        }
    }

	samFile.close();
    return 0;
}
