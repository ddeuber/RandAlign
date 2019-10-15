#include <iostream>
#include <vector>

#include "bwt.h"
#include "read_genes.h"

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


int main(int argc, char** argv) {

    if (argc != 4) {
        cout << "Usage: " << argv[0] << " <reference_gene_file> <reads_gene_file_1> <reads_gene_file      2>" << endl;
        exit(1);
    }

    // vector that contains pairs of positions and lengths on sequences of N's
    vector< pair<int, int> > holes;
    string reference = read_reference_gene(argv[1], holes);

    for(vector< pair<int, int> >::size_type i = 0; i != holes.size(); i++) {
        cout << holes[i].first << " " << holes[i].second << endl;
    }

    cout << "The total length of the string is " << reference.length() << endl;    

    cout << "0 " << convert_index_to_original_index(0, holes) << endl;
    cout << "647850 " << convert_index_to_original_index(647850, holes) << endl;
    cout << "2978253 " << convert_index_to_original_index(2978253, holes) << endl;
    cout << "2978255 " << convert_index_to_original_index(2978255, holes) << endl;


    exit(1);








    BWT* bwt = new BWT(reference, false);
   
    ReadGenes* rg = new ReadGenes(argv[2], argv[3]);

    string seed, reverse_complement_str;
    block * bl;
    while (true) {
        read_block* rb = rg->get_one_read();
        if (rb == NULL)
            break;
            
        cout << "---------------------------------------------" << endl;
        cout << rb->id << endl;

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
    }

    return 0;
}