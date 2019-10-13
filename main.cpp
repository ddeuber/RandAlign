#include <iostream>

#include "bwt.h"
#include "read_genes.h"
#include "samfile.h"
#include "randalign.h"


using namespace std;

int main(int argc, char** argv) {

    if (argc != 4) {
        cout << "Usage: " << argv[0] << "<reference_gene_file> <reads_gene_file_1> <reads_gene_file_2>" << endl;
        exit(1);
    }

	string refName;
    string reference = read_reference_gene(argv[1], refName);

    BWT* bwt = new BWT(reference, false);

	string input1(argv[2]);
	string output = input1.substr(0, input1.length()-4);
	output += ".generated.mod.sam";
	SAMFile samFile(output, refName, reference.length()-1);
	RandomizedAligner randAlign(bwt, &samFile);
   
    ReadGenes* rg = new ReadGenes(argv[2], argv[3]);

    string seed, reverse_complement_str;
    //block * bl;
    while (true) {
        read_block* rb = rg->get_one_read();
        if (rb == NULL)
            break;
            
        cout << "---------------------------------------------" << endl;
        cout << rb->id << endl;

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
