#ifndef READ_GENES_H
#define READ_GENES_H 

#include <string>
#include <fstream>

#include "bwt.h"

#define ERROR_EXIT_CODE 1

typedef struct {
    std::string id;
    std::string forward_read;
    std::string forward_read_quality;
    std::string backward_read;
    std::string backward_read_quality;
} read_block;

// reafactor these
std::string reverse_complement(std::string const &s);
void print_matches(block* bl, std::string id, std::string test, BWT* bwt);

std::string read_reference_gene(std::string filename);

class ReadGenes {
    private:
        std::string filename_1;
        std::string filename_2;
        std::ifstream file_1;
        std::ifstream file_2;
    public:
        ReadGenes(std::string filename_1, std::string filename_2);
        read_block* get_one_read();       
};

#endif
