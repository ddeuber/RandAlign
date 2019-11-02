# Efficient Search and Read Alignment

This project is part of the Computational Biomedicine I (Fall Semester 2019) at ETH.


| Name  | Email |
| ------------- | ------------- |
| Deuber David  | ddeuber@student.ethz.ch  |
| Adamos Solomou  | solomoua@student.ethz.ch  |
| Anagnostidis Sotiris  | sanagnos@student.ethz.ch  |
| Vasilakopoulos George  | gvasilak@student.ethz.ch  |

## Project structure

    .
    ├── SuffixTree.cpp                      # implementation of Ukkonen's algorithm to build the suffix tree
    ├── bwt.cpp                             # create the index used for the alignment
    ├── read_genes.cpp                      # read from the files containing the reference and the read genes
    ├── randalign.cpp                       # implementation of the alignment strategy
    ├── samfile.cpp                         # generates the outfile in the required form
    ├── main.cpp                            # main function
    └── README.md


## Getting Started

Install tbb dependencies if not already installed (Ubuntu):
 ```
 sudo apt-get install libtbb-dev
 ```

To compile:
 ```
 make PARALLEL=0 # to run serial version
 make PARALLEL=1 # to run parallel version
 ```

 Default compiler `clang`. Can be changed by assigning `CXX` during compiling.

To run:

 ```
 ./<appname> index <genome_file> # index: create index, store and exit
 ./<appname> index_run <genome_file> <read_file_1> <read_file_2> # index_run: create index, store and run
 ./<appname> recover <read_file_1> <read_file_2> # recover: recover the index and run
 ./<appname> run <genome_file> <read_file_1> <read_file_2> #run: just run alignments
 ```

We are assuming the read_files to have the form `<read_file>Cov1.fq` and `<read_file>Cov2.fq` <br>
The generated sam output file will be located in the same directory under the naming `<read_file>Cov.generated.mod.sam`


To evaluate our method we compared results generated from other read alignment tools (bowtie2), outputing in 100% of the cases the same alignments for the small genome and more than 98% of the cases the same alignments for the large genome.

There is also the option to modify the environmental varialbe `PARALLEL_READS` to fine tune the number of possible parallel reads, although the default value gives the best results.
