Install tbb dependencies if not already installed (Ubuntu):
- sudo apt-get install libtbb-dev

To compile:
- make PARALLEL=0 # to run serial version
- make PARALLEL=1 # to run parallel version

To run:

index: create index, store and exit
    ./<appname> index <genome_file> .
index_run: create index, store and run
    ./<appname> index_run <genome_file> <read_file_1> <read_file_2>
recover: recover the index and run
    ./<appname> recover <read_file_1> <read_file_2>
run: just run
    ./<appname> run <genome_file> <read_file_1> <read_file_2>