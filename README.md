Install tbb dependencies if not already installed:
- git clone https://github.com/wjakob/tbb.git
- cd tbb/build
- cmake ..
- make -j
- sudo make install
- LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib
- export LD_LIBRARY_PATH

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