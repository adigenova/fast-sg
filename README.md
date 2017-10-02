# FAST-SG
FAST-SG is an alignment-free algorithm for ultrafast scaffolding graph construction from short or long reads.

#Compilation Instructions

## Get FAST-SG code

	git clone https://github.com/adigenova/fastsg

## Get KMC (kmer couter)

### Mac users
Obtain pre-compiled binaries from:

	wget  https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.mac.tar.gz 

or compile yourself following the instruction provided in:

	https://github.com/refresh-bio/KMC

### Linux users	
Obtain pre-compiled binaries from:

	wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz

or compile yourself following the instruction provided in:

	https://github.com/refresh-bio/KMC

After download or compile KMC3, put the binaries inside FAST-SG directory, specifically in :

$FAST-SG/KMC/bin/kmc

$FAST-SG/KMC/bin/kmc_dump

## Compile FAST-SG
	make all
## Test FAST-SG (small test)
	make test

# Libraries required
Currently FAST-SG use the following C++ opensource libraries:
 
1.- kseqcpp (https://github.com/ctSkennerton/kseqcpp.git checkout cfa50bcd17bbcb3225d431df4a2c1396f58a0993)

2.- ntHash (https://github.com/bcgsc/ntHash.git checkout ff326a8c9ccf6186f42c1f49950c1ebaadbd7f7a)

3.- BBHash (https://github.com/rizkg/BBHash.git checkout 99c905828a58fa119979df5c26bdbea93f0a7696)

4,- quasi_dictionary (https://github.com/pierrepeterlongo/quasi_dictionary.git checkout 9e8c64b150b129035f92d010a12085bd6c9490f0)



