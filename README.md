# Fast-SG
Fast-SG is an alignment-free algorithm for ultrafast scaffolding graph construction from short or long reads.

# Compilation Instructions

## Get Fast-SG code

	git clone https://github.com/adigenova/fastsg

## Get KMC (kmer couter)

### Mac users
Obtain pre-compiled binaries from:

	wget  https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.mac.tar.gz 

or compile yourself following the instructions provided in:

	https://github.com/refresh-bio/KMC

### Linux users	
Obtain pre-compiled binaries from:

	wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz

or compile yourself following the instructions provided in:

	https://github.com/refresh-bio/KMC

After getting KMC3, put the binaries inside Fast-SG directory, specifically in:

$Fast-SG/KMC/bin/kmc

$Fast-SG/KMC/bin/kmc_dump

## Compile Fast-SG
	make all

c++ compiler; compilation was tested with g++ version 5.3 (Linux) and clang version 4.2 (Mac OSX).

## Test Fast-SG (small test)
	make test
	
## Libraries used
Currently Fast-SG uses the following C++ opensource libraries:
 
1.- kseqcpp (https://github.com/ctSkennerton/kseqcpp.git checkout cfa50bcd17bbcb3225d431df4a2c1396f58a0993)

2.- ntHash (https://github.com/bcgsc/ntHash.git checkout ff326a8c9ccf6186f42c1f49950c1ebaadbd7f7a)

3.- BBHash (https://github.com/rizkg/BBHash.git checkout 99c905828a58fa119979df5c26bdbea93f0a7696)

4,- quasi_dictionary (https://github.com/pierrepeterlongo/quasi_dictionary.git checkout 9e8c64b150b129035f92d010a12085bd6c9490f0)

# Usage instructions
FAST-SG.pl is the wrapper script used to run FastSG++.
## Mandatory arguments
FAST-SG requires 4 mandatory arguments:

1.- The k-mer size (-k) restricted to the range [12-256]

2.- The set of contig sequences (-r) in FASTA format.

3.- The output prefix (-p)

4.- The read configuration file (-r) having the following format (space separated):

  Short reads:

	#type libID Path(fwd) Path(rev) SAM(1:single 2:paired)	

	short lib1 example/ecoli.ill-sim.fwd.fq.gz example/ecoli.ill-sim.rev.fq.gz 1

  Long reads:	

	#type libID Path(long read) Insert-sizes SAM(1:single 2:paired)	

	long pac example/ECOLI-PACSEQ.subset.fasta.gz 1000,2000,3000,5000 1

	long ont example/ECOLI-ONT-1D.subset.fasta.gz 1000,2000,3000,5000 1

Example of a read configuration file can be found in examples/ecoli-reads.txt

## Running Fast-SG (example):

	./FAST-SG.pl  -k 15 -l example/ecoli-reads.txt -r example/ecoli-illumina.fa.gz -p test
	
	To obtain help and advanced options details execute:

	./FAST-SG.pl  --help
	
