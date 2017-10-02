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
FAST-SG  use 

