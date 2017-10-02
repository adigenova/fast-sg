
CXX ?= g++
CFLAGS = -O3 -std=c++11 -lpthread -lz
EXEC= FastSG++

ifeq ($(prof),1)
 CFLAGS+= -pg
endif
ifeq ($(deb),1)
 CFLAGS+= -O0 -DASSERTS -g
endif

.PHONY: test clean


all: $(EXEC)

FastSG++: FastSG.cpp modules
	$(CXX) -o $@  $< $(CFLAGS)

%.o: %.cpp %.h
	$(CXX) -o $@ -c $< $(CFLAGS)
clean:
	-rm -f  FastSG++ modules *.sam
#run small Fast-SG test and compare the resuting samfiles 
test: $(EXEC)
	./FAST-SG.pl  -k 15 -l example/ecoli-reads.txt -r example/ecoli-illumina.fa.gz -p test
	ls *.sam | xargs -n 1 sh  misc/compare-results.sh

modules :
	git clone https://github.com/ctSkennerton/kseqcpp.git && cd kseqcpp  && git checkout cfa50bcd17bbcb3225d431df4a2c1396f58a0993
	git clone https://github.com/bcgsc/ntHash.git && cd ntHash  && git checkout ff326a8c9ccf6186f42c1f49950c1ebaadbd7f7a
	git clone https://github.com/rizkg/BBHash.git && cd BBHash  && git checkout 99c905828a58fa119979df5c26bdbea93f0a7696
	git clone https://github.com/pierrepeterlongo/quasi_dictionary.git && cd quasi_dictionary  && git checkout 9e8c64b150b129035f92d010a12085bd6c9490f0
	touch $@

##code to download kmc in authomatic way
##kmc :
##ifeq ($(shell uname),Darwin)
##	curl https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.mac.tar.gz > KMC3.mac.tar.gz
## Run MacOS commands 
##else
## check for Linux and run other commands
##	wget https://github.com/refresh-bio/KMC/releases/download/v3.0.0/KMC3.linux.tar.gz	
##endif

