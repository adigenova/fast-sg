#!/usr/bin/env perl
###############################################################################
# Author: Alex Di Genova 
# Laboratory: ERABLE/Mathomics
# Copyright (c)
# year: 2017
###############################################################################
use Getopt::Std;
use strict;


my %opts = ();
getopts( "k:l:r:t:x:y:z:u:v:w:s:c:p:N:R:", \%opts) or die usage();

#define and check required variables
my $kmer_size=0;
my $read_file="";
my $contig_file="";
my $prefix=$opts{p};
my $nthread=1;
## in c++ we ask for values larger than >, so defaults became 10,6,4
my $vhs_ill=9;	
my $ms_ill=5;	
my $msr_ill=3;	
#advanced options for long reads 20, 16, 4, 200
my $vhs_long=19;
my $ms_long=15;
my $msr_long=3;
my $ssr_long=200;

	
if ( !defined $opts{k} or !defined $opts{l} or !defined $opts{r} or !defined $opts{p}) {
	print "\n\nERROR: Mandatory options need to be specified -k -l -r -p\n\n";
	usage();
}

#check required variables
my @kmers=();
#check k-mer value
if ($opts{k} =~ /^\d+?$/ or $opts{k}=~/^\d+-\d+:\d+$/) {
    	# integer
	if($opts{k}=~/^(\d+)-(\d+):(\d+)$/){
		#we store the range of kmers
		for(my $i=$1;$i<=$2; $i+=$3){
			push(@kmers,$i);
		}
	}else{
		push(@kmers,$opts{k});
	}
	foreach my $kmer_size(@kmers){	
		if($kmer_size < 12){
	    		print "\n\nERROR: Value expecified for kmer_size (-k) is shorter than the minimum allowed [12]\n\n";
	     		usage();	
		}
	
		if($kmer_size > 256){
	     		print "\n\nERROR: Value expecified for kmer_size (-k) is larger than the maximum allowed [256]\n\n";
	     		usage();	
		}
	}
}else{
		
	     print "\n\nERROR: Value expecified for kmer_size (-k) is not an integer value\n\n";
	     usage();	
}

#check read configuration
check_read_file($opts{l});
$read_file=$opts{l};
#check file of contigs
if(-s $opts{r}){
}else{
	print "ERROR: contig file $opts{r} don't exist or has 0 size\n";
	usage();	
}

$contig_file=$opts{r};

#check the number of thread
if(defined $opts{t}){
	#check if integer	
	if ($opts{t} =~ /^\d+?$/) {
		$nthread=$opts{t};
	}else{
	     print "\n\nERROR: Value expecified  for number of thread (-t) is not an integer value\n\n";
	     usage();	
	}
}


#check advanced short variables if necesary
if(defined $opts{x} or defined $opts{y} or defined $opts{z}){
	
	if(defined $opts{x} and defined $opts{y} and defined $opts{z}){
		#we need to check that they are number
		if ($opts{x} =~ /^\d+?$/ and $opts{y} =~ /^\d+?$/  and $opts{z}=~ /^\d+?$/){
				if($opts{x} > $opts{y} and $opts{y} >= $opts{z}){
					#overwrite the defaults options
					$vhs_ill=$opts{x};	
					$ms_ill=$opts{y};	
					$msr_ill=$opts{z};	
				}else{
	     			print "\n\nERROR: Values expecified for advanced short read options must follow x>y>z\n\n";
	    	 		usage();	
				}
			
		}else{
	     		print "\n\nERROR: Values expecified for advanced short read options must be integers\n\n";
	    	 	usage();	
		}
		
	}else{
	     print "\n\nERROR: For advanced short read options you need to define x,y and z values\n\n";
	     usage();	
	}
}


#check advanced long variables if necesary
if(defined $opts{u} or defined $opts{v} or defined $opts{w} or defined $opts{s}){
	
	if(defined $opts{u} and defined $opts{v} and defined $opts{w} and defined $opts{s}){
		#we need to check that they are number
		if ($opts{u} =~ /^\d+?$/ and $opts{v} =~ /^\d+?$/  and $opts{w}=~ /^\d+?$/ and $opts{s}=~ /^\d+?$/){
				if($opts{u} > $opts{v} and $opts{v} >= $opts{w}){
					#overwrite the defaults options
					$vhs_long=$opts{u};	
					$ms_long=$opts{v};	
					$msr_long=$opts{w};	
					$ssr_long=$opts{s};
					
				}else{
	     			print "\n\nERROR: Values expecified for advanced long read options must follow x>y>z\n\n";
	    	 		usage();	
				}
			
		}else{
	     		print "\n\nERROR: Values expecified for advanced long read options must be integers\n\n";
	    	 	usage();	
		}
		
	}else{
	     print "\n\nERROR: For advanced long read options you need to define u,v,w and s values\n\n";
	     usage();	
	}
}

#We need to check for fast-sg exec and kmc
my $fhome=dirname($0);
my $fbin=$fhome."/FastSG++";#fastsg binary
my $khome=$fhome."/KMC";#kmc home
my $kbin=$khome."/bin/kmc";#kmc binary
my $kdump=$khome."/bin/kmc_dump";#kmc binary

#kmc binary
if(defined $opts{c}){
	$khome=$opts{c};	
	$kbin=$khome."/bin/kmc";#kmc binary
	$kdump=$khome."/bin/kmc_dump";#kmc binary
}
#we check if KMC binaries are found
#is a regular file and is executable
if(-f $kbin and -x $kbin and -f $kdump and -x $kdump){
	#print "KMC  binary: $kbin and $kdump\n";	
}else{
	if(!defined $opts{c}){
		print "ERROR: KMC binaries not found, it need to be in :\n$kbin\n$kdump\n";
		exit 1;	
	}else{
		print "ERROR: You defined KMC directory in $opts{c} and we cannot found the binaries there:\n$kbin\n$kdump\n";
		exit 1;
	}
}

#is a regular file and is executable
if(-f $fbin and -x $fbin){
	#print "FAST-SG binary: $fbin\n";	
}else{
	print "ERROR: FAST-SG binary not found, it need to be in :\n$fbin\n";
	exit 1;	
}

#FAST-SG advanced parameters
my $fvars=join(" ",$vhs_ill,$ms_ill,$msr_ill,$vhs_long,$ms_long,$msr_long,$ssr_long);
#create make file with K variable
open(FILE,">".$prefix.".mk") or die "cannot create $prefix.mk\n";
#build  makefile
print FILE ".DELETE_ON_ERROR:\n\n";
print FILE "#program locations\nKMCC=$kbin\n";
print FILE "KMCD=$kdump\n";
print FILE "FSG=$fbin\n\n";
#KMC commands
print FILE "$prefix.kmc.K\${K}.txt:\n";
print FILE "\t\${KMCC} -r -k\${K} -v -m4 -fm -ci1 -cx1 -t$nthread $contig_file $prefix.kmc.K\${K} $prefix.kmctmp.K\${K} 2>$prefix.kmc.K\${K}.err >$prefix.kmc.K\${K}.log\n";
print FILE "\t\${KMCD} $prefix.kmc.K\${K} $prefix.kmc.K\${K}.txt\n";
#FAST-SG commands
print FILE "$prefix.fast-sg.K\${K}.done: $prefix.kmc.K\${K}.txt\n";
print FILE "\t\${FSG} \$^ $contig_file \${K} $nthread $read_file $fvars 2>\$@.err \n";
print FILE "\ttouch $prefix.fast-sg.K\${K}.done\n";
#magic intructions to makefile
print FILE "\nall: $prefix.kmc.K\${K}.txt $prefix.fast-sg.K\${K}.done \n";
print FILE "\nclean:\n\trm -f $prefix.kmc.K\${K}.txt $prefix.fast-sg.K\${K}.done\n";

#show commands to execute:
if(defined $opts{N}){
	foreach my $k (@kmers){
		print "\n\nFAST-SG COMMANDS K=$k :\n";
		system("make -f $prefix.mk K=$k all -n") == 0 or die "make failed $?"; 
	}
	exit 0;
}

#clean results to restart
if(defined $opts{R}){
	foreach my $k (@kmers){
		print "\n\nCLEAN FAST-SG results K=$k\n\n";
		system("make -f $prefix.mk K=$k clean ") == 0 or die "make failed $?"; 
	}
	print "CLEAN done\n";
	exit 0;
}

#run FAST-SG 
foreach my $k (@kmers){
	print "\n\nSTARTING ".localtime." FAST-SG K=$k \n\n";
	system("make -f $prefix.mk K=$k all ") == 0 or die "make failed $?"; 
	print "\n\nENDING ".localtime." FAST-SG with K=$k\n\n";
}

close(FILE);	
print "all done\n";


#funtions
sub check_read_file{
	my ($file)=@_;
	open(FILE,$file) or die "cannot open read file : $file\n";
	#illumina reads
	#type prefix fwdfile revfile sam_output_format [1:unparied 2:paired] 
	#short lib1 saureus.reads_1.fq.gz saureus.reads_2.fq.gz 0
	#long reads
	#type prefix file insert_sizes sam output format [1:unparied 2:paired] 
	#long ultra_ont_raw ULTRA-LONG-RENAME-FOLD.fa.gz 2000,4000,6000,8000,10000 1
	my $prefix=();
	while(my $line=<FILE>){
		chomp $line;
		my @data=split (" ",$line);
		$prefix->{$data[1]}++;	#for checking the prefix at the end;
		if($data[0] eq "short"){
			if(scalar(@data) ne "5"){
				print "ERROR in read file: uncorrect number of columns we need 5 and there are ".scalar(@data)."\n$line\n";
				exit 1;
			}	
			#we check for the existance of file
			#foward
			if(-s $data[2]){
			}else{
				print "ERROR in read file: foward file don't exist or have 0 size. \n$data[2]\n";
				exit 1;
			}
			#reverse
			if(-s $data[3]){
			}else{
				print "ERROR in read file: rev file don't exist or have 0 size. \n$data[3]\n";
				exit 1;
			}
			#we check that the samflag is number
			if ($data[4] =~ /^\d+?$/) {
			   }else{
	     			print "ERROR in read file: Value expecified $data[4] for sam ouput is not an integer value[1 or 2]\n";
				exit 1;
			  }
			
		}elsif($data[0] eq "long"){
			if(scalar(@data) ne "5"){
				print "ERROR in read file: uncorrect number of columns we need 5 and there are ".scalar(@data)."\n$line\n";
				exit 1;
			}	
			#we check for the existance of file
			#long read file
			if(-s $data[2]){
			}else{
				print "ERROR in read file: long file don't exist or have 0 size. \n$data[2]\n";
				exit 1;
			}
			#reverse
			my @inserts=split(",",$data[3]);	
			foreach my $i (@inserts){
				if($i =~/^\d+?$/){
				}else{
	     			print "ERROR in read file: Value expecified $i for insert size is not an integer value\n";
				exit 1;
				}
			}
			#we check that the samflag is number
			if ($data[4] =~ /^\d+?$/) {
			   }else{
	     			print "ERROR in read file: Value expecified $data[4] for sam ouput is not an integer value[1 or 2]\n";
				exit 1;
			  }

		}else{
			print "ERROR in read file:$data[0] unknow read type[short or long]\n";
			exit 1;
		}
		
	}

	#we check that all prefix are unique
	foreach my $p (keys %{$prefix}){
		if($prefix->{$p} > 1){
			print "ERROR in read file: prefix $p used more than once\n";
			exit 1;
		}
	}

	close(FILE);
}

#get dirname
sub dirname {
	my $prog = shift;
	return '.' unless ($prog =~ /\//);
	$prog =~ s/\/[^\s\/]+$//g;
	return $prog;
}


sub usage {
  die(qq/

Usage: FAST-SG.pl  [<arguments>]\n

Basic options:
	
	-k	*K-mer size [12-256] could be a single value or a range [15-40:5]
	-l	*Read configuration file
	-r	*Contig fasta file
	-p	*Output prefix
	-t	Number of thread to be used [1]
	

Advanced options:

	Illumina reads options:
	-x	Vector hit size [10]	
	-y	Mimimum score to report reads [6]	
	-z	Mimimum score for pair rescue [4]	
	
	The value of -x must be larger than -y and -y larger than -z [-x > -y > -z]. 
	To modify these values take into account the k-mer size, the length and error rate of the reads.
	Specified values must be larger than defaults values.	
	
	Long reads options:
	-u	Vector hit size [20]	
	-v	Mimimum score to report a pair of reads [15]	
	-w	Mimimum score for pair rescue [4]	
	-s	Length of syntethic read [200]
	
	The value of -u must be larger than -v and -v larger than -w [-u > -v > -w]. 
	To modify these values take into account the k-mer size, the length and error rate of the reads.
	Specified values must be larger than defaults values.

Other options:
	-c		Full path to KMC3 directory [FAST-SG-DIR\/KMC]
	-N 		Show pipeline commands [-N 1]
	-R 		Delete results to restart FAST-SG [-R 1]
	--help 		Show help message
	--version 	Display current version of FAST-SG

Note: Options with description starting with (*) are mandatory.
\n/);

}

#current version and help messages
$Getopt::Std::STANDARD_HELP_VERSION = 1;
our $VERSION = "1.0"; 


#help 
sub HELP_MESSAGE {
    usage();
}


