//
// Created by Alex Di Genova on 16/10/2016.
//
//external libraries
#include "kseqcpp/kseq.hpp"
#include "BBHash/BooPHF.h"
#include "ntHash/nthash.hpp"
#include "quasi_dictionary/src/probabilistic_set.h"


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/types.h>
#include <random>
#include <algorithm>
#include <fstream>
#include <fcntl.h>
#include <iostream>
#include <map>
#include <sstream>
#include <iterator>
#include <thread>



using namespace std;

//Custom fake hasher for Boophf because the values are already hashed into int64_t using nthash
class Custom_uint64_Hasher
{
public:
    uint64_t operator ()   (uint64_t key, uint64_t seed=0) const
    {
        return key;
    }
};

//we define the perfect hash funtion datatypes
typedef boomphf::mphf<  u_int64_t, Custom_uint64_Hasher  > boophf_t;
//to work in file instead of memory, code from Boophf
// iterator from disk file of u_int64_t with buffered read,   todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, u_int64_t>{
public:

    bfile_iterator()
            : _is(nullptr)
            , _pos(0) ,_inbuff (0), _cptread(0)
    {
        _buffsize = 10000;
        _buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
    }

    bfile_iterator(const bfile_iterator& cr)
    {
        _buffsize = cr._buffsize;
        _pos = cr._pos;
        _is = cr._is;
        _buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
        memcpy(_buffer,cr._buffer,_buffsize*sizeof(u_int64_t) );
        _inbuff = cr._inbuff;
        _cptread = cr._cptread;
        _elem = cr._elem;
    }

    bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
    {
        _buffsize = 10000;
        _buffer = (u_int64_t *) malloc(_buffsize*sizeof(u_int64_t));
        int reso = fseek(_is,0,SEEK_SET);
        advance();
    }

    ~bfile_iterator()
    {
        if(_buffer!=NULL)
            free(_buffer);
    }


    u_int64_t const& operator*()  {  return _elem;  }

    bfile_iterator& operator++()
    {
        advance();
        return *this;
    }

    friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
    {
        if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
        assert(lhs._is == rhs._is);
        return rhs._pos == lhs._pos;
    }

    friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
private:
    void advance()
    {
        _pos++;

        if(_cptread >= _inbuff)
        {
            int res = fread(_buffer,sizeof(u_int64_t),_buffsize,_is);
            _inbuff = res; _cptread = 0;

            if(res == 0)
            {
                _is = nullptr;
                _pos = 0;
                return;
            }
        }

        _elem = _buffer[_cptread];
        _cptread ++;
    }
    u_int64_t _elem;
    FILE * _is;
    unsigned long _pos;

    u_int64_t * _buffer; // for buffered read
    int _inbuff, _cptread;
    int _buffsize;
};

//binary file to store uint64_t in file
class file_binary{
public:

    file_binary(const char* filename)
    {
        _is = fopen(filename, "rb");
        if (!_is) {
            throw std::invalid_argument("Error opening " + std::string(filename));
        }
    }

    ~file_binary()
    {
        fclose(_is);
    }

    bfile_iterator begin() const
    {
        return bfile_iterator(_is);
    }

    bfile_iterator end() const {return bfile_iterator(); }

    size_t        size () const  {  return 0;  }//todo ?

private:
    FILE * _is;
};


//global array for  to compute revcomp
char comp_tab[] = {
        0,   1,       2,       3,       4,   5,       6,       7,       8,   9,  10,  11,      12,  13,  14,  15,
        16,  17,  18,  19,      20,  21,  22,  23,      24,  25,  26,  27,      28,  29,  30,  31,
        32,  33,  34,  35,      36,  37,  38,  39,      40,  41,  42,  43,      44,  45,  46,  47,
        48,  49,  50,  51,      52,  53,  54,  55,      56,  57,  58,  59,      60,  61,  62,  63,
        64, 'T', 'V', 'G', 'H', 'E', 'F', 'C', 'D', 'I', 'J', 'M', 'L', 'K', 'N', 'O',
        'P', 'Q', 'Y', 'S', 'A', 'A', 'B', 'W', 'X', 'R', 'Z',  91,      92,  93,  94,  95,
        64, 't', 'v', 'g', 'h', 'e', 'f', 'c', 'd', 'i', 'j', 'm', 'l', 'k', 'n', 'o',
        'p', 'q', 'y', 's', 'a', 'a', 'b', 'w', 'x', 'r', 'z', 123, 124, 125, 126, 127
};



/*structure to store the results of lookup in the PHF*/
struct hit{
    bool map;
    uint ctg;
    uint pos;
    bool strand;
    string kseqm;//kmer matching
    string seqname;
    int score;
    //constructor of the hit structure
    hit(){
        map=0;
        ctg=0;
        pos=0;
        strand=0;
        kseqm="";
        seqname="";
        score=0;
    }
};


//class declaration
class KmerDB;
class Illumina;
class LongReads;

//threads varibales
struct parallelref {

    pthread_mutex_t mut;
    pthread_mutex_t pin;
    vector<string> *pfiles;
    KmerDB *db;
    uint64_t *counter;
};
//thread function declaration
void * threaded_ref (void* args);//map in parallel bases from the reference
void * threaded_shortmap2 (void* args);//map in parallel short illumina reads
void * threaded_longmap (void* args);//map in parallel long PacBio/ONT reads
//Kmer database store Kmers in a MPHF and contigs coordinates for each key

class KmerDB
{
private :
    uint KmerSize;
    u_int64_t numberOfKmers;
    int bytesPerKmer;
    uint8_t* Kmers;//arrays of Kmers for exact mapping
    //reference related variables
    map< uint , string > id2seq;
    map< string , uint > seq2id;
    map< uint , uint > id2len;
    vector <uint> pos;
    vector <uint> ctg;
    vector <bool> strand;
    //MPHF related variables
    hash<string> hash_fn;//hash function for strings to numbers
    boophf_t * bphf ;//perfect hash function to ask in constant time to the array of Kmers
    string infile;
    probabilisticSet *ps;//array that store a finger print of the kmer for probabilistic search with FP rate of ~ 1/k*2

public:

    pthread_mutex_t mutex;//to access to read files
    pthread_mutex_t pinfo;//to print information from the threads

   KmerDB(string kmerfile,uint kmerSize,string rfasta,uint threads){
        //we read the file of kmers
        KmerSize=kmerSize;
       //allocate the array
        ifstream readFile(kmerfile);
        string line;
        uint64_t nelem = 0;
       FILE * key_file = NULL;
       printf("Hashing uniq kmers...\n");
       double t_begin,t_end; struct timeval timet;
       gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
       //to obtain a unique tmp number for the uniq file that we use
       infile=kmerfile+"."+to_string(getpid())+".hash64";
       key_file = fopen(infile.c_str(),"w+");
       u_int64_t current = 0;
       while(getline(readFile,line)){
            istringstream iss(line);
            vector<string> tokens ( (istream_iterator<string>(iss)),istream_iterator< string > ());
            current=getFhval(tokens[0].c_str(),kmerSize);;//we hash the current kmer
            fwrite(&current, sizeof(u_int64_t), 1, key_file);
            pos.push_back(0);//we made the initialization of associated variables
            ctg.push_back(0);//
            strand.push_back(0);//
            nelem++;
       }
       readFile.close();
       //we close and then we remove the key file
       fclose(key_file);
       //cout << "File has " << nelem<<" kmers"<<endl;
       gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
       double elapsed = t_end - t_begin;
       numberOfKmers=nelem;
       printf("time to hash %llu kmers in %.2fs\n", numberOfKmers,elapsed);
       //create the probabilistic structure with an error rate of 1/2^16 = 0,00001525878906 or 99,9999847412 effectiveneess
       ps=new probabilisticSet(numberOfKmers,16);
       //we build the perfect hash funtion
       build_mphf(threads);
       //we load the coordiantes for each kmer stored in the MPHF
       load_reference_data(rfasta,threads);
    }

    uint64_t get_number_kmers(void){
        return numberOfKmers;
    }

    void build_mphf(uint number_cpu){
        double t_begin,t_end; struct timeval timet;
        bphf = NULL;
        printf("Construct a BooPHF with  %lli elements  \n",numberOfKmers);
        gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        double gammaFactor = 4.0; // lowest bit/elem is achieved with gamma=1, higher values lead to larger mphf but faster construction/query
        //build the mphf from file
        auto data_iterator = file_binary(infile.c_str());
        bphf = new boomphf::mphf<u_int64_t,Custom_uint64_Hasher>(numberOfKmers,data_iterator,number_cpu,gammaFactor);
        //we store a fingerprint for each key in the database
        for (auto it = data_iterator.begin(); it != data_iterator.end(); ++it ){
            uint64_t key=*(it);
            //cout << *(it) <<endl;
            uint64_t idx = bphf->lookup(key);
            if(idx > numberOfKmers){
                cout <<"Error building MPHF function key: "<<key<<" having an index larger than total number of kmers"<<endl;
                exit(1);
            }else{
                //add a finger print to the probabilistic set
                ps->add(idx,key);
            }
        }
        //we remove temporary file
        if( remove(infile.c_str()) != 0 ) {
            cout << "Error deleting file "<< infile <<endl;
        }
        gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
        double elapsed = t_end - t_begin;
        printf("BooPHF constructed perfect hash for %llu keys in %.2fs\n", numberOfKmers,elapsed);
        printf("boophf  bits/elem : %f\n",(float) (bphf->totalBitSize())/numberOfKmers);
    }

    vector<uint64_t> mapkmer(string k){
        vector<uint64_t>map;
        uint64_t  ifwd = bphf->lookup(getFhval(k.c_str(),KmerSize));
        if(ifwd < numberOfKmers){
            //if(compareKmer(k,hdata[ifwd])){
            if(ps->exists(ifwd,getFhval(k.c_str(),KmerSize))){
                map.push_back(ifwd);//kmer id
                map.push_back(0);//kmers corresponde to foward
                return map;
            }
        }
        //we try the reverse kmer
        string rk=revcomp(k);
        uint64_t irev = bphf->lookup(getFhval(rk.c_str(),KmerSize));
        if(irev < numberOfKmers){
            //if(compareKmer(rk,hdata[irev])){
            if(ps->exists(irev,getFhval(rk.c_str(),KmerSize))){
                map.push_back(irev);//kmer id
                map.push_back(1);//kmers correspond to reverse
                return map;
            }
        }
        //no mapping in fwd or rev kmer dont belong to the datase
        map.push_back(ULLONG_MAX);
        return map;
    }
    //we ask if a value exist in the database
    vector<uint64_t> mapkmer(uint64_t hfwd, uint64_t hrev){
        vector<uint64_t> map;
        uint64_t  ifwd = bphf->lookup(hfwd);
        if(ifwd < numberOfKmers){
            if(ps->exists(ifwd,hfwd)){
                map.push_back(ifwd);//kmer id
                map.push_back(0);//kmers corresponde to foward
                return map;
            }
        }
        //we try the reverse kmer
        uint64_t irev = bphf->lookup(hrev);
        if(irev < numberOfKmers){
            if(ps->exists(irev,hrev)){
                map.push_back(irev);//kmer id
                map.push_back(1);//kmers correspond to reverse
                return map;
            }
        }
        //no mapping in fwd or rev kmer dont belong to the datase
        map.push_back(ULLONG_MAX);
        return map;
    }

    //heng li revcomp this is for printing only
    string revcomp(string kmer){
        int c0, c1;
        int klen=kmer.length();
        for (int i = 0; i < klen>>1; ++i) { // reverse complement sequence
            c0 = comp_tab[(int)kmer[i]];
            c1 = comp_tab[(int)kmer[klen - 1 - i]];
            kmer[i] = c1;
            kmer[klen - 1 - i] = c0;
        }
        if (klen & 1) // complement the remaining base
            kmer[klen>>1] = comp_tab[(int)kmer[klen>>1]];
        return kmer;
    }

    //I need to paralize this step because with big genomes take a lot of time
    void load_reference_data(string fastafile, uint ntreads){
        //to reads the fasta file we expect the fastafile uncompressed
        kseq seq;
        gzFile fp1 = gzopen(fastafile.c_str(),"r");
        FunctorZlib r1;
        kstream<gzFile, FunctorZlib> ks1(fp1, r1);
        int l=0;
        uint seqid=0;
        double t_begin,t_end; struct timeval timet;
        uint64_t query=0;
        gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        uint64_t length=0;
        ofstream fastq;
        vector<string> pfiles;//vector that save the name of fastq files;
        int part=0;
        string pid=to_string(::getpid());
        fastq.open("refseq"+to_string(part)+"."+pid+".txt");
        pfiles.push_back("refseq"+to_string(part)+"."+pid+".txt");
        while((l = ks1.read(seq)) >= 0) {
            if(seq.seq.length() < KmerSize){
                continue;
            }
            if(length > 10000000){
                length=0;
                fastq.close();
                part++;
                fastq.open("refseq"+to_string(part)+"."+pid+".txt");
                pfiles.push_back("refseq"+to_string(part)+"."+pid+".txt");
            }
            //we write to the file and then we store the variables
            fastq << seqid <<" "<<seq.seq<<endl;
            id2seq[seqid]=seq.name;//store  contig name
            id2len[seqid]=seq.seq.length();//store contig length
            seq2id[seq.name]=seqid;
            seqid++;
            length+=seq.seq.length();
        }

        gzclose(fp1);
        //multithreading version of database building
        //we create the threads for concurrency in the alignments step
        pthread_t *tab_threads= new pthread_t [ntreads];
        //we create the mutex to reads sequence from fasta file in an atomic way
        pthread_mutex_init(&mutex, NULL);
        pthread_mutex_init(&pinfo, NULL);
        parallelref tmp;
        tmp.mut=mutex;//MUTEX to control access to files
        tmp.pin=pinfo;//MUTEX to control the ouput fo files
        tmp.pfiles=&pfiles;
        //tmp.ks1=&ks2;
        tmp.db=this;//pointer to MPHF
        tmp.counter=&query;
        //create the threads
        for(int ii=0;ii<ntreads;ii++) {
            pthread_create(&tab_threads[ii], NULL, threaded_ref, &tmp);
        }
        //wait for thread to finish
        for(int ii=0;ii<ntreads;ii++)
        {
            pthread_join(tab_threads[ii], NULL);
        }

        gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
        double elapsed = t_end - t_begin;
        printf("Total time population database %.2fs\n",elapsed);
        printf("A total of %llu kmers were indexed \n", query);
    }

    uint getKmerSize(void){
        return KmerSize;
    }
    //return the position contig strand of the hit
    hit lookup_kmer(string k){
        vector< uint64_t > map=mapkmer(k);
        hit h;
        if(map[0] < numberOfKmers){
            h.map=1;
            h.ctg=ctg[map[0]];//contig position
            h.pos=pos[map[0]];//mapping position
            h.strand=0;//means foward
            h.kseqm=k;
            if(map[1] != strand[map[0]]){
                h.strand=1;//means reverse strand
                h.kseqm=revcomp(k);//we return the mapped kmer in reverse
            }
            //due to the probabilistic structure the pos can be 0
            if(h.pos == 0){
                h.map=0;
            }
        }
        return h;
    }
    //return the position in contig strands of the hit
    hit lookup_kmer(uint64_t hfwd, uint64_t hrev,string k){
        vector< uint64_t > map=mapkmer(hfwd,hrev);
        hit h;
        if(map[0] < numberOfKmers){
            h.map=1;
            h.ctg=ctg[map[0]];//contig position
            h.pos=pos[map[0]];//mapping position
            h.strand=0;//means foward
            h.kseqm=k;
            if(map[1] != strand[map[0]]){
                h.strand=1;//means reverse strand
                h.kseqm=revcomp(k);//we return the mapped kmer in reverse
                //we now that the kmer was mapped in reverse but we put in the sam file the kmer in foward, just to wing a little in speed
                //h.kseqm=k;
            }
            //due to the probabilistic structure the pos can be 0 we says that the kmer is unmap
            /*if(h.pos == 0){
                h.map=0; reduce the speed a lot
            }*/
        }
        return h;
    }
    //return the position in contig strands of the hit
    hit lookup_kmer(uint64_t hfwd, uint64_t hrev){
        vector< uint64_t > map=mapkmer(hfwd,hrev);
        hit h;
        if(map[0] < numberOfKmers){
            h.map=1;
            h.ctg=ctg[map[0]];//contig position
            h.pos=pos[map[0]];//mapping position
            h.strand=0;//means foward
            h.kseqm="";
            if(map[1] != strand[map[0]]){
                h.strand=1;//means reverse strand
                h.kseqm="";//we return the mapped kmer in reverse
            }
            //due to the probabilistic structure the pos can be 0 we says that the kmer is unmap
            if(h.pos == 0){
                h.map=0;
            }
        }
        return h;
    }
    //return the number of sequences stored in the database
   uint get_number_seq(void){
        return id2seq.size();
    }
    //return the name of the sequence from the id
   string get_seq_name(uint id){
       return id2seq[id];
   }
    //return the length of the sequence from the id
   uint get_seq_len(uint id){
       return id2len[id];
   }
    uint get_seq_id(string name){
        return seq2id[name];
    }
    //set references values
    void set_references_values(uint ctgid, uint posctg, bool strandk, uint64_t index){
            ctg[index]=ctgid;
            pos[index]=posctg;
            strand[index]=strandk;
    }

    void dump_database(void){
        for (int i = 0; i < ctg.size() ; i++) {
            cerr << i<<" "<<ctg[i]<<" "<<pos[i]<<" "<<strand[i]<<endl;
        }
    }
};


//we defined the parallel funtion to ask for reference basess in parallel
//todo:add variable to count how many uniq kmers were mapped
void * threaded_ref (void* args)
{
    //casting of variables
    parallelref *tmp = (parallelref*) args;
    pthread_mutex_t *mutex = &tmp->mut;
    pthread_mutex_t *pinfo = &tmp->pin;
    KmerDB *dbt=(KmerDB *)tmp->db;
    //kstream<int, FunctorRead> *ks1= (kstream<int, FunctorRead> *) tmp->ks1;
    uint64_t *gcounter=(uint64_t *)tmp->counter;
    vector <string> *files = ( vector<string> *) tmp->pfiles;
    //we start the pthred execution
    bool dojob=1;//variable to control thread execution
    //we check the time spent
    uint kmer_size;
    kmer_size = dbt->getKmerSize();
    uint64_t same=0;
    uint64_t query=0;
    double t_begin,t_end; struct timeval timet;
    gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
    string file_to_process="";
    while(dojob) {
        //we try to get a file to process
        pthread_mutex_lock(mutex);
        if(files->size() > 0) {
            file_to_process = files->back();
            files->pop_back();
        }else{
            dojob=0;
        }
        pthread_mutex_unlock(mutex);


        if(dojob) {

            ifstream readFile2(file_to_process);
            string line="";
            while(getline(readFile2,line)) {
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                //we iterate the stored bases
                int max_i = tokens[1].length() - kmer_size;
                uint seqid = (uint)atoi(tokens[0].c_str());
                //string seq=tokens[1];
                uint64_t hfwd = 0, hrev = 0;
                //we get both kmers in foward and reverse mers
                uint64_t hVal = 0;
                //int max_i= seq.seq.length() - KmerSize;
                hVal = NTPC64(tokens[1].substr(0, kmer_size).c_str(), kmer_size, hfwd, hrev); // initial hash value
                for (uint i = 0; i <= max_i; i++) {
                    vector<uint64_t> idx = dbt->mapkmer(hfwd, hrev);
                    if (idx[0] < ULLONG_MAX) {
                        dbt->set_references_values(seqid, i + 1, bool(idx[1]), idx[0]);
                        same++;
                    }
                    //get next kmer hash
                    hVal = NTPC64(tokens[1][i], tokens[1][i + kmer_size], kmer_size, hfwd, hrev); //recursive hash k+L
                    query++;
                }
            }
            //we remove the tmp file
            if( remove(file_to_process.c_str()) != 0 ) {
                pthread_mutex_lock(pinfo);
                cout << "Error deleting file "<< file_to_process <<endl;
                pthread_mutex_unlock(pinfo);
            }

        }
    }
    gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
    double elapsed3 = t_end - t_begin;
    //we print some stats from thread information
    pthread_mutex_lock(pinfo);
    cout << "End of mapping in thread\n";
    printf("Query of %llu kmer  in %.2fs rate %.2f per second\n", query,elapsed3,query/elapsed3);
    //adding variable to count the number of matched kmers
    (*gcounter)+=same;
    pthread_mutex_unlock(pinfo);
    return NULL;
}


//thread related variables
struct targ{
    pthread_mutex_t mut;
    pthread_mutex_t pin;
    pthread_mutex_t *pinm;
    //to control samout
    ofstream *sam;
    uint *samflag;
    uint *samp;
    //file pointer to fastq/fasta file for reading
    kstream<gzFile, FunctorZlib> *ks1;
    kstream<gzFile, FunctorZlib> *ks2;
    //parallel files names
    vector<string> *pfiles;
    //to control access to kmer database
    KmerDB* db;
    //to control access to illumina object
    Illumina* ill;
    //to control access to longRead object
    LongReads* lon;

};



//funtion to sort a pair of kmers in a vector
bool compare_by_ctg_pos(const hit& ihs, const hit& rhs) {
    return ((ihs.ctg < rhs.ctg) || (ihs.ctg == rhs.ctg && ihs.pos < rhs.pos));
}

//funtion to sort a pair of kmers  by score
bool compare_score(const hit& ihs, const hit& rhs) {
    return ((ihs.score > rhs.score));
}



//class to handle illumina paired reads
class Illumina{
private:
    //variables related to library
    int inferred_insert_size;
    int inferred_orientation;
    int inferred_std;
    //advanced options for illumina
    int vhs_ill;
    int ms_ill;
    int msr_ill;

    vector<bool> obs_ori;
    string fwdfile;
    string revfile;
    string prefix;
    int compresed;
    //return the reverse complement of a given orientation
    uint get_rev_ori(uint ori){
        switch(ori)
        {
            case 0: //FF -> RR
                return 3;
            case 1://RF -> FR
                return 2;
            case 2://FR->RF
                return 1;
            case 3://RR->FF
                return 0;
            default:
                cout << "error in orientation "<<ori<<" unknow value"<<endl;
                exit(1);
        }
    }
    //return the orientation given two contigs
    uint get_orientation(int  a, int b){
        switch(a){
            case 0:
                switch(b){
                    case 0:
                        return 0; //FF
                    case 1:
                        return 2; //FR
                    default:
                        cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                        exit(1);
                }
            case 1:
                switch(b){
                    case 0://RF
                        return 1;
                    case 1:
                        return 3;//RR
                    default:
                        cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                        exit(1);
                }
            default:
                cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                exit(1);
        }
    }
    vector<bool> get_ori_fwdrev(void){
        vector<bool> ori(2,0);
        switch(inferred_orientation){
            case 0:
                return ori;
            case 1:
                ori[0]=1;
                return ori;
            case 2:
                ori[1]=1;
                return ori;
            case 3:
                ori[0]=1;
                ori[1]=1;
                return ori;
            default:
                cout << "error in orientation "<<inferred_orientation<<" unknow value"<<endl;
                exit(1);
        }
    }
    //funtion to estimate the orienation and insert size given from a number of pairs
    void get_pairs_dist_orientation(uint number_pairs){
        //we open the fastq files
        kseq fwd;
        gzFile fp1 = gzopen(fwdfile.c_str(),"r");
        FunctorZlib r1;
        kstream<gzFile, FunctorZlib> ks1(fp1, r1);

        kseq rev;
        gzFile fp2 = gzopen(revfile.c_str(),"r");
        FunctorZlib r2;
        kstream<gzFile, FunctorZlib> ks2(fp2, r2);

        int l=0,c=0;
        uint64_t spairs=0;//mapped kmers pairs
        uint64_t tpairs=0;//mapped total pairs
        uint64_t asking=0;
        uint64_t dist=0;//distance
        map< uint, uint64_t > ori_insert; //to infer orientation
        map< uint, uint64_t > ori_count; //
        map<uint, vector<int> > values;
        cout << "Estimating insert-size and orientation using "<<number_pairs<<" number of read pairs"<<endl;
        double t_begin,t_end; struct timeval timet;
        gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
        uint kmer_size;
        kmer_size = kdb->getKmerSize();
        hit indexs,indexf;

        while(((l = ks1.read(fwd)) >= 0) && ((c = ks2.read(rev)) >= 0) && (tpairs < number_pairs)) {
            int max_i= fwd.seq.length() - kmer_size;
            int max_j= rev.seq.length() - kmer_size;
            //we skypt pairs of reads shorter than the kmer length
            if(max_i < 0 || max_j < 0){
                continue;
            }
            tpairs++;

            uint64_t hVal1=0,hfwd1=0,hrev1=0;
            hVal1 = NTPC64(fwd.seq.substr(0, kmer_size).c_str(), kmer_size, hfwd1, hrev1); // initial hash value for fwd
            vector<hit> tfwd;
            int max_h=0;
            for (int i = 0; i <= max_i; i++) {
                indexf = kdb->lookup_kmer(hfwd1, hrev1);
                if (indexf.map) {
                    tfwd.push_back(indexf);
                    max_h++;
                    if(max_h > 9){
                        break;
                    }
                }
                hVal1 = NTPC64(fwd.seq[i], fwd.seq[i + kmer_size], kmer_size, hfwd1, hrev1); //recursive hash k+L
                asking++;
            }
            //we try to map the reverse if we mapped the forward read
            if(tfwd.size() > 9){
                uint64_t hVal2=0,hfwd2=0,hrev2=0;
                hVal2 = NTPC64(rev.seq.substr(0, kmer_size).c_str(), kmer_size, hfwd2, hrev2); // initial hash value for rev
                vector<hit> trev;
                int max_h=0;
                for (int i = 0; i <= max_j; i++) {
                    indexf = kdb->lookup_kmer(hfwd2, hrev2);
                    if (indexf.map) {
                        trev.push_back(indexf);
                        max_h++;
                        if (max_h > 9) {
                            break;
                        }
                    }
                    hVal2 = NTPC64(rev.seq[i], rev.seq[i + kmer_size], kmer_size, hfwd2, hrev2); //recursive hash k+L
                    asking++;
                }
                //we got the maximum score using the function
                if(trev.size() > 9){
                    //we sort the hits by contig and post
                    sort(tfwd.begin(),tfwd.end(),compare_by_ctg_pos);
                    sort(trev.begin(),trev.end(),compare_by_ctg_pos);
                    //now we get the mayority or in the case of no mayority the max
                    hit bfwd=compute_maj_max(tfwd,fwd.seq.length());
                    hit brev=compute_maj_max(trev,rev.seq.length());
                    //we fill the insert size and the orientation if we had a good pairs of reads
                    if(bfwd.score > 7  && brev.score > 7){
                        if(bfwd.ctg == brev.ctg){
                            //we check if they are in foward or reverse
                            uint64_t insert_size=0;
                            uint ori=get_orientation((int)bfwd.strand,(int)brev.strand);
                            if(bfwd.pos > brev.pos){
                                insert_size=(bfwd.pos-brev.pos);
                                ori=get_rev_ori(ori);
                            }else{
                                insert_size=(brev.pos-bfwd.pos);
                            }
                            //we define a maximum  insert size
                            ori_count[ori]++;
                            ori_insert[ori] += insert_size;
                            dist += insert_size; //both contigs are in reverse
                            values[ori].push_back(int(insert_size));
                            spairs++;
                        }

                    }

                }

            }

        }
        gzclose(fp1);//closing file descriptors
        gzclose(fp2);//closing file descriptors

        gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
        double elapsed3 = t_end - t_begin;
        cout << "End of insert-size and orientation inference\n";
        printf("Query of %llu kmer-pairs  in %.2fs\n", asking,elapsed3);
        printf("Query of %llu pairs  in %.2fs\n", tpairs,elapsed3);
        printf("found %llu kmer-pairs at average dist=%f\n",spairs,float(dist/spairs));
        //we infer the orientaton of the librery
        uint64_t max_count=0;
        for (map< uint , uint64_t>::iterator it=ori_count.begin(); it!=ori_count.end(); ++it){
            cout << "Orientation "<<it->first << " " << it->second << " average_insert_size= "<< int(ori_insert[it->first]/it->second) << endl;
            if(max_count < it->second){
                max_count=it->second;
                inferred_orientation=it->first;
            }
        }

        //to compute average and insert size we will use [Q2,Q3] of sorted array from observations
        sort(values[inferred_orientation].begin(),values[inferred_orientation].end());
        //we discard 10% outlayers from begin and end
        int lower=int(0.1 * values[inferred_orientation].size());
        int upper=int(0.9 * values[inferred_orientation].size());
        uint64_t q2q3_avg=0;
        for (int i = lower; i < upper ; i++) {
                     q2q3_avg+=values[inferred_orientation][i];
        }
        q2q3_avg=uint64_t(q2q3_avg/(upper - lower));
        //cout << "inferred average from [Q2,Q3] " << q2q3_avg<<" lower="<<lower<<" upper="<<upper<<" total obs="<<upper-lower<<endl;
        inferred_insert_size=int(q2q3_avg);
        //compute std
        double sqtotal=0;
        for(int i=lower;i<upper; i++){
            // cout << values[inferred_orientation][i]<<" "<<inferred_insert_size<<" "<<inferred_insert_size-values[inferred_orientation][i]<<endl;
            sqtotal+=pow((inferred_insert_size-values[inferred_orientation][i]),2);
            //cout << sqtotal <<endl;
        }
        //cout << "inferred std from [Q2,Q3] "<< sqrt(sqtotal2/(upper-lower)) <<endl;
        //we split the orientation in fwd and rev
        obs_ori=get_ori_fwdrev();
        inferred_std=int(sqrt(sqtotal/(upper-lower)));
        //the maximal std that we allow is 15% of the inferred insert_size
        if(((float)inferred_std/(float)inferred_insert_size) > 0.15){
            cout << "adjusting inferred_std to 10% of inferred insert_size due that std > 0.15 "<< inferred_std <<endl;
            //we adjust the inferred_std to 10% of the average insert_size distribution
            inferred_std=int((float)inferred_insert_size * 0.1);
            cout << "new adjusted std = "<< inferred_std <<endl;
        }
        cout << "Inferred Orientation = " << inferred_orientation <<endl;
        cout << "Inferred average insert size = " << inferred_insert_size<<endl;
        cout <<"Inferred std of insert size = " << inferred_std<<endl;
        cout << "Number of pairs to infer avg/std = "<<upper-lower<<endl;
    }

public:
    pthread_mutex_t mutex;//to access to read files
    pthread_mutex_t pinfo;//to print information from the threads
    //variables related to KmerDB
    KmerDB *kdb;//pointer to KmerDB for get mapping positions
    //constructor of the class
    Illumina(string fprefix, string file1, string file2, int vhs_ill1, int ms_ill1, int msr_ill1, KmerDB *kmer_db){
        fwdfile=file1;
        revfile=file2;
        if((fwdfile.find(".gz")!=string::npos && revfile.find(".gz") !=string::npos) || (fwdfile.find(".gzip")!=string::npos && revfile.find(".gzip")!=string::npos)){
            cout << fwdfile <<" "<<revfile<<" has file extension of gzip files(*.gz or *.gzip)\n";
            compresed=1;//not compressed
        }else{
            if((fwdfile.find(".fq")!=string::npos && revfile.find(".fq") !=string::npos) || (fwdfile.find(".fastq")!=string::npos && revfile.find(".fastq")!=string::npos)){
                cout << fwdfile <<" "<<revfile<<" are text based *.fq or *.fastq\n";
                compresed=2;//gzip compressed
            }else{
                cout << fwdfile <<" "<<revfile<<" has not the file extension .gz or .fq or *.fastq\n";
                cout << "not supported fastq file\n";
                exit(1);
            }
        }
        prefix=fprefix;
        //we obtain a connector to KmerDB
        kdb=kmer_db;
        cout << "Total Kmers stored in datase "<<kdb->get_number_kmers()<<endl;
        //set the advanced options
        vhs_ill=vhs_ill1;
        ms_ill=ms_ill1;
        msr_ill=msr_ill1;
        //we compute the insert size distribution and orientation using 10000 pairs
        get_pairs_dist_orientation(100000);
    }
    void mapshortreads(uint samout, uint ntreads){
        gzFile fp1 = gzopen(fwdfile.c_str(),"r");
        FunctorZlib r1;
        kstream<gzFile, FunctorZlib> ks1(fp1, r1);

        gzFile fp2 = gzopen(revfile.c_str(),"r");
        FunctorZlib r2;
        kstream<gzFile, FunctorZlib> ks2(fp2, r2);
        //we close the fastq file because we already split the files into other structures
        //we create the threads for concurrency in the alignments step
        pthread_t *tab_threads= new pthread_t [ntreads];
        //we create the mutex to reads sequence from fasta file in an atomic way
        pthread_mutex_init(&mutex, NULL);
        pthread_mutex_init(&pinfo, NULL);
        ofstream samfile;
        uint kmer_size=kdb->getKmerSize();

        samfile.open(prefix+".FastSG_K"+to_string(kmer_size)+".sam");
        //we need to output sam header to the samfile
        samfile <<"@HD\tVN:1.0\tSO:unsorted"<<endl;
        for (int i=0; i < kdb->get_number_seq(); i++){
                string seqn=kdb->get_seq_name(i);
                uint seql=kdb->get_seq_len(i);
                samfile <<"@SQ\tSN:" << seqn <<"\tLN:"<< seql <<endl;
        }
        samfile <<"@PG\tID:" << prefix << "\tPN:FAST-SG\tVN:1.0\tCL:FAST-SG pair-end modes"<<endl;

        targ tmp;//struct to share data whit threads
        tmp.mut=mutex;//MUTEX to control access to files
        tmp.pin=pinfo;//MUTEX to control the ouput fo files
        tmp.db=kdb;//pointer to MPHF

        tmp.ks1=&ks1;//pointer to fastq files fwd
        tmp.ks2=&ks2;//pointer to fastq files reverse
        tmp.ill=this;
        //we need to create samoutput
        uint sampc=1;

        tmp.sam=&samfile;
        tmp.samflag=&samout;
        tmp.samp=&sampc;//we init the counter in 0

        //create the threads
        for(int ii=0;ii<ntreads;ii++) {
            pthread_create(&tab_threads[ii], NULL, threaded_shortmap2, &tmp);
        }
        //wait for thread to finish
        for(int ii=0;ii<ntreads;ii++)
        {
            pthread_join(tab_threads[ii], NULL);
        }
        //closing samfile
        samfile.close();
        gzclose(fp1);
        gzclose(fp2);
        //print some mapping information
        cout << "Number of total mapped pairs: " << sampc << endl;

    }
    //function to compute the majority or the max of a sorted hit vector
    hit compute_maj_max(vector<hit> &tfwd, int seqlen){
        //find the majority strict
        hit mfwd;
        int mc=0;
        hit maxfwd;
        int max=0;
        for (int i = 0; i < tfwd.size() ; i++) {
            if(mc==0){
                mfwd=tfwd[i];
                mc=1;
            } else{
                //means that both kmers comes from the same region
                if((mfwd.ctg == tfwd[i].ctg) && (tfwd[i].pos-mfwd.pos <= seqlen)){
                    mc++;
                }else{
                    mc--;
                }
            }
            //we try to keep a reference to max value obtained
            if(mc >= max){
                max=mc;
                maxfwd=mfwd;
            }
        }
        //we set the score of the best single end hit
        maxfwd.score=max;
        //we return the strict mayority or the max element
        return maxfwd;
    }
    //function to compute the score given the vector of mers
    vector<hit> compute_score(vector<hit> &tfwd, int seqlen){
        //we sort the array by ctg star
        sort(tfwd.begin(),tfwd.end(),compare_by_ctg_pos);
        hit mfwd=tfwd[0];
        int mc=0;
        vector<hit> index;
        int add_last=0;
        for (int i = 0; i < tfwd.size() ; i++) {
            //means that both kmers comes from the same region
            if((mfwd.ctg == tfwd[i].ctg) && (tfwd[i].pos-mfwd.pos <= seqlen)){
                    mc++;
            }else{  //we add the last element to the array
                    mfwd.score=mc;
                    index.push_back(mfwd);
                    //we add a new mfwd
                    mc=1;
                    mfwd=tfwd[i];
                    mfwd.score=1;
                    if(i == tfwd.size() - 1){
                            add_last=1;
                    }
            }
        }
        //we add the last kmer to the set because is has an score of 1
        if(add_last){
            index.push_back(mfwd);
        }
        //we add the kmer to the set because it has score max
        if(mc == tfwd.size()){
            mfwd.score=mc;
            index.push_back(mfwd);
        }
        //we sort by score in order to print the output
        sort(index.begin(),index.end(),compare_score);
        return index;
    }

    //funtion to compute the spam distance between two pairs of kmers assuming the rigth orientation
    int compute_distance_fragment2(hit &fwd, hit &rev){

        int64_t distance_fwd = 0;
        int64_t distance_rev = 0;
        int64_t spam=0;
        //for fragment within the same contig
        int p_ori=get_orientation(int(fwd.strand),int(rev.strand));
        int rp_ori=get_rev_ori(p_ori);//compute the reverse order

        //easy case two reads are in the same contig
        if (fwd.ctg == rev.ctg) {
                //we need to check that the orientation is compatible with the current observed orientation
            if(p_ori == inferred_orientation || rp_ori == inferred_orientation) {
                return (abs(int(fwd.pos - rev.pos)));
               }
        }

        //the reads connect two contig we check if they are in the correct orientation
        if(p_ori == inferred_orientation || rp_ori == inferred_orientation) {
                switch (inferred_orientation) {
                    case 0: //-> -> FWD
                        if( fwd.strand == 0 ){
                            distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                            if( rev.strand == 0 ){
                                // read ori is + +, contig ori should be + +
                                distance_rev = rev.pos;
                            }
                            else{
                                // read ori is + -; contig ori should be + -
                                distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                            }
                        }
                        else{
                            distance_fwd = fwd.pos;
                            if( rev.strand == 0 ){
                                // read ori is - +; contig ori should be - +
                                distance_rev = rev.pos;
                            }
                            else{
                                // read ori is - -; contig ori should be + +
                                distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                            }
                        }
                        spam = distance_fwd + distance_rev;
                        return spam;
                    case 2: //-> <- IN
                            //positive strand
                        if( fwd.strand == 0 ){
                            distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;//+
                            distance_rev = rev.pos;//-
                        }
                        else{
                            distance_fwd = fwd.pos; //-
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos; //+
                        }
                        spam = distance_fwd + distance_rev;
                        return spam;
                    case 1://<- -> OUT
                        // orientation is <-...-> (-,+)
                        if( fwd.strand == 0 ){
                            distance_fwd = fwd.pos;
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                        }
                        else{
                            distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                            distance_rev = rev.pos;
                        }
                        spam = distance_fwd + distance_rev;
                        return spam;
                    case 3: //<- <- REV
                        if( fwd.strand == 0 ){
                            distance_fwd = fwd.pos;
                            if( rev.strand == 0 ){
                                // read ori is + +, contig ori should be - -
                                distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                            }
                            else{
                                // read ori is + -; contig ori should be - +
                                distance_rev = rev.pos;
                            }
                        }
                        else{
                            distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                            if( rev.strand == 0 ){
                                // read ori is - +; contig ori should be + -
                                distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                            }
                            else{
                                // read ori is - -; contig ori should be + +
                                distance_rev = rev.pos;
                            }
                        }
                        spam = distance_fwd + distance_rev;
                        return spam;
                    default:
                        cout << "error in orientation " << inferred_orientation << " unknow value" << endl;
                        exit(1);
                }

        }

        return int(spam);
    }

    //function to rescue a pair with a low score using its mate we return the hit
    hit rescue_mate(hit &bhit, vector<hit> &mates,bool fwdrev ){
        int distance=0; //assuming foward first
        hit tmp;
        int min_spam=inferred_insert_size - 4*inferred_std;
        if(min_spam < 0){
            min_spam=0;
        }
        int max_spam=inferred_insert_size + 4*inferred_std;

        //we rescue only  if one pair satisfice the constrains of orientation and distance
        vector<int> rpairs;

        for (int i = 0; i <mates.size() ; i++) {
            distance=0;
                //best hit is in foward read
                if (fwdrev) {
                    //to try to rescue we need at least 4 unique mers supporting the position
                    if(mates[i].score > msr_ill) {
                        distance = compute_distance_fragment2(bhit, mates[i]);
                    }
                } else {//best hit is in reverse read
                    if(mates[i].score > msr_ill) {
                        distance = compute_distance_fragment2(mates[i], bhit);
                    }
                }

            if (distance <= max_spam && distance >= min_spam) {
                rpairs.push_back(i);
            }
        }
        //rescue was sucess because we found only one position that satisfices orientation and distance
        if(rpairs.size() == 1){
            //we found only one posible pair to perform the rescue
            tmp=mates[rpairs[0]];
            //if within a contig we give an score of 12 = 60
            if(bhit.ctg == tmp.ctg) {
                tmp.score=12;
            }//else{
               //we give the score of the rescue in this case
                //tmp.score=11;
            //}
            return tmp;
        }
        return tmp;
    }

    //small functions to return the insert size distribution
    int get_inferred_insert_size(void){
        return inferred_insert_size;
    }
    //infer std
    int get_inferred_std(void){
        return inferred_std;
    }
    //small funtion to return the clone name of a read
    string get_clone_name(string s){
        string s2=s;
        if (s.length() > 2 && s[s.length()-2] == '/' )
            s2=s.substr(0,s.length()-2);
        return s2;
    }
    //funtions for get advanced values
    int get_vhs(void){
        return vhs_ill;
    }
    int get_ms(void){
        return ms_ill;
    }
    int get_msr(void){
        return msr_ill;
    }
};



//first try to found one kmer in the foward reads if success try to find another kmer in the second read, if the two kmers are found we print the hit to the samfile
void * threaded_shortmap2 (void* args)
{
    //casting of variables
    targ *tmp = (targ*) args;
    pthread_mutex_t *mutex = &tmp->mut;
    pthread_mutex_t *pinfo = &tmp->pin;
    KmerDB *dbt=(KmerDB *)tmp->db;
    Illumina *ill=(Illumina *)tmp->ill;

    //compressed
    kstream<gzFile, FunctorZlib> *ks1= (kstream<gzFile, FunctorZlib> *) tmp->ks1;
    kstream<gzFile, FunctorZlib> *ks2= (kstream<gzFile, FunctorZlib> *) tmp->ks2;


    uint *samout=(uint *) tmp->samflag;
    ofstream *sam;
    uint *samp;
    uint *single;

    sam=(ofstream *) tmp->sam;
    samp=(uint *) tmp->samp;
    //init of locals variables
    uint number_seq=0;//variable to control the number of sequence that we read before the mapping
    uint kmer_size;
    kmer_size = (uint)dbt->getKmerSize();
    string qual(kmer_size,'I');
    string file_to_process="";
    hit indexs,indexf;
    vector<hit> fwdh;//local storage of hits
    vector<hit> revh;//local storage of hits
    //we start the mapping
    uint64_t tpairs=0;//total pairs mapped in thread
    uint64_t tmap=0;//total asked kmers
    uint64_t pmap=0;//total mapped-kmers pairs
    uint64_t tmapp=0;//total mapped-pairs
    //hashing related variables
    uint64_t hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
    //we get both kmers in foward and reverse mers
    uint64_t hVal1 = 0, hVal2 = 0;
    //determine min/max allowed spam for a pair of reads
    int min_spam=ill->get_inferred_insert_size() - ill->get_inferred_std() * 4;
    int max_spam=ill->get_inferred_insert_size() + ill->get_inferred_std() * 4;
    if(min_spam < 0){
        min_spam=0;
    }
    //advanded options
    int vhs=ill->get_vhs();
    int ms=ill->get_ms();
    int msr=ill->get_msr(); //for future use

    //cout << "printing advanced options from thread\n" << vhs << " " << ms <<" "<<msr<<endl;

    int l=0,c=0;
    kseq fwd;
    kseq rev;
    vector<string> localseq;
    bool dojob=1;//variable to control thread execution
    string kmersizeS=to_string(kmer_size);
    //we check the time spent
    double t_begin,t_end; struct timeval timet;
    gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
    while(dojob) {
        //we try to lock the mutex to reads sequences from files
        number_seq=0;
        //we got acess to the global fasta file
        pthread_mutex_lock(mutex);
        //we read 1mb reads from  fwd and rev fastq/fasta file keep reads names and assing an unique id
        while ( (number_seq < 500000) && ((l = ks1->read(fwd)) >= 0) && ((c = ks2->read(rev)) >= 0) ){
            localseq.push_back(to_string(number_seq)+" "+fwd.seq+" "+rev.seq+" "+fwd.name+" "+rev.name);
            number_seq++;
        }
        pthread_mutex_unlock(mutex);
        //we ask if the need to map reads
        if(localseq.size() > 0){
            dojob=1;
        }else{
            dojob=0;
        }
        //we do the mapping
        if(dojob){
            //max/min spam values

            for(int l=0;l<localseq.size();l++){
                string line=localseq[l];
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                //we check that the split was ok
                assert(tokens.size() == 5);

                int max_i = tokens[1].length() - kmer_size;
                int max_j = tokens[2].length() - kmer_size;
                tpairs++;
                //we skyp pairs of reads shorter than the kmer length
                if(max_i < 0 || max_j < 0){
                    continue;
                }
                //init hash values before mapping of kmers
                hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
                //we get both kmers in foward and reverse mers
                hVal1 = 0, hVal2 = 0;
                //init foward read
                hVal1 = NTPC64(tokens[1].substr(0, kmer_size).c_str(), kmer_size, hfwd1, hrev1); // initial hash value for fwd
                vector<hit> tfwd;
                int max_h=0;
                for (int i = 0; i <= max_i; i++) {
                    string fkmer = tokens[1].substr(i, kmer_size);
                    indexf = dbt->lookup_kmer(hfwd1, hrev1, fkmer);
                    tmap++;
                    if (indexf.map) {
                        tfwd.push_back(indexf);
                        max_h++;
                        if(max_h > vhs){
                            break;
                        }

                    }
                    hVal1 = NTPC64(tokens[1][i], tokens[1][i + kmer_size], kmer_size, hfwd1, hrev1); //recursive hash k+L
                }
                //if we found one kmer in foward we try the reverse, we can improve asking if we find at least 4 kmers
                if(tfwd.size() > 0){
                    //init reverse read
                    hVal2 = NTPC64(tokens[2].substr(0, kmer_size).c_str(), kmer_size, hfwd2, hrev2); //initial hash value for re
                    vector<hit> trev;
                    int max_h=0;
                    for (int i = 0; i <= max_j; i++) {
                        string skmer = tokens[2].substr(i, kmer_size);
                        indexs = dbt->lookup_kmer(hfwd2, hrev2, skmer);
                        tmap++;
                        if (indexs.map) {
                            trev.push_back(indexs);
                            pmap++;
                            max_h++;
                            if(max_h > vhs){
                                break;
                            }
                        }

                        hVal2 = NTPC64(tokens[2][i], tokens[2][i + kmer_size], kmer_size, hfwd2, hrev2); //recursive hash k+L
                    }
                    //we can improve just asking if we find at least 4 kmers
                    if(trev.size() > 0) {
                        vector<hit> tfwd2=ill->compute_score(tfwd,tokens[1].length());
                        vector<hit> trev2=ill->compute_score(trev,tokens[2].length());
                        hit maxfwd=tfwd2[0];
                        hit maxrev=trev2[0];
                        maxfwd.seqname=ill->get_clone_name(tokens[3]);
                        maxrev.seqname=ill->get_clone_name(tokens[4]);

                        hit rescue;
                        //if two reads are majority we don't touch them and we trush in their positions
                        if(maxfwd.score > ms && maxrev.score > ms){
                            fwdh.push_back(maxfwd);
                            revh.push_back(maxrev);
                            tmapp++;
                        }else{
                            //we rescue if at least one of each other is majority
                            if(maxfwd.score > ms){
                                rescue = ill->rescue_mate(maxfwd, trev2, 1);
                                rescue.seqname = maxrev.seqname;
                                if (rescue.score > 0) {
                                    fwdh.push_back(maxfwd);
                                    revh.push_back(rescue);
                                    tmapp++;
                                }
                            }else{
                                  if(maxrev.score > ms) {
                                      rescue = ill->rescue_mate(maxrev, tfwd2, 0);
                                      rescue.seqname = maxfwd.seqname;
                                      if (rescue.score > 0) {
                                          fwdh.push_back(rescue);
                                          revh.push_back(maxrev);
                                          tmapp++;
                                      }
                                  }
                            }
                        }
                        //clear tmp memory for container
                        tfwd2.erase(tfwd2.begin(),tfwd2.end());
                        trev2.erase(trev2.begin(),trev2.end());
                    }
                    //clear tmp memory for container
                    tfwd.erase(tfwd.begin(),tfwd.end());
                    trev.erase(trev.begin(),trev.end());
                }
            }
            //we clear the localseq variable
            localseq.erase(localseq.begin(),localseq.end());
            //we print the result in sam  format using singlest or paired reads, it depends of the scaffolder
            pthread_mutex_lock(pinfo);
            //means write resuls as single end reads
            if(*samout == 1){
                //we need to write the results to the bam file if is needed
                    //we print reads as single end reads for BOSS/SCAFFMATCH
                    for (int j = 0; j < fwdh.size(); ++j) {
                        string fsam_flag = "0";
                        string rsam_flag = "0";
                        //if the kmer hit the reverse strand I need to ouput the reverse kmer in the samfile
                        if (fwdh[j].strand) {
                            fsam_flag = "16";
                        }
                        if (revh[j].strand) {
                            rsam_flag = "16";
                        }
                        string tmp1 =
                                fwdh[j].seqname + "\t" + fsam_flag + "\t" + dbt->get_seq_name(fwdh[j].ctg) + "\t" +
                                to_string(fwdh[j].pos) + "\t" + to_string(fwdh[j].score * 10) + "\t" +
                                        kmersizeS + "M\t*\t0\t0\t" + fwdh[j].kseqm + "\t" + qual +
                                        "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        string tmp2 =
                                revh[j].seqname + "\t" + rsam_flag + "\t" + dbt->get_seq_name(revh[j].ctg) + "\t" +
                                to_string(revh[j].pos) + "\t" + to_string(revh[j].score * 10) + "\t" +
                                        kmersizeS + "M\t*\t0\t0\t" + revh[j].kseqm + "\t" + qual +
                                        "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";

                        (*samp)++;//we increment the pointer of number of mapped kmers pairs
                        sam->write(tmp1.c_str(), tmp1.length());
                        sam->write(tmp2.c_str(), tmp2.length());
                    }
                }else{
                    //we print reads as single end reads pair-ends for BESST/OPERA scaffolders
                    for (int j = 0; j < fwdh.size(); ++j) {
                        int fsam_flag = 0;
                        int rsam_flag = 0;
                        //if the kmer hit the reverse strand I need to ouput the reverse kmer in the samfile
                        if (fwdh[j].strand) {
                            fsam_flag += 16;
                            rsam_flag+=32;
                        }
                        if (revh[j].strand) {
                            rsam_flag = 16;
                            fsam_flag+=32;
                        }

                        fsam_flag+=65;//paired and read first in pair
                        rsam_flag+=129;//paired and read second in pair

                        string tmp1="";
                        string tmp2="";
                        //if within contigs we can compute the distance
                        if(fwdh[j].ctg == revh[j].ctg ) {
                            int distance1=(int)(revh[j].pos - fwdh[j].pos);
                             tmp1 = fwdh[j].seqname + "\t" + to_string(fsam_flag) + "\t" + dbt->get_seq_name(fwdh[j].ctg) + "\t" +
                                    to_string(fwdh[j].pos) + "\t" + to_string(fwdh[j].score * 10) + "\t" +
                                     kmersizeS + "M\t=\t"+to_string(revh[j].pos)+"\t"+to_string(distance1)+"\t" + fwdh[j].kseqm + "\t" + qual +
                                     "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                            int distance2=(int)(fwdh[j].pos - revh[j].pos);
                             tmp2 = revh[j].seqname + "\t" + to_string(rsam_flag) + "\t" + dbt->get_seq_name(revh[j].ctg) + "\t" +
                                    to_string(revh[j].pos) + "\t" + to_string(revh[j].score * 10) + "\t" +
                                     kmersizeS + "M\t=\t"+to_string(fwdh[j].pos)+"\t"+to_string(distance2)+"\t" + revh[j].kseqm + "\t" + qual +
                                     "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        }else{

                            tmp1 = fwdh[j].seqname + "\t" + to_string(fsam_flag) + "\t" + dbt->get_seq_name(fwdh[j].ctg) + "\t" +
                                   to_string(fwdh[j].pos) + "\t" + to_string(fwdh[j].score * 10) + "\t" +
                                    kmersizeS + "M\t"+dbt->get_seq_name(revh[j].ctg)+"\t"+to_string(revh[j].pos)+"\t0\t" + fwdh[j].kseqm + "\t" + qual +
                                    "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";

                            tmp2 = revh[j].seqname + "\t" + to_string(rsam_flag) + "\t" + dbt->get_seq_name(revh[j].ctg) + "\t" +
                                   to_string(revh[j].pos) + "\t" + to_string(revh[j].score * 10) + "\t" +
                                    kmersizeS + "M\t"+dbt->get_seq_name(fwdh[j].ctg)+"\t"+to_string(fwdh[j].pos)+"\t0\t" + revh[j].kseqm + "\t" + qual +
                                    "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        }

                        (*samp)++;//we increment the pointer of number of mapped kmers pairs
                        sam->write(tmp1.c_str(), tmp1.length());
                        sam->write(tmp2.c_str(), tmp2.length());
                    }
                }

            pthread_mutex_unlock(pinfo);
            //we need to clear the local variables of the thread before the next execution
            //clear mapped hits
            fwdh.erase(fwdh.begin(),fwdh.end());
            revh.erase(revh.begin(),revh.end());
        }

    }
    gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
    double elapsed3 = t_end - t_begin;
    //we print some stats from thread information
    pthread_mutex_lock(pinfo);
    cout << "End of mapping in thread\n";
    printf("Query of %llu kmer-pairs  in %.2fs rate %.2f per second\n", tmap,elapsed3,tmap/elapsed3);
    printf("Query of %llu pairs  in %.2fs\n", tpairs,elapsed3);
    printf("Number of total mapped kmer-pairs %llu\n", pmap);
    printf("Number of total mapped reads-pairs %llu\n", tmapp);
    pthread_mutex_unlock(pinfo);
    return NULL;
}

//class to handle long reads data
class LongReads{
private:
    string seqfile;
    string prefix;
    vector<int> inserts;//vector that store the desired insert sizes
    vector<int> inferred_insert_sizes;
    int compresed;
    //advanced long read variables
    int vhs_long;
    int ms_long;
    int msr_long;
    int ssr_long;

    uint get_rev_ori(uint ori){
        switch(ori)
        {
            case 0: //FF -> RR
                return 3;
            case 1://RF -> FR
                return 2;
            case 2://FR->RF
                return 1;
            case 3://RR->FF
                return 0;
            default:
                cout << "error in orientation "<<ori<<" unknow value"<<endl;
                exit(1);
        }
    }
    //return the orientation given two contigs
    uint get_orientation(int  a, int b){
        switch(a){
            case 0:
                switch(b){
                    case 0:
                        return 0; //FF
                    case 1:
                        return 2; //FR
                    default:
                        cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                        exit(1);
                }
            case 1:
                switch(b){
                    case 0://RF
                        return 1;
                    case 1:
                        return 3;//RR
                    default:
                        cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                        exit(1);
                }
            default:
                cout << "error in orientation "<<a<<" "<<b<<" unknow value"<<endl;
                exit(1);
        }
    }



    void get_pairs_dist_orientation(uint number_pairs) {
        //we open the fastq files
        kseq fwd;
        gzFile fp1 = gzopen(seqfile.c_str(), "r");
        FunctorZlib r1;
        kstream<gzFile, FunctorZlib> ks1(fp1, r1);
        int l = 0, c = 0;
        uint64_t spairs = 0;//mapped kmers pairs
        uint64_t dist = 0;//distance
        //vector to store the observed insert sizes
        vector<int> values;
        cout << "Estimating insert-sizes using " << number_pairs << " number of long reads" << endl;
        double t_begin, t_end;
        struct timeval timet;
        gettimeofday(&timet, NULL);
        t_begin = timet.tv_sec + (timet.tv_usec / 1000000.0);

        //init of locals variables
        uint number_seq=0;//variable to control the number of sequence that we read before the mapping
        uint kmer_size;
        kmer_size = (uint)kdb->getKmerSize();
        hit indexs,indexf;
        vector<hit> fwdh;//local storage of hits
        vector<hit> revh;//local storage of hits
        //we start the mapping
        uint64_t tpairs=0;//total pairs mapped in thread
        uint64_t tmap=0;//total asked kmers
        uint64_t pmap=0;//total mapped-kmers pairs
        uint64_t tmapp=0;//total mapped-pairs
        uint64_t totalseqs=0;
        //hashing related variables
        uint64_t hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
        //we get both kmers in foward and reverse mers
        uint64_t hVal1 = 0, hVal2 = 0;

        vector<string> localseq;//storage of local seq

        //we read number_pairs reads
        while ((number_seq < number_pairs) && ((l = ks1.read(fwd)) >= 0) ) {
		//we continue reading if seq is shorter kmer_size
		if(fwd.seq.length() < kmer_size){
                	continue;
            	}
            //we store the long read in memory
            localseq.push_back(to_string(number_seq) + " " + fwd.seq + " " + fwd.name );
            number_seq++;
            totalseqs++;
        }

        //we need to iterate the seq of sequence in each distance
        for (int i = 0; i <inserts.size() ; i++) {
            int d=inserts[i];//this is the target insert size
            //now we iter each sequence picking pairs of mers
            for (int j = 0; j < localseq.size(); ++j) {
                string line=localseq[j];
                istringstream iss(line);
                vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                //we check that the split was ok //0=id,1=seq,2=name
                assert(tokens.size() == 3);
                //we took two part of sequences at a given distance(d)
                //we use sliding windows of 50bp after try the first map
                int max_d = tokens[1].length() - d;
                //int max_j = tokens[1].length() - d - 150;//max
                //tpairs++;
                //we skyp a sequence if is shorter than the length of d
                if(max_d < 0 ){
                    continue;
                }
                //we need to split the long read in read pairs of length 150 whit a moving window of 100bp
                for (int k = 0; k < max_d ; k+=100) {
                    string fwds=tokens[1].substr(k,200); //0+150
                    string revs=tokens[1].substr(k+d-200,200);//d-150
                    //we attempt to map the pair of subreads
                    int max_i = fwds.length() - kmer_size;
                    int max_j = revs.length() - kmer_size;
                    tpairs++;
                    //we skyp pairs of reads shorter than the kmer length
                    if(max_i < 0 || max_j < 0){
                        continue;
                    }
                    //init hash values before mapping of kmers
                    hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
                    //we get both kmers in foward and reverse mers
                    hVal1 = 0, hVal2 = 0;
                    //init foward read
                    hVal1 = NTPC64(fwds.substr(0, kmer_size).c_str(), kmer_size, hfwd1, hrev1); // initial hash value for fwd
                    vector<hit> tfwd;
                    int max_h=0;
                    for (int f = 0; f <= max_i; f++) {
                        string fkmer = fwds.substr(f, kmer_size);
                        indexf = kdb->lookup_kmer(hfwd1, hrev1, fkmer);
                        tmap++;

                        if (indexf.map) {
                            tfwd.push_back(indexf);
                            max_h++;
                            if(max_h > 9){
                                break;
                            }

                        }
                        hVal1 = NTPC64(fwds[f], fwds[f + kmer_size], kmer_size, hfwd1, hrev1); //recursive hash k+L
                    }
                    //if we found at least 4 kmer in foward we try the reverse
                    if(tfwd.size() > 5) {
                        //we reversecomp the revs to produce a pair in -> <- orientation
                        revs=this->revcomp(revs);
                        //init reverse read
                        hVal2 = NTPC64(revs.substr(0, kmer_size).c_str(), kmer_size, hfwd2,
                                       hrev2); //initial hash value for re
                        vector <hit> trev;
                        int max_h = 0;
                        for (int r = 0; r <= max_j; r++) {
                            string skmer = revs.substr(r, kmer_size);
                            indexs = kdb->lookup_kmer(hfwd2, hrev2, skmer);
                            tmap++;

                            if (indexs.map) {
                                trev.push_back(indexs);
                                pmap++;
                                max_h++;
                                if (max_h > 9) {
                                    break;
                                }
                            }
                            hVal2 = NTPC64(revs[r], revs[r + kmer_size], kmer_size, hfwd2, hrev2); //recursive hash k+L
                        }
                        //now we need to decide if we trusths in the positions
                        //we can improve just asking if we find at least 4 kmers
                        if(trev.size() > 5) {
                            vector<hit> tfwd2=this->compute_score(tfwd,fwds.length());
                            vector<hit> trev2=this->compute_score(trev,revs.length());
                            hit maxfwd=tfwd2[0];
                            hit maxrev=trev2[0];
                            //maxfwd.seqname=ill->get_clone_name(tokens[3]);
                            //maxrev.seqname=ill->get_clone_name(tokens[4]);
                            maxfwd.seqname=tokens[2]+"_"+to_string(d)+"_"+to_string(k);
                            maxrev.seqname=tokens[2]+"_"+to_string(d)+"_"+to_string(k);

                            //we will check if they are at the expected distance and orientation assumming a 10% standart dev
                            if(maxfwd.score > 5 && maxrev.score > 5 ){
                                //fwdh.push_back(maxfwd);
                                //revh.push_back(maxrev);
                                //I need to compute the distance within a contig pair for the current insert size
                                //easy case two reads are in the same contig
                                if (maxfwd.ctg == maxrev.ctg) {
                                    //we need to check that the orientation is compatible with the current observed orientation
                                        values.push_back(abs(int(maxfwd.pos - maxrev.pos)));
                                }
                                tmapp++;
                               }
                            //clear tmp memory for container
                            tfwd2.erase(tfwd2.begin(),tfwd2.end());
                            trev2.erase(trev2.begin(),trev2.end());
                            }
                            //clear tmp memory for container
                           trev.erase(trev.begin(),trev.end());
                           tfwd.erase(tfwd.begin(),tfwd.end());
                        }

                    }

                }
                //here we need to determine the insert size for this variable
                //cout <<d <<" "<<values.size()<<endl;
            if(values.size() >=300) {
                //to compute average and insert size we will use [Q2,Q3] of sorted array from observations
                sort(values.begin(), values.end());
                //we discard 10% outlayers from begin and end
                int lower = int(0.1 * values.size());
                int upper = int(0.9 * values.size());
                uint64_t q2q3_avg = 0;
                for (int i = lower; i < upper; i++) {
                    q2q3_avg += values[i];
                }
                q2q3_avg = uint64_t(q2q3_avg / (upper - lower));
                cout << "insert-size  d="<<d<<" observed average insert-size " << q2q3_avg<<" lower="<<lower<<" upper="<<upper<<" total obs="<<upper-lower<<endl;
                //I need to update the value of the distances
                this->inferred_insert_sizes[i]=int(q2q3_avg);
            }else{
                cout <<"Few observations to infere insert size distribution for d="<< d <<" "<<values.size()<<endl;
            }
                //we clean the container to the insert-size
                values.erase(values.begin(),values.end());

            }

        gzclose(fp1);//closing file descriptors
        gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
        double elapsed3 = t_end - t_begin;
        cout << "End of insert-size inference long reads\n";
        printf("Total time spent in insert size inference %.2fs\n", elapsed3);


    }

public:
    KmerDB *kdb;//pointer to KmerDB for get mapping positions
    pthread_mutex_t lmutex;//to access to read files
    pthread_mutex_t lpinfo;//to print information from the threads
    pthread_mutex_t *lpinfom;//to coordinates sam writing

    //constructor of the class
    LongReads(string fprefix, string file1, int vhs_long1, int ms_long1, int msr_long1, int ssr_long1, vector<int> &inserts, KmerDB *kmer_db){
        seqfile=file1;
        this->inserts=inserts;
        //we create an equivalent array
        this->inferred_insert_sizes=inserts;

        if((seqfile.find(".gz")!=string::npos ) || (seqfile.find(".gzip")!=string::npos)){
            cout << seqfile <<" "<<" has file extension of gzip files(*.gz or *.gzip)\n";
            compresed=1;//not compressed
        }else{
            if((seqfile.find(".fq")!=string::npos ) || (seqfile.find(".fastq")!=string::npos)){
                cout << seqfile <<" "<<" are text based *.fq or *.fastq\n";
                compresed=2;//gzip compressed
            }else{
                cout << seqfile <<" "<<" has not the file extension .gz or .fq or *.fastq\n";
                cout << "not supported fastq file\n";
                exit(1);
            }
        }
        prefix=fprefix;
        //we obtain a connector to KmerDB
        kdb=kmer_db;

        //we set the advanced variables
        vhs_long=vhs_long1;
        ms_long=ms_long1;
        msr_long=msr_long1;
        ssr_long=ssr_long1;

        cout << "Total Kmers stored in datase "<<kdb->get_number_kmers()<<endl;
        //we compute the insert size distribution and orientation using 1000 long reads
        get_pairs_dist_orientation(1000);

        for (int i = 0; i <this->inferred_insert_sizes.size() ; ++i) {
                cout << "Observed average insert size for d="<<this->inserts[i]<<" was "<<this->inferred_insert_sizes[i]<<endl;
        }
    }


    void maplongreads(uint samout, uint ntreads){
        //we open a file handler for the long reads
        gzFile fp1 = gzopen(seqfile.c_str(),"r");
        FunctorZlib r1;
        kstream<gzFile, FunctorZlib> ks1(fp1, r1);

        //we close the fastq file because we already split the files into other structures
        //we create the threads for concurrency in the alignments step
        pthread_t *tab_threads= new pthread_t [ntreads];
        //we create the mutex to reads sequence from fasta file in an atomic way
        pthread_mutex_init(&lmutex, NULL);
        pthread_mutex_init(&lpinfo, NULL);
        //we init a set of lpinfo to coordinate the writing to the several output files
        lpinfom = new pthread_mutex_t [inserts.size()];
        for (int l = 0; l <inserts.size() ; ++l) {
            pthread_mutex_init(&lpinfom[l],NULL);
        }

        //pthread_mutex_init(&lpinfo, NULL);
        ofstream *samfiles=new ofstream[inserts.size()];
        uint kmer_size=kdb->getKmerSize();
        for (int j = 0; j <inserts.size() ; ++j) {
            samfiles[j].open(prefix +".I"+to_string(inserts[j])+ ".FastSG_K" + to_string(kmer_size) + ".sam");
            cout << prefix +".I"+to_string(inserts[j])+ ".FastSG_K" + to_string(kmer_size) + ".sam" <<endl;
            //we need to output sam header to the samfile
            samfiles[j] << "@HD\tVN:1.0\tSO:unsorted" << endl;
            for (int i = 0; i < kdb->get_number_seq(); i++) {
                string seqn = kdb->get_seq_name(i);
                uint seql = kdb->get_seq_len(i);
                samfiles[j] << "@SQ\tSN:" << seqn << "\tLN:" << seql << endl;
            }
            samfiles[j] << "@PG\tID:" << prefix+".I"+to_string(inserts[j]) << "\tPN:FAST-SG\tVN:1.0\tCL:FAST-SG pair-end modes" << endl;
        }


        targ tmp;//struct to share data whit threads
       //Mutex
        tmp.mut=lmutex;//MUTEX to control access to files
        tmp.pinm=lpinfom;//MUTEX to control the ouput fo files
        tmp.pin=lpinfo;//Mutex to control printing of thread information
        tmp.db=kdb;//pointer to MPHF

        tmp.ks1=&ks1;//pointer to fastq file with long reads
        tmp.lon=this;//pointer to LongRead object
        //we need to create samoutput
        uint sampc=1;
        tmp.sam=samfiles;//pointer to array of samfiles
        tmp.samflag=&samout;//
        tmp.samp=&sampc;//we init the counter in 0

        //create the threads
        for(int ii=0;ii<ntreads;ii++) {
            pthread_create(&tab_threads[ii], NULL, threaded_longmap, &tmp);
        }
        //wait for thread to finish
        for(int ii=0;ii<ntreads;ii++)
        {
            pthread_join(tab_threads[ii], NULL);
        }
        //we close the seq file
        gzclose(fp1);
        //closing samfile
        //print some mapping information
        cout << "Number of total mapped pairs: " << sampc << endl;
        //we close the samfiles
        for (int k = 0; k <inserts.size() ; k++) {
               samfiles[k].close();
        }
    }

    //target inserts sizes
    vector<int> get_inserts(void){
        return this->inserts;
    }
    //observed inserts sizes
    vector<int> get_obs_inserts(void){
        return this->inferred_insert_sizes;
    }
    string revcomp(string kmer){
        int c0, c1;
        int klen=kmer.length();
        for (int i = 0; i < klen>>1; ++i) { // reverse complement sequence
            c0 = comp_tab[(int)kmer[i]];
            c1 = comp_tab[(int)kmer[klen - 1 - i]];
            kmer[i] = c1;
            kmer[klen - 1 - i] = c0;
        }
        if (klen & 1) // complement the remaining base
            kmer[klen>>1] = comp_tab[(int)kmer[klen>>1]];
        return kmer;
    }

    int compute_distance_fragment2(hit &fwd, hit &rev){

        int64_t distance_fwd = 0;
        int64_t distance_rev = 0;
        int64_t spam=0;
        int inferred_orientation=2;//orientation is always FR for long reads pairs
        //for fragment within the same contig
        int p_ori=get_orientation(int(fwd.strand),int(rev.strand));
        int rp_ori=get_rev_ori(p_ori);//compute the reverse order

        //easy case two reads are in the same contig
        if (fwd.ctg == rev.ctg) {
            //we need to check that the orientation is compatible with the current observed orientation
            if(p_ori == inferred_orientation || rp_ori == inferred_orientation) {
                return (abs(int(fwd.pos - rev.pos)));
            }
        }

        //the reads connect two contig we check if they are in the correct orientation
        if(p_ori == inferred_orientation || rp_ori == inferred_orientation) {
            switch (inferred_orientation) {
                case 0: //-> -> FWD
                    if( fwd.strand == 0 ){
                        distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                        if( rev.strand == 0 ){
                            // read ori is + +, contig ori should be + +
                            distance_rev = rev.pos;
                        }
                        else{
                            // read ori is + -; contig ori should be + -
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                        }
                    }
                    else{
                        distance_fwd = fwd.pos;
                        if( rev.strand == 0 ){
                            // read ori is - +; contig ori should be - +
                            distance_rev = rev.pos;
                        }
                        else{
                            // read ori is - -; contig ori should be + +
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                        }
                    }
                    spam = distance_fwd + distance_rev;
                    return spam;
                case 2: //-> <- IN
                    //positive strand
                    if( fwd.strand == 0 ){
                        distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;//+
                        distance_rev = rev.pos;//-
                    }
                    else{
                        distance_fwd = fwd.pos; //-
                        distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos; //+
                    }
                    spam = distance_fwd + distance_rev;
                    return spam;
                case 1://<- -> OUT
                    // orientation is <-...-> (-,+)
                    if( fwd.strand == 0 ){
                        distance_fwd = fwd.pos;
                        distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                    }
                    else{
                        distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                        distance_rev = rev.pos;
                    }
                    spam = distance_fwd + distance_rev;
                    return spam;
                case 3: //<- <- REV
                    if( fwd.strand == 0 ){
                        distance_fwd = fwd.pos;
                        if( rev.strand == 0 ){
                            // read ori is + +, contig ori should be - -
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                        }
                        else{
                            // read ori is + -; contig ori should be - +
                            distance_rev = rev.pos;
                        }
                    }
                    else{
                        distance_fwd = kdb->get_seq_len(fwd.ctg) - fwd.pos;
                        if( rev.strand == 0 ){
                            // read ori is - +; contig ori should be + -
                            distance_rev = kdb->get_seq_len(rev.ctg) - rev.pos;
                        }
                        else{
                            // read ori is - -; contig ori should be + +
                            distance_rev = rev.pos;
                        }
                    }
                    spam = distance_fwd + distance_rev;
                    return spam;
                default:
                    cout << "error in orientation " << inferred_orientation << " unknow value" << endl;
                    exit(1);
            }

        }

        return int(spam);
    }

    vector<hit> compute_score(vector<hit> &tfwd, int seqlen){
        //we sort the array by ctg star
        sort(tfwd.begin(),tfwd.end(),compare_by_ctg_pos);
        hit mfwd=tfwd[0];
        int mc=0;
        vector<hit> index;
        int add_last=0;
        for (int i = 0; i < tfwd.size() ; i++) {
            //means that both kmers comes from the same region
            if((mfwd.ctg == tfwd[i].ctg) && (tfwd[i].pos-mfwd.pos <= seqlen)){
                mc++;
            }else{  //we add the last element to the array
                mfwd.score=mc;
                index.push_back(mfwd);
                //we add a new mfwd
                mc=1;
                mfwd=tfwd[i];
                mfwd.score=1;
                if(i == tfwd.size() - 1){
                    add_last=1;
                }
            }
        }
        //we add the last kmer to the set because is has an score of 1
        if(add_last){
            index.push_back(mfwd);
        }
        //we add the kmer to the set because it has score max
        if(mc == tfwd.size()){
            mfwd.score=mc;
            index.push_back(mfwd);
        }
        //we sort by score in order to print the output
        sort(index.begin(),index.end(),compare_score);
        return index;
    }

    hit rescue_mate(hit &bhit, vector<hit> &mates,bool fwdrev,int d){
        int distance=0; //assuming foward first
        hit tmp;
        //d=d-150; //is distance less the read length
        int min_spam=d - int(d*0.2);//for the moment we allow 20% from the observed distance
        if(min_spam < 0){
            min_spam=0;
        }
        int max_spam=d + int(d*0.2);//for the moment we allow 20% from the observed distance

        //we rescue only  if one pair satisfice the constrains of orientation and distance
        vector<int> rpairs;

        for (int i = 0; i <mates.size() ; i++) {
            distance=0;
            //best hit is in foward read
            if (fwdrev) {
                //to try to rescue we need at least 4 unique mers supporting the position
                if(mates[i].score > msr_long) {
                    distance = compute_distance_fragment2(bhit, mates[i]);
                }
            } else {//best hit is in reverse read
                if(mates[i].score > msr_long) {
                    distance = compute_distance_fragment2(mates[i], bhit);
                }
            }

            if (distance <= max_spam && distance >= min_spam) {
                rpairs.push_back(i);
            }
        }
        //rescue was sucess because we found only one position that satisfices orientation and distance
        if(rpairs.size() == 1){
            //we found only one posible pair to perform the rescue
            tmp=mates[rpairs[0]];
            //if within a contig we give an score of 12 = 60
            if(bhit.ctg == tmp.ctg) {
                tmp.score=22;
            }else{
            //we give the score of the rescue in this case
            tmp.score=21;//only to check if rescue is being efective
            }
            return tmp;
        }
        return tmp;
    }

//function to check the distance of a given pair of mers
 bool check_dist(hit &bhit, hit &ahit, int d){
        int distance=0; //assuming foward first
        hit tmp;
        //d=d-150; //distance less the read length 
        int min_spam=d - int(d*0.2);//for the moment we allow 20% from the observed distance
        if(min_spam < 0){
            min_spam=0;
        }
        int max_spam=d + int(d*0.2);//for the moment we allow 20% from the observed distance
	//we compute the distance beetween the pair of reads
       distance = compute_distance_fragment2(bhit, ahit);
	//we check if the distance live in the expected range of min-max
       if (distance <= max_spam && distance >= min_spam) {
			return true;
	  }else{
			return false;
	  }
 }
    //funtions for get advanced values
    int get_vhs(void){
        return vhs_long;
    }
    int get_ms(void){
        return ms_long;
    }
    int get_msr(void){
        return msr_long;
    }
    int get_ssr(void){
        return ssr_long;
    }

};


//thread function to map long partial corrected reads to a set of contigs using several inserts sizes
void * threaded_longmap (void* args)
{
    //casting of variables
    targ *tmp = (targ*) args;
    pthread_mutex_t *mutex = &tmp->mut;//a single mutex to control access to fastq file between threads
    pthread_mutex_t *pinfo = &tmp->pin;//to print information related to mapped pairs etc
    pthread_mutex_t *pinfom= tmp->pinm; //array of mutex to control writing to samfiles

    KmerDB *dbt=(KmerDB *)tmp->db;
    LongReads *lon=(LongReads *)tmp->lon;

    //compressed
    kstream<gzFile, FunctorZlib> *ks1= (kstream<gzFile, FunctorZlib> *) tmp->ks1;

    uint *samout=(uint *) tmp->samflag;
    ofstream *sam;
    uint *samp;
    uint *single;

    sam=(ofstream *) tmp->sam;//array of sams files
    samp=(uint *) tmp->samp;//flag that indicate if is a wirting in pair or single mode sam
    vector<int> inserts=lon->get_inserts();//vector of insert sizes
    vector<int> obs_inserts=lon->get_obs_inserts();

    //advanded options
    int vhs=lon->get_vhs();//vector size
    int ms=lon->get_ms(); //minimum score
    int msr=lon->get_msr(); //for future use
    int ssr=lon->get_ssr(); //length reads

    //init of locals variables
    uint number_seq=0;//variable to control the number of sequence that we read before the mapping
    uint kmer_size;
    kmer_size = (uint)dbt->getKmerSize();
    string qual(kmer_size,'I');
    hit indexs,indexf;
    vector<hit> fwdh;//local storage of hits
    vector<hit> revh;//local storage of hits
    //we start the mapping
    uint64_t tpairs=0;//total pairs mapped in thread
    uint64_t tmap=0;//total asked kmers
    uint64_t pmap=0;//total mapped-kmers pairs
    uint64_t tmapp=0;//total mapped-pairs
    uint64_t totalseqs=0;
    //hashing related variables
    uint64_t hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
    //we get both kmers in foward and reverse mers
    uint64_t hVal1 = 0, hVal2 = 0;
    int l=0,c=0;
    kseq fwd;
    vector<string> localseq;//storage of local seq
    bool dojob=1;//variable to control thread execution
    string kmersizeS=to_string(kmer_size);
    //we check the time spent
    double t_begin,t_end; struct timeval timet;
    gettimeofday(&timet, NULL); t_begin = timet.tv_sec +(timet.tv_usec/1000000.0);
    while(dojob) {
        //we try to lock the mutex to reads sequences from files
        number_seq = 0;
        //we got acess to the global fasta file
        pthread_mutex_lock(mutex);
        //we read 1mb reads from  fwd and rev fastq/fasta file keep reads names and assing an unique id
        while ((number_seq < 1000) && ((l = ks1->read(fwd)) >= 0) ) {
		//we continue reading if seq is shorter kmer_size
		if(fwd.seq.length() < kmer_size){
                	continue;
            	}
            //we store the long read in memory
            localseq.push_back(to_string(number_seq) + " " + fwd.seq + " " + fwd.name );
            number_seq++;
            totalseqs++;
        }
        pthread_mutex_unlock(mutex);

        //we ask if the need to map reads
        if (localseq.size() > 0) {
            dojob = 1;
        } else {
            dojob = 0;
        }
        //here we perform the mapping of the batch of sequences
        if(dojob){
            //we need to iterate the seq of sequence in each distance
            for (int i = 0; i <inserts.size() ; i++) {
                   int d=inserts[i];//this is the target insert size
                    //now we iter each sequence picking pairs of mers
                for (int j = 0; j < localseq.size(); ++j) {
                    string line=localseq[j];
                    istringstream iss(line);
                    vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());
                    //we check that the split was ok //0=id,1=seq,2=name
                    assert(tokens.size() == 3);
                    //we took two part of sequences at a given distance(d)
                    //we use sliding windows of 50bp after try the first map
                    int max_d = tokens[1].length() - d;
                    //int max_j = tokens[1].length() - d - 150;//max
                    //tpairs++;
                    //we skyp a sequence if is shorter than the length of d
                    if(max_d < 0 ){
                        continue;
                    }
                    //we need to split the long read in read pairs of length 200 whit a moving window of 100bp
                    for (int k = 0; k < max_d ; k+=100) {
                        string fwds=tokens[1].substr(k,ssr); //0+200 ssr fwd read
                        string revs=tokens[1].substr(k+d-ssr,ssr);//d-200 ssr rev read
                        //test write to samout
                        //sam[i] <<tokens[2]<<" "<<fwds<<" "<<revs<<" "<<k<<" "<<k+d-150<<" "<<d<<" "<<max_d<<" "<<tokens[1].length()<<endl;
                        //we attempt to map the pair of subreads
                        int max_i = fwds.length() - kmer_size;
                        int max_j = revs.length() - kmer_size;
                        tpairs++;
                        //we skyp pairs of reads shorter than the kmer length
                        if(max_i < 0 || max_j < 0){
                            continue;
                        }
                        //init hash values before mapping of kmers
                        hfwd1 = 0, hrev1 = 0, hfwd2 = 0, hrev2 = 0;
                        //we get both kmers in foward and reverse mers
                        hVal1 = 0, hVal2 = 0;
                        //init foward read
                        hVal1 = NTPC64(fwds.substr(0, kmer_size).c_str(), kmer_size, hfwd1, hrev1); // initial hash value for fwd
                        vector<hit> tfwd;
                        int max_h=0;
                        for (int f = 0; f <= max_i; f++) {
                            string fkmer = fwds.substr(f, kmer_size);
                            indexf = dbt->lookup_kmer(hfwd1, hrev1, fkmer);
                            tmap++;
                            if (indexf.map) {
                                tfwd.push_back(indexf);
                                max_h++;
                                if(max_h > vhs){
                                    break;
                                }

                            }
                            hVal1 = NTPC64(fwds[f], fwds[f + kmer_size], kmer_size, hfwd1, hrev1); //recursive hash k+L
                        }
                        //if we found at least 4 kmer in foward we try the reverse
                        if(tfwd.size() > 10) {
                            //we reversecomp the revs to produce a pair in -> <- orientation
                            revs=lon->revcomp(revs);
                            //init reverse read
                            hVal2 = NTPC64(revs.substr(0, kmer_size).c_str(), kmer_size, hfwd2,
                                           hrev2); //initial hash value for re
                            vector <hit> trev;
                            int max_h = 0;
                            for (int r = 0; r <= max_j; r++) {
                                string skmer = revs.substr(r, kmer_size);
                                indexs = dbt->lookup_kmer(hfwd2, hrev2, skmer);
                                tmap++;
                                if (indexs.map) {
                                    trev.push_back(indexs);
                                    pmap++;
                                    max_h++;
                                    if (max_h > vhs) {
                                        break;
                                    }
                                }
                                hVal2 = NTPC64(revs[r], revs[r + kmer_size], kmer_size, hfwd2, hrev2); //recursive hash k+L
                            }
                            //now we need to decide if we trusths in the positions
                            //we can improve just asking if we find at least 4 kmers
                            if(trev.size() > 10) {
                                vector<hit> tfwd2=lon->compute_score(tfwd,fwds.length());
                                vector<hit> trev2=lon->compute_score(trev,revs.length());
                                hit maxfwd=tfwd2[0];
                                hit maxrev=trev2[0];
                                //maxfwd.seqname=ill->get_clone_name(tokens[3]);
                                //maxrev.seqname=ill->get_clone_name(tokens[4]);

                                maxfwd.seqname=tokens[2]+"_"+to_string(d)+"_"+to_string(k);
                                maxrev.seqname=tokens[2]+"_"+to_string(d)+"_"+to_string(k);


                                hit rescue;
                                //if two reads are majority we don't touch them and we trush in their positions

                                //we will check if they are at the expected distance and orientation assumming a 10% standart dev
                                //if(maxfwd.score > 5 && maxrev.score > 5 && lon->check_dist(maxfwd,maxrev,obs_inserts[i])){
                                if(maxfwd.score > ms && maxrev.score > ms ){
                                    fwdh.push_back(maxfwd);
                                    revh.push_back(maxrev);
                                    tmapp++;
                                }else{
                                    //we rescue if at least one of each other is majority
                                    if(maxfwd.score > ms){
                                        rescue = lon->rescue_mate(maxfwd, trev2, 1, obs_inserts[i]);
                                        rescue.seqname = maxrev.seqname;
                                        if (rescue.score > 0) {
                                            fwdh.push_back(maxfwd);
                                            revh.push_back(rescue);
                                            tmapp++;
                                        }
                                    }else{
                                        if(maxrev.score > ms) {
                                            rescue = lon->rescue_mate(maxrev, tfwd2, 0, obs_inserts[i]);
                                            rescue.seqname = maxfwd.seqname;
                                            if (rescue.score > 0) {
                                                fwdh.push_back(rescue);
                                                revh.push_back(maxrev);
                                                tmapp++;
                                            }
                                        }
                                    }
                                }
                                //clear tmp memory for container
                                tfwd2.erase(tfwd2.begin(),tfwd2.end());
                                trev2.erase(trev2.begin(),trev2.end());
                            }
                            //clear tmp memory for container
                            tfwd.erase(tfwd.begin(),tfwd.end());
                            trev.erase(trev.begin(),trev.end());

                        }

                    }

                }
                //here we need to write the results to samfiles
                //we print the result in sam  format using singlest or paired reads, it depends of the scaffolder
                pthread_mutex_lock(&pinfom[i]);
                //means write resuls as single end reads
                if(*samout == 1){
                    //we need to write the results to the bam file if is needed
                    //we print reads as single end reads for BOSS/SCAFFMATCH
                    for (int jj = 0; jj < fwdh.size(); ++jj) {
                        string fsam_flag = "0";
                        string rsam_flag = "0";
                        //if the kmer hit the reverse strand I need to ouput the reverse kmer in the samfile
                        if (fwdh[jj].strand) {
                            fsam_flag = "16";
                        }
                        if (revh[jj].strand) {
                            rsam_flag = "16";
                        }
                        string tmp1 =
                                fwdh[jj].seqname + "\t" + fsam_flag + "\t" + dbt->get_seq_name(fwdh[jj].ctg) + "\t" +
                                to_string(fwdh[jj].pos) + "\t" + to_string(fwdh[jj].score * 10) + "\t" +
                                kmersizeS + "M\t*\t0\t0\t" + fwdh[jj].kseqm + "\t" + qual +
                                "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        string tmp2 =
                                revh[jj].seqname + "\t" + rsam_flag + "\t" + dbt->get_seq_name(revh[jj].ctg) + "\t" +
                                to_string(revh[jj].pos) + "\t" + to_string(revh[jj].score * 10) + "\t" +
                                kmersizeS + "M\t*\t0\t0\t" + revh[jj].kseqm + "\t" + qual +
                                "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";

                        (*samp)++;//we increment the pointer of number of mapped kmers pairs
                        sam[i].write(tmp1.c_str(), tmp1.length());
                        sam[i].write(tmp2.c_str(), tmp2.length());
                    }
                }else{
                    //we print reads as single end reads pair-ends for BESST/OPERA scaffolders
                    for (int jj = 0; jj < fwdh.size(); ++jj) {
                        int fsam_flag = 0;
                        int rsam_flag = 0;
                        //if the kmer hit the reverse strand I need to ouput the reverse kmer in the samfile
                        if (fwdh[jj].strand) {
                            fsam_flag += 16;
                            rsam_flag+=32;
                        }
                        if (revh[jj].strand) {
                            rsam_flag = 16;
                            fsam_flag+=32;
                        }

                        fsam_flag+=65;//paired and read first in pair
                        rsam_flag+=129;//paired and read second in pair

                        string tmp1="";
                        string tmp2="";
                        //if within contigs we can compute the distance
                        if(fwdh[jj].ctg == revh[jj].ctg ) {
                            int distance1=(int)(revh[jj].pos - fwdh[jj].pos);
                            tmp1 = fwdh[jj].seqname + "\t" + to_string(fsam_flag) + "\t" + dbt->get_seq_name(fwdh[jj].ctg) + "\t" +
                                   to_string(fwdh[jj].pos) + "\t" + to_string(fwdh[jj].score * 10) + "\t" +
                                   kmersizeS + "M\t=\t"+to_string(revh[jj].pos)+"\t"+to_string(distance1)+"\t" + fwdh[jj].kseqm + "\t" + qual +
                                   "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                            int distance2=(int)(fwdh[jj].pos - revh[jj].pos);
                            tmp2 = revh[jj].seqname + "\t" + to_string(rsam_flag) + "\t" + dbt->get_seq_name(revh[jj].ctg) + "\t" +
                                   to_string(revh[jj].pos) + "\t" + to_string(revh[jj].score * 10) + "\t" +
                                   kmersizeS + "M\t=\t"+to_string(fwdh[jj].pos)+"\t"+to_string(distance2)+"\t" + revh[jj].kseqm + "\t" + qual +
                                   "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        }else{

                            tmp1 = fwdh[jj].seqname + "\t" + to_string(fsam_flag) + "\t" + dbt->get_seq_name(fwdh[jj].ctg) + "\t" +
                                   to_string(fwdh[jj].pos) + "\t" + to_string(fwdh[jj].score * 10) + "\t" +
                                   kmersizeS + "M\t"+dbt->get_seq_name(revh[jj].ctg)+"\t"+to_string(revh[jj].pos)+"\t0\t" + fwdh[jj].kseqm + "\t" + qual +
                                   "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";

                            tmp2 = revh[jj].seqname + "\t" + to_string(rsam_flag) + "\t" + dbt->get_seq_name(revh[jj].ctg) + "\t" +
                                   to_string(revh[jj].pos) + "\t" + to_string(revh[jj].score * 10) + "\t" +
                                   kmersizeS + "M\t"+dbt->get_seq_name(fwdh[jj].ctg)+"\t"+to_string(fwdh[jj].pos)+"\t0\t" + revh[jj].kseqm + "\t" + qual +
                                   "\tAS:i:0\tXN:i:0\tXM:i:0\tXO:i:0\tXG:i:0\tNM:i:0\tMD:Z:"+kmersizeS+"\tYT:Z:UU\n";
                        }

                        (*samp)++;//we increment the pointer of number of mapped kmers pairs
                        sam[i].write(tmp1.c_str(), tmp1.length());
                        sam[i].write(tmp2.c_str(), tmp2.length());
                    }
                }

                pthread_mutex_unlock(&pinfom[i]);
                //we need to clear the local variables of the thread before the next execution
                //clear mapped hits
                fwdh.erase(fwdh.begin(),fwdh.end());
                revh.erase(revh.begin(),revh.end());

            }
            //we need to erase the file the container of long-reads
            //we clear the localseq variable
            localseq.erase(localseq.begin(),localseq.end());
        }

    }
    gettimeofday(&timet, NULL); t_end = timet.tv_sec +(timet.tv_usec/1000000.0);
    double elapsed3 = t_end - t_begin;
    //we print some stats from thread information
    pthread_mutex_lock(pinfo);
    cout << "End of mapping in thread\n";
    printf("Query of %llu kmer-pairs  in %.2fs rate %.2f per second\n", tmap,elapsed3,tmap/elapsed3);
    printf("Query of %llu pairs  in %.2fs\n", tpairs,elapsed3);
    printf("Number of total mapped kmer-pairs %llu\n", pmap);
    printf("Number of total mapped reads-pairs %llu\n", tmapp);
    printf("Number of total long reads %llu\n", totalseqs);
    pthread_mutex_unlock(pinfo);

    
    return NULL;
}





int main (int argc, char* argv[]){
    if(argc !=13 ){
        printf("Usage :\n");
          //10 6 4 20 15 4 200
        printf("%s <uniq_kmers> <rfasta> <kmer_size> <ncpu> <libs> <advanced illumina> <advanced long>\n",argv[0]);
        return EXIT_FAILURE;
    }
    string uniqkmers = argv[1];
    string rfasta=argv[2];
    uint kmersize=(uint)atoi(argv[3]);
    uint cpus=(uint)atoi(argv[4]);
    string libfile=argv[5];
    //advanced short values
    int vhs_ill=atoi(argv[6]);
    int ms_ill=atoi(argv[7]);
    int msr_ill=atoi(argv[8]);
    //advanced long values
    int vhs_long=atoi(argv[9]);
    int ms_long=atoi(argv[10]);
    int msr_long=atoi(argv[11]);
    int ssr_long=atoi(argv[12]);

    string line="";
    //kmerDB we use  a perfect hash function for fast kmer query
    KmerDB *b = new KmerDB(uniqkmers,kmersize,rfasta,cpus);
    //b->dump_database();
    //read pair of reads including a prefix separated by space
    //start illumina map
    ifstream illuminalibs(libfile);
    //string line="";
    Illumina *c;//pointer to the illumina object
    LongReads *d;//we create the long read object
    vector<int> inserts; //variable to store the insert-sizes

    while(getline(illuminalibs,line)) {
        istringstream iss(line);
        vector<string> lib((istream_iterator<string>(iss)), istream_iterator<string>());
        //short library
        if(lib[0].compare("short") == 0) {
            //we create the illumina object to do the mapping of the libs
            c = new Illumina(lib[1], lib[2], lib[3], vhs_ill, ms_ill, msr_ill, b);
            //1 means singlest and other value means paired report of sam
            c->mapshortreads((uint) atoi(lib[4].c_str()), cpus);
        }else{
            //is a long library
            istringstream iss2(lib[3]);
            string token;
            //we get the inserts sizes
            while(getline(iss2, token, ',')) {
                inserts.push_back(atoi(token.c_str()));
            }
            //we create the illumina object to do the mapping of the libs
            d = new LongReads(lib[1],lib[2], vhs_long, ms_long, msr_long, ssr_long, inserts,b);
            d->maplongreads((uint)atoi(lib[4].c_str()),cpus);
            inserts.erase(inserts.begin(),inserts.end());//we clean the insert-size vector
        }
    }
    //end

}
