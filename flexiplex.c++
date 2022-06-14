// Copyright 2022 Nadia Davidson 
// This program is distributed under the MIT License.
// We also ask that you cite this software in publications
// where you made use of it for any part of the data analysis.

#include <iostream>
#include <fstream>
#include <istream>
#include <string>
#include <sstream>
#include <unordered_set>
#include <unordered_map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>

#include "edlib.h"

//#include <gperftools/profiler.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << "usage: flexiplex [options] reads_input"  << endl;
  cerr << "  reads_input: a .fastq or .fasta file" << endl;
  cerr << "  options: " << endl;
  cerr << "     -k known_list  Text file of expected barcodes (one row per barcode)." << endl; 
  cerr << "                    if not provided, flexiplex will accept any sequence as a barcodes. " << endl; 
  cerr << "     -p primer  Left flank sequence to search for (default: CTACACGACGCTCTTCCGATCT)" << endl;
  cerr << "     -T polyT   Right flank sequence to search for (default: TTTTTTTTT)" << endl;
  cerr << "     -b N   Barcode length (default: 16)" << endl;
  cerr << "     -u N   UMI length (default: 12)" << endl;
  cerr << "     -e N   Maximum edit distance to barcode (default 2)" << endl;
  cerr << "     -f N   Maximum edit distance to primer+ployT (default 10)" << endl;
  cerr << "     -r true/false   Whether to remove barcode sequence from reads" << endl;
  cerr << endl;
}


// compliment nucleotides - used to reverse compliment string
char compliment(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//Inplace reverse compliment
void reverse_compliment(string & seq){
   reverse(seq.begin(),seq.end());
   transform(seq.begin(),seq.end(),seq.begin(),compliment);
}

//Holds the search string patterns
struct SearchSeq {
  string primer;
  string polyA;
  string umi_seq;
  string temp_barcode;
} ;

//Holds the found barcode and assocaited information 
struct Barcode {
  string barcode;
  string umi;
  int editd;
  int next_editd;
  int flank_editd;
  int flank_start;
  int flank_end;
} ;


//Given a string 'seq' search for substring with primer and polyT sequence followed by
//a targeted search in the region for barcode
//Seaquence seearch is performed using edlib
Barcode get_barcode(string & seq,
		    unordered_set<string> *known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd,
		    SearchSeq & ss){
  
  //initialise struct variables for return:
  Barcode barcode;
  barcode.editd=100; barcode.next_editd=100; barcode.flank_editd=100;

  //initialise edlib configuration
  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 5};
  string search_string; 
  EdlibAlignResult result; 

  //search for primer and ployT (barcode and umi as wildcards)
  search_string=ss.primer+ss.temp_barcode+ss.umi_seq+ss.polyA;
  result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK | result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found - return
  } //fill in info about found primer and polyT location
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=result.endLocations[0];
  
  if(known_barcodes->size()==0){ //if not checking against known list return sequence after the primer
    barcode.barcode=seq.substr(result.startLocations[0]+ss.primer.size(),ss.temp_barcode.length());
    barcode.editd=0;
    edlibFreeAlignResult(result);
    return(barcode);
  }

  //otherwise, check again a list of known barcodes
  int OFFSET=5; //include this many bases extra to the left and right incase of indels errors
  int BC_START=result.startLocations[0]+ss.primer.length()-OFFSET;
  string barcode_seq=seq.substr(BC_START,ss.temp_barcode.length()+2*OFFSET);
  //reconfigure edlib for the targetted search
  edlibConf = edlibNewAlignConfig(barcode_max_editd,EDLIB_MODE_HW,EDLIB_TASK_LOC,NULL,0);
  unordered_set<string>::iterator known_barcodes_itr=known_barcodes->begin();
  for(; known_barcodes_itr!=known_barcodes->end(); known_barcodes_itr++){
    search_string=(*known_barcodes_itr);
    result = edlibAlign(search_string.c_str(), search_string.length(), barcode_seq.c_str(), barcode_seq.length(),edlibConf);
    if(result.status == EDLIB_STATUS_OK && result.numLocations!=0 && (result.editDistance < barcode.editd)) {
      barcode.next_editd=barcode.editd; //if this barcode has the lowest edit distance update best match..
      barcode.editd=result.editDistance;
      barcode.barcode=*known_barcodes_itr;
      barcode.umi=seq.substr(BC_START+result.endLocations[0]+1,ss.umi_seq.length()); //get the umi
      if(barcode.editd==0){ //perfect match, stop looking and return
	edlibFreeAlignResult(result);
	return(barcode); 
      }
    }
  }
  
  edlibFreeAlignResult(result);
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
vector<Barcode> big_barcode_search(string & sequence, unordered_set<string> & known_barcodes,
				   int max_flank_editd, int max_editd, SearchSeq & ss){
  vector<Barcode> return_vec; //vector of all the barcodes found

  //search for barcode
  Barcode result=get_barcode(sequence,&known_barcodes,max_flank_editd,max_editd,ss);
  if(result.editd<=max_editd) //add to return vector if edit distance small enough
    return_vec.push_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if(return_vec.size()>0){
    string masked_sequence = sequence;
    for(int i=0; i<return_vec.size(); i++){
      int flank_length=return_vec.at(i).flank_end-return_vec.at(i).flank_start;
      masked_sequence.replace(return_vec.at(i).flank_start,flank_length,string(flank_length,'X'));
    } //recursively call this function until no more barcodes are found
    vector<Barcode> masked_res;
    masked_res=big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd,ss);
    return_vec.insert(return_vec.end(),masked_res.begin(),masked_res.end()); //add to result
  }
    
  return(return_vec);
}

// utility function to check true/false input options
bool get_bool_opt_arg(string value){
  transform(value.begin(), value.end(), value.begin(), ::tolower);
  if( value.compare("true")==0 | value.compare("t")==0 | value.compare("1")==0){
    return true;
  } else if (value.compare("false")!=0 | value.compare("f")!=0 | value.compare("0")!=0){
    return false;
  } else {
    cerr << "Unknown argument to boolean option" << endl;
    print_usage();
    exit(1);
  } 
}

// print information about barcodes
void print_stats(string read_id, vector<Barcode> & vec_bc, ostream & out_stream){
  for(int b=0; b<vec_bc.size() ; b++){
    out_stream << read_id << "\t"
	       << vec_bc.at(b).barcode << "\t"
	       << vec_bc.at(b).flank_editd << "\t"
	       << vec_bc.at(b).editd << "\t"
	       << vec_bc.at(b).next_editd << "\t"
	       << endl;
  }
}

//print fastq or fasta lines..
void print_read(string read_id, string read, string qual, vector<Barcode> & vec_bc, ostream & out_stream){
  //flag for read format
  bool is_fastq=!(qual==""); //no quality scores passed = fasta

  //loop over the barcodes found... usually will just be one
  for(int b=0; b<vec_bc.size() ; b++){

    //format the new read id. Using FLAMES format.
    stringstream ss;
    ss << (b+1) << "of" << vec_bc.size() ;
    string new_read_id=vec_bc.at(b).barcode+"_"+vec_bc.at(b).umi+"#"+read_id+ss.str();
    
    // work out the end base in case multiple barcodes
    int read_length=read.length()-vec_bc.at(b).flank_end;
    for(int f=0; f<vec_bc.size(); f++){
      int temp_read_length=vec_bc.at(f).flank_start-vec_bc.at(b).flank_end;
      if(temp_read_length>0 && temp_read_length<read_length)
	read_length=temp_read_length;
    }
    
    //output to the new read file
    if(is_fastq)
      out_stream << "@" << new_read_id << endl;
    else
      out_stream << ">" << new_read_id << endl;
    out_stream << read.substr(vec_bc.at(b).flank_end,read_length) << endl;
    if(is_fastq){
      out_stream << "+" << new_read_id << endl;
      out_stream << qual.substr(vec_bc.at(b).flank_end,read_length) << endl;
    }
  }
}


// MAIN is here!!
int main(int argc, char **argv){

  //Variables to store user options
  //Set these to their defaults
  int expected_cells=0; //(d)
  int edit_distance=2; //(e)
  int flank_edit_distance=10; //(f)
  string out_bc_filename=""; //(o)
  string out_stat_filename=""; //(t)
  bool split_file_by_barcode=false; //(s)
  bool remove_barcodes=true; //(r)
  
  SearchSeq search_patterns;
  search_patterns.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_patterns.polyA = string(9,'T'); //(T)
  search_patterns.umi_seq = string(12,'?'); //(length u)
  search_patterns.temp_barcode = string(16,'?'); //(length b)

  //Set of known barcodes 
  unordered_set<string> known_barcodes;

  /*** Pass command line option *****/
  int c;
  int params=1;
  ifstream file;
  string line;

  while((c =  getopt(argc, argv, "k:e:d:s:o:t:r:f:p:T:u:b:")) != EOF){
    switch(c){
    case 'k': { //k=list of known barcodes
      string file_name(optarg);
      cerr << "Reading barcode file "<< file_name << endl;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      if(!(file.good())){
	cerr << "Unable to open file: " << file_name << endl;
	exit(1);
      }
      while ( getline (file,line) ){
	string bc;
	istringstream line_stream(line);
	line_stream >> bc;
	known_barcodes.insert(bc); 
      }
      cerr << "Barcodes read: " << known_barcodes.size() << endl;
      file.close();
      cerr << "Done reading barcode file " << file_name << endl;
      params+=2;
      break;
    }
    case 'r':{
      remove_barcodes=get_bool_opt_arg(optarg);
      params+=2;
      break;
    }
    case 'e':{
      edit_distance=atoi(optarg);
      cerr << "Setting max barcode edit distance to "<< edit_distance << endl;
      params+=2;
      break;
    }
    case 'f':{
      flank_edit_distance=atoi(optarg);
      cerr << "Setting max flanking sequence edit distance to "<< flank_edit_distance << endl;
      params+=2;
      break;
    }
    case 'p':{
      search_patterns.primer=optarg;
      cerr << "Setting primer to search for: " << search_patterns.primer << endl;
      params+=2;
      break;
    }
    case 'T':{
      search_patterns.polyA=optarg;
      cerr << "Setting polyT to search for: " << search_patterns.polyA << endl;
      params+=2;
      break;
    }
    case 'u':{
      int ul=atoi(optarg);
      search_patterns.umi_seq=string(ul,'?');
      cerr << "Setting UMI length to " << ul << endl;
      params+=2;
      break;
    }
    case 'b':{
      int bl=atoi(optarg);
      search_patterns.temp_barcode=string(bl,'?');
      cerr << "Setting barcode length to " << bl << endl;
      params+=2;
      break;
    }
    case '?': //other unknown options
      cerr << "Unknown option.. stopping" << endl;
      print_usage();
      exit(1);
      break; 
    }
  }

  istream * in;
  //check that a read file is given
  if(params>=argc){
    cerr << "No filename given... getting reads from stdin..." << endl;
    in=&std::cin;
  } else {
    // check that the reads fileis okay
    string reads_file=argv[params];
    file.open(reads_file);
    if(!(file.good())){
      cerr << "Unable to open file " << reads_file << endl;
      print_usage();
      exit(1);
    }
    in=&file;
  }
    
  //set the output filenames
  out_stat_filename="flexiplex-reads_barcodes.txt";
  out_bc_filename="flexiplex-barcodes_counts.txt"; 

  //  ProfilerStart("my.prof");
  
  /********* FIND BARCODE IN READS ********/
  string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;
  
  ofstream out_stat_file;
  if(known_barcodes.size()>0){
    out_stat_file.open(out_stat_filename);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tNextBestBarcodeEditDist"<<endl;
  }
  cerr << "Searching for barcodes..." << endl;
  bool is_fastq=true;
  unordered_map< string, int > barcode_counts; 
  while ( getline (*in,line) ){
    istringstream line_stream(line);
    string read_id;
    line_stream >> read_id;
    int lines_to_skip=0;
    if (read_id[0] == '>'){ is_fastq=false;
    } else if (read_id[0] == '@'){ //fasta
    } else {
      cerr << "Unknown read format... exiting" << endl; exit(1);
    }
    read_id.erase(0,1);
    if(getline (*in,line)){
      r_count++; //progress counter
      if(r_count % 100000 == 0)
	cerr << r_count/1000 << " thousand reads processed.." << endl;
      
      //forward search
      vector<Barcode> vec_bc_for=big_barcode_search(line,known_barcodes,
						    flank_edit_distance,edit_distance,search_patterns);
      for(int b=0; b<vec_bc_for.size(); b++)
	barcode_counts[vec_bc_for.at(b).barcode]++;
      string rev_line=line;
      reverse_compliment(rev_line);
      //Check the reverse compliment of the read
      vector<Barcode> vec_bc_rev=big_barcode_search(rev_line,known_barcodes,
						    flank_edit_distance,edit_distance,search_patterns);
      for(int b=0; b<vec_bc_rev.size(); b++)
	barcode_counts[vec_bc_rev.at(b).barcode]++;

      if((vec_bc_for.size()+vec_bc_rev.size())>0)
	bc_count++;
      if((vec_bc_for.size()+vec_bc_rev.size())>1 ){
	multi_bc_count++;
      }

      if(known_barcodes.size()==0)
	continue; // if we are just looking for all possible barcodes don't output reads etc.
      
      print_stats(read_id, vec_bc_for, out_stat_file);
      print_stats(read_id, vec_bc_rev, out_stat_file);

      string qual_scores="";
      if(is_fastq)
	for(int s=0; is_fastq && s<2; s++) getline (file,qual_scores);
      print_read(read_id+"_+",line,qual_scores,vec_bc_for, std::cout);
      reverse(qual_scores.begin(),qual_scores.end());
      print_read(read_id+"_-",rev_line,qual_scores,vec_bc_rev, std::cout);
    }
  }
  file.close();

  cerr << "Number of reads processed: " << r_count << endl;
  cerr << "Number of reads where a barcode was found: " << bc_count << endl;
  cerr << "Number of reads where more than one barcode was found: " << multi_bc_count << endl;
  cerr << "All done!" << endl;

  if(known_barcodes.size()>0){
    out_stat_file.close();
    return(0);
  }
  
  if(barcode_counts.size()==0)
    return(0);
  
  typedef std::pair<std::string, int> pair;
  vector<pair> bc_vec;
  copy(barcode_counts.begin(),barcode_counts.end(), back_inserter<vector<pair>>(bc_vec));
  sort(bc_vec.begin(), bc_vec.end(),[](const pair &l, const pair &r){
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });
  vector<int> hist(bc_vec[0].second);
  ofstream out_bc_file;
  out_bc_file.open(out_bc_filename);
  for (auto const &bc_pair: bc_vec){
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << endl;
    hist[bc_pair.second-1]++;
  }
  out_bc_file.close();

  cout << "Reads\tBarcodes" << endl;
  for(int i=hist.size()-1; i>=0; i--)
    cout << i+1 << "\t" << hist[i] << endl;
    
}
