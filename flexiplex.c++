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

const static string VERSION="0.9";

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << "usage: flexiplex [options] [reads_input]"  << endl;
  cerr << "  reads_input: a .fastq or .fasta file. Will read from stdin if empty." << endl;
  cerr << "  options: " << endl;
  cerr << "     -k known_list   Either 1) a text file of expected barcodes in the first column," << endl; 
  cerr << "                     one row per barcode, or 2) acomma separate string of barcodes. " << endl;
  cerr << "                     Without this option, flexiplex will search and report possible barcodes." << endl;
  cerr << "                     The generated list can be used for known_list in subsequent runs." << endl; 
  cerr << "     -r true/false   Replace read ID with barcodes+UMI, remove search strings" << endl;
  cerr << "                     including flanking sequenence and split read if multiple" << endl;
  cerr << "                     barcodes found (default: true)." << endl;
  cerr << "     -s true/false   Sort reads into separate files by barcode (default: false)" << endl;
  cerr << "     -p primer   Left flank sequence to search for (default: CTACACGACGCTCTTCCGATCT)." << endl;
  cerr << "     -T polyT    Right flank sequence to search for (default: TTTTTTTTT)." << endl;
  cerr << "     -n prefix   Prefix for output filenames." << endl;
  cerr << "     -b N   Barcode length (default: 16)." << endl;
  cerr << "     -u N   UMI length (default: 12)." << endl;
  cerr << "     -e N   Maximum edit distance to barcode (default 2)." << endl;
  cerr << "     -f N   Maximum edit distance to primer+ployT (default 10)." << endl;
  cerr << "     -h     Print this usage information." << endl;
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

//Holds the found barcode and associated information 
struct Barcode {
  string barcode;
  string umi;
  int editd;
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
  barcode.editd=100; barcode.flank_editd=100;

  //initialise edlib configuration
  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_LOC, additionalEqualities, 5};

  //search for primer and ployT (barcode and umi as wildcards)
  string search_string=ss.primer+ss.temp_barcode+ss.umi_seq+ss.polyA;
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
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
	       << vec_bc.at(b).umi << "\t"
	       << endl;
  }
}

void print_line(string id, string read, string quals, ostream & out_stream){

  //flag for read format
  bool is_fastq=!(quals==""); //no quality scores passed = fasta

  //output to the new read file
    if(is_fastq)
      out_stream << "@" << id << endl;
    else
      out_stream << ">" << id << endl;
    out_stream << read << endl;
    if(is_fastq){
      out_stream << "+" << id << endl;
      out_stream << quals << endl;
    }
}

//print fastq or fasta lines..
void print_read(string read_id,string read, string qual,
		vector<Barcode> & vec_bc, string prefix,
		bool split, unordered_set<string> & found_barcodes,
		bool trim_barcodes){
  //loop over the barcodes found... usually will just be one
  for(int b=0; b<vec_bc.size() ; b++){
    
    //format the new read id. Using FLAMES format.
    stringstream ss;
    ss << (b+1) << "of" << vec_bc.size() ;
    string barcode=vec_bc.at(b).barcode;
    string new_read_id=barcode+"_"+vec_bc.at(b).umi+"#"+read_id+ss.str();
    
    // work out the start and end base in case multiple barcodes
    int read_start=vec_bc.at(b).flank_end;
    int read_length=read.length()-read_start;
    for(int f=0; f<vec_bc.size(); f++){
      int temp_read_length=vec_bc.at(f).flank_start-read_start;
      if(temp_read_length>0 && temp_read_length<read_length)
	read_length=temp_read_length;
    }
    string qual_new=""; //don't trim the quality scores if it's a fasta file
    if(qual!="") qual_new=qual.substr(read_start,read_length);
    string read_new=read.substr(read_start,read_length);

    if(b==0 && !trim_barcodes){ //override if read shouldn't be cut
      new_read_id=read_id;
      read_new=read;
      qual_new=qual;
      b=vec_bc.size(); //force loop to exit after this iteration
    }
    
    if(split){ //to a file if spliting by barcode
      string outname=prefix+"_"+barcode+".";
      if(qual=="") outname+="fasta"; else outname+="fastq";
      ofstream outstream;
      if(found_barcodes.insert(barcode).second)
	outstream.open(outname); //override file if this is the first read for the barcode
      else
	outstream.open(outname,ofstream::app); //remove file if this is the first read for the barcode
      print_line(new_read_id,read_new,qual_new,outstream);
      outstream.close();
    } else {
      print_line(new_read_id,read_new,qual_new,std::cout);
    }
  }
}


// MAIN is here!!
int main(int argc, char **argv){

  cerr << "FLEXIPLEX " << VERSION << endl;

  //Variables to store user options
  //Set these to their defaults
  int expected_cells=0; //(d)
  int edit_distance=2; //(e)
  int flank_edit_distance=10; //(f)

  //set the output filenames
  string out_stat_filename="reads_barcodes.txt";
  string out_bc_filename="barcodes_counts.txt";
  string out_filename_prefix="flexiplex"; //(n)

  bool split_file_by_barcode=false; //(s)
  bool remove_barcodes=true; //(r)
  
  SearchSeq search_patterns;
  search_patterns.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_patterns.polyA = string(9,'T'); //(T)
  search_patterns.umi_seq = string(12,'?'); //(length u)
  search_patterns.temp_barcode = string(16,'?'); //(length b)
  
  //Set of known barcodes 
  unordered_set<string> known_barcodes;
  unordered_set<string> found_barcodes;

  /*** Pass command line option *****/
  int c;
  int params=1;
  ifstream file;
  string line;

  while((c =  getopt(argc, argv, "k:r:p:T:b:u:e:f:n:s:h")) != EOF){
    switch(c){
    case 'k': { //k=list of known barcodes
      string file_name(optarg);
      string bc;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      cerr << "Setting known barcodes from "<< file_name << endl;
      if(!(file.good())){ //if the string given isn't a file
	stringstream bc_list(file_name); string s;
	while (getline(bc_list, bc, ',')) //tokenize
	  known_barcodes.insert(bc);
      } else {
	// otherwise get the barcodes from the file..
	while ( getline (file,line) ){
	  istringstream line_stream(line);
	  line_stream >> bc;
	  known_barcodes.insert(bc); 
	}
	file.close();
      }
      cerr << "Number of known barcodes: " << known_barcodes.size() << endl;
      if(known_barcodes.size()==0){
	print_usage();
	exit(1); //case barcode file is empty
      }
      //set barcode length automatically from known barcodes..
      int bl=(known_barcodes.begin())->length();
      search_patterns.temp_barcode=string(bl,'?');
      cerr << "Setting barcode length automatically to " << bl << endl;
      params+=2;
      break;     
    }
    case 'r':{
      remove_barcodes=get_bool_opt_arg(optarg);
      cerr << "Setting read IDs to be replaced: "<< remove_barcodes << endl;
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
    case 'h':{
      print_usage();
      exit(1);
    }
    case 'n':{
      out_filename_prefix=optarg;
      cerr << "Setting output filename prefix to: " << out_filename_prefix << endl;
      params+=2;
      break;
    }
    case 's':{
      split_file_by_barcode=get_bool_opt_arg(optarg);
      cerr << "Split read output into separate files by barcode: " << split_file_by_barcode << endl;
      int max_split_bc=50;
      if(known_barcodes.size()>max_split_bc){
	cerr << "Too many barcodes to split into separate files: "<< known_barcodes.size()
	     << "> "<< max_split_bc<< endl;
	split_file_by_barcode=false;
      }
      params+=2;
      break;
    }
    case '?': //other unknown options
      cerr << "Unknown option.. stopping" << endl;
      print_usage();
      exit(1);
    }
  }

  cerr << "For usage information type: flexiplex -h" << endl;
  
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
    
  //  ProfilerStart("my.prof");
  
  /********* FIND BARCODE IN READS ********/
  string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;
  
  ofstream out_stat_file;
  out_stat_filename=out_filename_prefix+"_"+out_stat_filename;
  out_bc_filename=out_filename_prefix+"_"+out_bc_filename;
  params+=2;

  if(known_barcodes.size()>0){
    out_stat_file.open(out_stat_filename);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI"<<endl;
  }
  cerr << "Searching for barcodes..." << endl;
  bool is_fastq=true;
  unordered_map< string, int > barcode_counts; 
  string read_id_line;
  if(getline (*in,read_id_line)){ //check the first line for file type
    if(read_id_line[0]=='>'){ is_fastq=false;
    } else if (read_id_line[0] == '@'){ //fasta
    } else {
      cerr << "Unknown read format... exiting" << endl; exit(1);
    }
  }
  while ( getline (*in,line) ){
    string read_id;
    istringstream line_stream(read_id_line);
    line_stream >> read_id;
    read_id.erase(0,1);
    
    string qual_scores="";
    if(!is_fastq){ //fastq (account for multi-lines per read)
      string buffer; 
      while(getline(*in,buffer) && buffer[0]!='>')
	line+=buffer;
      read_id_line=buffer;
    } else { //fastq (get quality scores)
      for(int s=0; s<2; s++) getline(*in,qual_scores);
      getline(*in,read_id_line);
    }
    
    r_count++; //progress counter
    if(r_count % 100000 == 0)
      cerr << r_count/((double) 1000000 ) << " million reads processed.." << endl;
    
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
    
    print_read(read_id+"_+",line,qual_scores,vec_bc_for,out_filename_prefix,
	       split_file_by_barcode,found_barcodes,remove_barcodes);
    reverse(qual_scores.begin(),qual_scores.end());
    if(remove_barcodes || vec_bc_for.size()==0) //case we just want to print read once if multiple bc found.
      print_read(read_id+"_-",rev_line,qual_scores,vec_bc_rev,out_filename_prefix,
		 split_file_by_barcode,found_barcodes,remove_barcodes);
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
