// Copyright 2021 Nadia Davidson 
// This program is distributed under the GNU
// General Public License. We also ask that you cite this software in
// publications where you made use of it for any part of the data
// analysis.

/** 
 ** 
 **
 ** Author: Nadia Davidson
 ** Modified: 
 **/ 

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

#include "edlib.h"

#include <gperftools/profiler.h>

using namespace std;

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << "usage: flexiplex [options] reads_input"  << endl;
  cerr << "  reads_input: a .fastq or .fasta file" << endl;
  cerr << "  options: " << endl;
  cerr << "     -w white_list Text file of expected barcdes (one row per barcode)" << endl; 
  cerr << "     -p primer  Primer sequence to search for (default: CTACACGACGCTCTTCCGATCT)" << endl;
  cerr << "     -T polyT   polyT sequence to search for (default: TTTTTTTTT)" << endl;
  cerr << "     -b N   Barcode length (default: 16)" << endl;
  cerr << "     -u N   UMI length (default: 12)" << endl;
  cerr << "     -d N   Fetch barcode list from reads_input, assuming approx. this many cells (default 0 - off)" << endl;
  cerr << "     -e N   Maximum edit distance to barcode (default 2)" << endl;
  cerr << "     -f N   Maximum edit distance to primer+ployT (default 10)" << endl;
  //  cerr << "     -t stat_output  Name of file to output table of barcode-read assignments" << endl;
  //  cerr << "     -o read_output  Name of the file to output adjusted reads" << endl;
  cerr << "     -s true/false   Whether to split reads up when multiple barcodes seen" << endl;
  cerr << "     -r true/false   Whether to remove barcode sequence from reads" << endl;
  cerr << endl;
}


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

struct SearchSeq {
  string primer;
  string polyA;
  string umi_seq;
  string temp_barcode;
} ;

struct Barcode {
  string barcode;
  string umi;
  int editd;
  int next_editd;
  int flank_editd;
  int flank_start;
  int flank_end;
} ;


//search for substring where barcode and UMI should be located using edlib
Barcode get_barcode(string & seq,
		    unordered_set<string> *known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd,
		    SearchSeq & ss){
  
  //initialise struct variables:
  Barcode barcode;
  barcode.editd=100; barcode.next_editd=100; barcode.flank_editd=100;

  EdlibEqualityPair additionalEqualities[5] = {{'?','A'},{'?','C'},{'?','G'},{'?','T'},{'?','N'}};
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 5};
  string search_string; 
  EdlibAlignResult result; 

  // search for primer and ployT
  search_string=ss.primer+ss.temp_barcode+ss.umi_seq+ss.polyA;
  result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK | result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found
  }
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=result.endLocations[0];
  
  if(known_barcodes==NULL){ //if not checking against known list return..
    barcode.barcode=seq.substr(result.startLocations[0]+ss.primer.size(),ss.temp_barcode.length());
    edlibFreeAlignResult(result);
    return(barcode);
  }

  //result = edlibAlign(primer.c_str(),primer.length(),seq.c_str(), seq.length(), edlibConf);
  int OFFSET=5;//barcode_max_editd/2;
  int BC_START=result.startLocations[0]+ss.primer.length()-OFFSET;
  //result.endLocations[0]+1;//result.startLocations[0];//+primer.length()-OFFSET;
  string barcode_seq=seq.substr(BC_START,ss.temp_barcode.length()+2*OFFSET);
  
  //other search over the list of known barcodes.
  /**unordered_set<string>::iterator known_barcodes_itr=known_barcodes->begin();
  string barcodes_concat="";
  for(; known_barcodes_itr!=known_barcodes->end(); known_barcodes_itr++)
    barcodes_concat+=primer+(*known_barcodes_itr)+string(temp_barcode.length(),'X');
  edlibConf = edlibNewAlignConfig(flank_max_editd+barcode_max_editd+2*OFFSET,EDLIB_MODE_HW,EDLIB_TASK_LOC,NULL,0);
  result = edlibAlign(barcode_seq.c_str(), barcode_seq.length(),barcodes_concat.c_str(),barcodes_concat.length(),edlibConf);
  if(result.status == EDLIB_STATUS_OK && result.numLocations!=0){
    //    barcode.editd=result.editDistance-OFFSET;
    int idx=result.endLocations[0]/(temp_barcode.length()+temp_barcode.length()+primer.length());
    known_barcodes_itr=known_barcodes->begin();
    advance(known_barcodes_itr,idx);
    barcode.barcode=*(known_barcodes_itr);
       cout << barcode_seq << " " <<  barcode.barcode << " " << result.endLocations[0] << " " << result.startLocations[0]
       << " " << result.editDistance << " " << "!"<< endl;
    //realign to get edit distance
    result=edlibAlign(barcode.barcode.c_str(),barcode.barcode.length(),barcode_seq.c_str(), barcode_seq.length(),edlibConf);
    if(result.status == EDLIB_STATUS_OK && result.numLocations!=0){
      barcode.editd=result.editDistance; //add check..
       cout << barcode_seq << " " <<  barcode.barcode << " " << result.endLocations[0] << " " << result.startLocations[0]
       << " " << barcode.editd << " " << "!"<< endl;
    } else {
      cout << "No alignment: "<< barcode_seq << endl;
    }
    }**/
    
  edlibConf = edlibNewAlignConfig(barcode_max_editd,EDLIB_MODE_HW,EDLIB_TASK_LOC,NULL,0);
  unordered_set<string>::iterator known_barcodes_itr=known_barcodes->begin();
  for(; known_barcodes_itr!=known_barcodes->end(); known_barcodes_itr++){
    search_string=(*known_barcodes_itr);
    result = edlibAlign(search_string.c_str(), search_string.length(), barcode_seq.c_str(), barcode_seq.length(),edlibConf);
    if(result.status == EDLIB_STATUS_OK && result.numLocations!=0 && (result.editDistance <= barcode.editd)) {
      barcode.next_editd=barcode.editd;
      barcode.editd=result.editDistance;
      barcode.barcode=*known_barcodes_itr;
      barcode.umi=seq.substr(BC_START+result.endLocations[0]+1,ss.umi_seq.length());
      if(barcode.editd==0){ //perfect match, stop looking and return
	edlibFreeAlignResult(result);
	return(barcode); 
      }
    }
  }
  
  edlibFreeAlignResult(result);
  return(barcode);
}

vector<Barcode> big_barcode_search(string & sequence, unordered_set<string> & known_barcodes,
				   int max_flank_editd, int max_editd, SearchSeq & ss){
  vector<Barcode> return_vec;
  
  Barcode result=get_barcode(sequence,&known_barcodes,max_flank_editd,max_editd,ss);
  if(result.editd<=max_editd)
    return_vec.push_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if(return_vec.size()>0){
    string masked_sequence = sequence;
    for(int i=0; i<return_vec.size(); i++){
      int flank_length=return_vec.at(i).flank_end-return_vec.at(i).flank_start;
      masked_sequence.replace(return_vec.at(i).flank_start,flank_length,string(flank_length,'X'));
    }
    vector<Barcode> masked_res;
    masked_res=big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd,ss);
    return_vec.insert(return_vec.end(),masked_res.begin(),masked_res.end()); //add to result
  }
    
  return(return_vec);
}


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


// MAIN is here
int main(int argc, char **argv){

  //Variables to store user options
  //Set these to their defaults
  int expected_cells=0; //(d)
  int edit_distance=2; //(e)
  int flank_edit_distance=10; //(f)
  string out_read_filename=""; //(o)
  string out_stat_filename=""; //(t)
  bool split_file_by_barcode=false; //(s)
  bool remove_barcodes=true; //(r)
  string pattern=""; //(p)

  SearchSeq search_patterns;
  search_patterns.primer = "CTACACGACGCTCTTCCGATCT"; //(p)
  search_patterns.polyA = "TTTTTTTTT"; //(T)
  search_patterns.umi_seq = "????????????"; //(length u)
  search_patterns.temp_barcode = "????????????????"; //(length b)

  //Set of known barcodes 
  unordered_set<string> known_barcodes;

  /*** Pass command line option *****/
  int c;
  int params=1;
  ifstream file;
  string line;

  while((c =  getopt(argc, argv, "w:e:d:s:o:t:r:f:p:T:u:b:")) != EOF){
    switch(c){
    case 'w': { //w=white list of barcodes from short-reads
      string file_name(optarg);
      cerr << "Reading barcode file "<< file_name << endl;
      /**** READ BARCODE WHITE LIST FROM FILE ******/
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
      break; }
    case 'd':{
      expected_cells=atoi(optarg);
      params+=2;
      break;
    }
    case 's':{
      split_file_by_barcode=get_bool_opt_arg(optarg);
      params+=2;
      break;
    }
    case 'r':{
      remove_barcodes=get_bool_opt_arg(optarg);
      params+=2;
      break;
    }
      /**    case 't':{
      out_stat_filename=optarg;
      cerr << "Setting stats table file prefix to: " << out_stat_filename << endl;
      params+=2;
      break;
    }      
    case 'o':{
      out_read_filename=optarg;
      cerr << "Setting output read file prefix to: " << out_read_filename << endl;
      params+=2;
      break;
      }**/
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
      cerr << "Setting UMI length to " << bl << endl;
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
  
  //check that a read file is given
  if(params>=argc){
    cerr << "Reads input file missing" << endl;
    print_usage();
    exit(1);
  }
    
     
  // check that the reads file is okay
  string reads_file=argv[params];
  file.open(reads_file);
  if(!(file.good())){
    cerr << "Unable to open file " << reads_file << endl;
    print_usage();
    exit(1);
  }

  // more checks...
  if((expected_cells==0) & (known_barcodes.size()==0)){
    cerr << "Unclear how barcodes should be obtained. "
	 << "Either set -d > 0 or provide a barcode white list" << endl;
    print_usage();
    exit(1);
  }
  
  //set the output filenames
  out_stat_filename=reads_file+"-barcodes_assigned.txt"; 
  out_read_filename=reads_file+"-barcodes_assigned"; //add .fastq or .fasta later..
  
  // *** GET THE BARCODE FROM THE READ DATA DIRECTLY ****/
  if(expected_cells>0){
    cerr << "Getting barcodes directly from read data" << endl;
    // A few constanst related to capturing barcdoes from the data:
    int BEGIN_SEARCH_LENGTH=100; //length to search for quick check
    int MAX_NOVEL_BARCODE=expected_cells*10000; //Sample until max. 10x number expected
    int novel_bc_counter=0;
    unordered_map<string,int> novel_barcodes;
    //Now if option to get barcodes from data is set get these:
    while ( getline (file,line) && novel_bc_counter < MAX_NOVEL_BARCODE){
      //is the file a fastq or a fasta?
      istringstream line_stream(line);
      string read_id;
      line_stream >> read_id;
      bool is_fastq=true;
      if (read_id[0] == '>'){ is_fastq=false;
      } else if (read_id[0] == '@'){ //fasta
      } else {
	cerr << "Unknown read format... exiting" << endl; exit(1);
      }
      getline(file,line);
      string line_begin=line.substr(0,BEGIN_SEARCH_LENGTH);
      //get sequences with perfect adapter seq.
      Barcode barcode=get_barcode(line_begin,NULL,0,0,search_patterns);//fix
      if(barcode.flank_editd==0){
	novel_barcodes[barcode.barcode]++;
	novel_bc_counter++;
      }
      for(int s=0; is_fastq && s<2; s++) getline (file,line);
    }
    
    // now process the novel barcodes:
    // loop over the barcode list twice.
    // first to get the threshhold, then to apply it.
    unordered_map<string,int>::iterator novel_barcodes_itr=novel_barcodes.begin();
    vector<int> temp_value;
    
    for(; novel_barcodes_itr!=novel_barcodes.end(); novel_barcodes_itr++)
      temp_value.push_back(novel_barcodes_itr->second);
    sort(temp_value.begin(),temp_value.end());
    int temp_i; int temp_cut;
    if(temp_value.size() > 0){
      temp_i=temp_value.size()-expected_cells;
      temp_cut=temp_value.at(max(0,temp_i));
    }
    //loop again
    cerr << "Barcodes from white list: " << known_barcodes.size() << endl;
    novel_barcodes_itr=novel_barcodes.begin();
    for(; novel_barcodes_itr!=novel_barcodes.end(); novel_barcodes_itr++){
      if(novel_barcodes_itr->second >= temp_cut)
	known_barcodes.insert(novel_barcodes_itr->first);
    }
    cerr << "Total Barcodes (from white list and read): " << known_barcodes.size() << endl;
    cerr << "Done getting barcodes directly from read data" << endl;

    file.close(); //close and reopen file for reading
    file.open(reads_file);
  }

  //  ProfilerStart("my.prof");
  
  /********* FIND BARCODE IN READS (again, this time allowing mismatches *****/
  string sequence;
  int bc_count=0;
  int r_count=0;
  int multi_bc_count=0;

  ofstream out_stat_file;
  out_stat_file.open(out_stat_filename);
  out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tNextBestBarcodeEditDist"<<endl;

  cerr << "Assigning barcodes..." << endl;
  bool is_fastq=true;
  while ( getline (file,line) ){
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
    if(getline (file,line)){
      r_count++;
      //forward search
      vector<Barcode> vec_bc_for=big_barcode_search(line,known_barcodes,
						    flank_edit_distance,edit_distance,search_patterns);
      string rev_line=line;
      reverse_compliment(rev_line);
      //Check the reverse compliment of the read
      vector<Barcode> vec_bc_rev=big_barcode_search(rev_line,known_barcodes,
						    flank_edit_distance,edit_distance,search_patterns);
      if((vec_bc_for.size()+vec_bc_rev.size())>0)
	bc_count++;
      if((vec_bc_for.size()+vec_bc_rev.size())>1 ){
	multi_bc_count++;
	//	continue;
      }
      if(r_count % 100000 == 0)
	cerr << r_count/1000 << " thousand reads processed.." << endl;

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
  out_stat_file.close();
  
  cerr << "Number of reads processed: " << r_count << endl;
  cerr << "Number of reads assigned a barcode: " << bc_count << endl;
  cerr << "Number of reads assigned more than one barcode: " << multi_bc_count << endl;
  cerr << "All done!" << endl;

  //  ProfilerStop();
  
}
