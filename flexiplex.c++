// Copyright 2023 Nadia Davidson 
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
#include <map>
#include <stdlib.h>
#include <unistd.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <numeric>
#include <cstring>

#include "edlib.h"

using namespace std;

// Append .1 to version for dev code, remove for release
// e.g. 1.00.1 (dev) goes to 1.01 (release)
const static string VERSION="1.02.1"; 

struct PredefinedStruct {
  string description;
  string params;
};


// predefined settings for known barcode/search schemes
// new setting added to this map will automatically be available 
static const map< string, PredefinedStruct> predefinedMap = { 
  {"10x3v3", {"10x version 3 chemistry 3'", //option string and description for help information
	      "-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT -f 8 -e 2"}}, //settings
  {"10x3v2",{"10x version 2 chemistry 3'",
	     "-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTTTTTTT -f 8 -e 2"}},
  {"10x5v2",{"10x version 2 chemistry 5'",
	     "-x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ?????????? -x TTTCTTATATGGG -f 8 -e 2"}},
  {"grep",{"Simple grep-like search (edit distance up to 2)",
	   "-f 2 -k ? -b \'\' -u \'\' -i false"}}  
};

// the help information which is printed when a user puts in the wrong
// combination of command line options.
void print_usage(){
  cerr << "usage: flexiplex [options] [reads_input]\n\n";
  
  cerr << "  reads_input: a .fastq or .fasta file. Will read from stdin if empty.\n\n";

  cerr << "  options: \n";
  cerr << "     -k known_list   Either 1) a text file of expected barcodes in the first column,\n"; 
  cerr << "                     one row per barcode, or 2) a comma separate string of barcodes.\n";
  cerr << "                     Without this option, flexiplex will search and report possible barcodes.\n";
  cerr << "                     The generated list can be used for known_list in subsequent runs.\n"; 
  cerr << "     -i true/false   Replace read ID with barcodes+UMI, remove search strings\n";
  cerr << "                     including flanking sequenence and split read if multiple\n";
  cerr << "                     barcodes found (default: true).\n";
  cerr << "     -s true/false   Sort reads into separate files by barcode (default: false)\n";
  cerr << "     -n prefix   Prefix for output filenames.\n";
  cerr << "     -e N   Maximum edit distance to barcode (default 2).\n";
  cerr << "     -f N   Maximum edit distance to primer+polyT (default 8).\n";
  cerr << "     -p N   Number of threads (default: 1).\n\n";
  
  cerr << "  Specifying adaptor / barcode structure : \n";
  cerr << "     -x sequence Append flanking sequence to search for\n";
  cerr << "     -b sequence Append the barcode pattern to search for\n";
  cerr << "     -u sequence Append the UMI pattern to search for\n";
  cerr << "     Notes:\n";
  cerr << "          The order of these options matters\n";
  cerr << "          ? - can be used as a wildcard\n";
  cerr << "     When no search pattern x,b,u option is provided, the following default pattern is used: \n";
  cerr << "          primer: CTACACGACGCTCTTCCGATCT\n";
  cerr << "          barcode: ????????????????\n";
  cerr << "          UMI: ????????????\n";
  cerr << "          polyT: TTTTTTTTT\n";
  cerr << "     which is the same as providing: \n";
  cerr << "         -x CTACACGACGCTCTTCCGATCT -b ???????????????? -u ???????????? -x TTTTTTTTT\n\n";
  
  cerr << "  Predefined search schemes:\n";
  auto pds_itr = predefinedMap.begin();
  for(; pds_itr!=predefinedMap.end(); pds_itr++){
    cerr << "    -d " << pds_itr->first << "\t\t" << pds_itr->second.description << ", equivalent to:\n";
    cerr << "\t\t\t\t"<< pds_itr->second.params << "\n";
  }
  cerr << "\n     -h     Print this usage information.\n\n";
  cerr << "Have a different barcode scheme you would like Flexiplex to work with? Post a request at:\n" ;
  cerr << "https://github.com/DavidsonGroup/flexiplex/issues\n" ;
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

//Holds the found barcode and associated information 
struct Barcode {
  string barcode;
  string umi;
  int editd;
  int flank_editd;
  int flank_start;
  int flank_end;
  bool unambiguous;
} ;

struct SearchResult {
  string read_id;
  string qual_scores;
  string line;
  string rev_line;
  vector<Barcode> vec_bc_for;
  vector<Barcode> vec_bc_rev;
};

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){
  
  std::size_t len1 = s1.size()+1, len2 = s2.size()+1;
  const char * s1_c = s1.c_str(); const char * s2_c = s2.c_str();

  vector< unsigned int> dist_holder(len1*len2);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0]=0; //[0][0]
  for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
  for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];
    
  int best=len2;
  end=len1-1;

  //loop over the distance matrix elements and calculate running distance
  for (unsigned int j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (unsigned int i = 1; i < len1; ++i) {
      int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1; // are the bases the same?

      // if yes, no need to increment distance
      if (sub == 0) {
        dist_holder[i * len2 + j] = dist_holder[(i - 1) * len2 + (j - 1)];
      }

      // otherwise add insertion, deletion or substitution
      else {
        // clang-format off
        dist_holder[i*len2+j] = std::min({ //j+i*len2  //set[i][j]
          dist_holder[(i-1)*len2+j] + 1, //[i-1][j]
          dist_holder[i*len2+(j-1)] + 1, //[i][j-1]
          dist_holder[(i-1)*len2+(j-1)] + 1}); // ((s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1) });
        // clang-format on
      }
      if (dist_holder[i * len2 + j] <= max_editd) {
        any_below_threshold = true;
      }

      // if this is the last row in j
      if (j == (len2 - 1) && dist_holder[i * len2 + j] < best) {
        // check if this is the best running score
        best = dist_holder[i * len2 + j];
        end = i; // update the end position of alignment
      }
    }

    // early exit to save time.
    if(!any_below_threshold) {
      return(100);
    }
  }

  return best; // return edit distance
}


// Given a string 'seq' search for substring with primer and polyT sequence followed by
// a targeted search in the region for barcode
// Sequence seearch is performed using edlib

Barcode get_barcode(string & seq,
		    unordered_set<string> *known_barcodes,
		    int flank_max_editd,
		    int barcode_max_editd,
        const std::vector<std::pair<std::string, std::string>> &search_pattern) {
  
  const int OFFSET=5; //wiggle room in bases of the expected barcode start site to search.

  //initialise struct variables for return:
  Barcode barcode;
  barcode.editd=100; barcode.flank_editd=100; barcode.unambiguous = false;

  //initialise edlib configuration
  // Use IUPAC codes
  EdlibEqualityPair additionalEqualities[28] = {
    {'R', 'A'}, {'R', 'G'},
    {'K', 'G'}, {'K', 'T'},
    {'S', 'G'}, {'S', 'C'},
    {'Y', 'C'}, {'Y', 'T'},
    {'M', 'A'}, {'M', 'C'},
    {'W', 'A'}, {'W', 'T'},
    {'B', 'C'}, {'B', 'G'}, {'B', 'T'},
    {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
    {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
    {'D', 'A'}, {'D', 'G'}, {'D', 'T'},
    {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
  };
  EdlibAlignConfig edlibConf = {flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 28};

  // concatenate patterns in search_pattern in insertion order
  std::string search_string;
  for (const auto &pair : search_pattern) {
    search_string += pair.second;
  }

  //search for the concatenated pattern
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK | result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode); // no match found - return
  } //fill in info about found primer and polyT location
  barcode.flank_editd=result.editDistance;
  barcode.flank_start=result.startLocations[0];
  barcode.flank_end=result.endLocations[0];

  // Extract sub-patterns from aligment directly
  std::vector<long unsigned int> subpattern_lengths;
  for (const auto &pair : search_pattern) {
    subpattern_lengths.push_back(pair.second.length());
  }

  std::vector<long unsigned int> subpattern_ends;
  subpattern_ends.resize(subpattern_lengths.size());
  std::partial_sum(subpattern_lengths.begin(), subpattern_lengths.end(), subpattern_ends.begin());

  vector<int> read_to_subpatterns;
  read_to_subpatterns.reserve(subpattern_ends.size() + 1);
  read_to_subpatterns.emplace_back(barcode.flank_start);

  // initialise pointers
  int i_read = barcode.flank_start;
  int i_pattern = 0;
  int i_subpattern = 0;

  // walk through edlib aligment
  // 0 for match
  // 1 for insertion to target
  // 2 for insertion to query
  // 3 for mismatch 
  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != 1) {
      i_read++;
    }
    if (value != 2) {
      i_pattern++;
    }
    if (i_pattern >= subpattern_ends[i_subpattern]) {
      read_to_subpatterns.emplace_back(i_read);
      i_subpattern++;
    }
  }

  edlibFreeAlignResult(result);
  

  // Work out the index of BC and UMI in the pattern
  // TODO: Should this be done in main to be more efficient?
  // TODO: Handle edge cases where BC / UMI is not speficied in the pattern
  int bc_index = -1, umi_index = -1;
  auto it_pattern =
      std::find_if(search_pattern.begin(), search_pattern.end(),
                   [](const std::pair<std::string, std::string> &pair) {
                     return pair.first == "BC";
                   });
  if (it_pattern != search_pattern.end()) {
    bc_index = std::distance(search_pattern.begin(), it_pattern);
  } else {
    // error
  }
  it_pattern =
      std::find_if(search_pattern.begin(), search_pattern.end(),
                   [](const std::pair<std::string, std::string> &pair) {
                     return pair.first == "UMI";
                   });
  if (it_pattern != search_pattern.end()) {
    umi_index = std::distance(search_pattern.begin(), it_pattern);
  } else {
    // error
  }

  //if not checking against known list of barcodes, return sequence after the primer
  //also check for a perfect match straight up as this will save computer later.
  // 
  // read_to_subpatterns[subpattern_index] gives start of subpattern in the read
  std::string exact_bc = seq.substr(
      read_to_subpatterns[bc_index],
      search_pattern[bc_index].second.length());
  if (known_barcodes->size()==0 || (known_barcodes->find(exact_bc) != known_barcodes->end())){ 
    barcode.barcode=exact_bc;
    barcode.editd=0;
    barcode.unambiguous=true;
    if (umi_index == -1) {
      barcode.umi = "";
    } else {
      barcode.umi =
          seq.substr(read_to_subpatterns[umi_index],
                     search_pattern[umi_index].second.length());
    }
    return(barcode);
  }
  
  // otherwise widen our search space and the look for matches with errors

  int left_bound = max(
    read_to_subpatterns[bc_index] - OFFSET, // widen the search by using an OFFSET
    0                                       // set a maximum starting character index of 0
  );
  int max_length = search_pattern[bc_index].second.length() + 2 * OFFSET;

  std::string barcode_seq = seq.substr(left_bound, max_length);
 
  //iterate over all the known barcodes, checking each sequentially
  unordered_set<string>::iterator known_barcodes_itr=known_barcodes->begin();
  unsigned int editDistance, endDistance;

  for(; known_barcodes_itr != known_barcodes->end(); known_barcodes_itr++){
    search_string = *known_barcodes_itr; //known barcode to check against
    editDistance = edit_distance(barcode_seq, search_string, endDistance, barcode_max_editd);

    if (editDistance == barcode.editd) {
      barcode.unambiguous = false;    
    } else if (editDistance < barcode.editd && editDistance <= barcode_max_editd) { // if best so far, update
      barcode.unambiguous = true;
      barcode.editd = editDistance;
      barcode.barcode = *known_barcodes_itr;

      if (umi_index == -1) {
        barcode.umi = "";
      } else if (umi_index == bc_index + 1) {
        // left_bound: start of barcode_seq
        // + endDistance: end of barcode
        // i.e. start of UMI
        barcode.umi = seq.substr(
          left_bound + endDistance,
          search_pattern[umi_index].second.length()
        ); // assumes no error in UMI seq.
      } else if (umi_index == bc_index - 1) {
        // Use the start of BC according to edit_distance(barcode_seq, ...) and go backwords
        // read_to_subpatterns[bc_index] - OFFSET + endDistance - search_pattern[bc_index].second.length():
        // = start of BC
        barcode.umi = seq.substr(
          left_bound + endDistance 
                     - search_pattern[bc_index].second.length() 
                     - search_pattern[umi_index].second.length(),
          search_pattern[umi_index].second.length()
        );
      } else {
        // BC and UMI not next to eachother, grab UMI according to aligment
        barcode.umi = seq.substr(
          read_to_subpatterns[umi_index],
          search_pattern[umi_index].second.length()
        ); // assumes no error in UMI seq.
      }
      
      //if perfect match is found we're done.
      if (editDistance == 0) {
      	return(barcode);
      }
    }

  }
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
vector<Barcode> big_barcode_search(string & sequence, unordered_set<string> & known_barcodes, int max_flank_editd, int max_editd, const std::vector<std::pair<std::string, std::string>> &search_pattern) {
  vector<Barcode> return_vec; //vector of all the barcodes found

  //search for barcode
  Barcode result=get_barcode(sequence,&known_barcodes,max_flank_editd,max_editd, search_pattern); //,ss);
  if(result.editd<=max_editd && result.unambiguous) //add to return vector if edit distance small enough
    return_vec.push_back(result);
  
  //if a result was found, mask out the flanking sequence and search again in case there are more.
  if(return_vec.size()>0){
    string masked_sequence = sequence;
    for(int i=0; i<return_vec.size(); i++){
      int flank_length=return_vec.at(i).flank_end-return_vec.at(i).flank_start;
      masked_sequence.replace(return_vec.at(i).flank_start,flank_length,string(flank_length,'X'));
    } //recursively call this function until no more barcodes are found
    vector<Barcode> masked_res;
    masked_res=big_barcode_search(masked_sequence,known_barcodes,max_flank_editd,max_editd, search_pattern); //,ss);
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
    cerr << "Unknown argument to boolean option\n";
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
	       << vec_bc.at(b).umi << "\n";
  }
}

void print_line(string id, string read, string quals, ostream & out_stream){

  //flag for read format
  bool is_fastq=!(quals==""); //no quality scores passed = fasta

  //output to the new read file
    if(is_fastq)
      out_stream << "@" << id << "\n";
    else
      out_stream << ">" << id << "\n";
    out_stream << read << "\n";
    if(is_fastq){
      out_stream << "+" << id << "\n";
      out_stream << quals << "\n";
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

// separated out from main so that this can be run with threads
void search_read(vector<SearchResult> & reads, unordered_set<string> & known_barcodes,
			 int flank_edit_distance, int edit_distance,
       const std::vector<std::pair<std::string, std::string>> &search_pattern) {
  
  for(int r=0; r<reads.size(); r++){
    
    //forward search
    reads[r].vec_bc_for=big_barcode_search(reads[r].line,
					   known_barcodes,
					   flank_edit_distance,
					   edit_distance,
             search_pattern);
    reads[r].rev_line=reads[r].line;
    reverse_compliment(reads[r].rev_line);
    //Check the reverse compliment of the read
    reads[r].vec_bc_rev=big_barcode_search(reads[r].rev_line,
					   known_barcodes,
					   flank_edit_distance,
					   edit_distance,
             search_pattern);
  }
}

// MAIN is here!!
int main(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);

  cerr << "FLEXIPLEX " << VERSION << "\n";

  // Variables to store user options
  // Set these to their defaults
  int expected_cells = 0;      //(d)
  int edit_distance = 2;       //(e)
  int flank_edit_distance = 8; //(f)

  // set the output filenames
  string out_stat_filename = "reads_barcodes.txt";
  string out_bc_filename = "barcodes_counts.txt";
  string out_filename_prefix = "flexiplex"; //(n)

  bool split_file_by_barcode = false; //(s)
  bool remove_barcodes = true;        //(r)

  std::vector<std::pair<std::string, std::string>> search_pattern;

  // Set of known barcodes
  unordered_set<string> known_barcodes;
  unordered_set<string> found_barcodes;

  // threads
  int n_threads = 1;

  /*** Pass command line option *****/
  int c;
  int params = 1;
  ifstream file;
  string line;

  vector<char *> myArgs(argv, argv + argc);

  while ((c = getopt(myArgs.size(), myArgs.data(),
                     "d:k:i:b:u:x:e:f:n:s:hp:")) != EOF) {
    switch (c) {
    case 'd': { // d=predefined list of settings for various search/barcode
                // schemes
      if (predefinedMap.find(optarg) == predefinedMap.end()) {
        cerr << "Unknown argument, " << optarg
             << ", passed to option -d. See list below for options.\n";
        print_usage();
        exit(1);
      } // else
      cerr << "Using predefined settings for " << optarg << ".\n";
      istringstream settingsStream(predefinedMap.find(optarg)->second.params);
      string token;
      vector<char *> newArgv;
      while (settingsStream >> token) { // code created with the help of chatGPT
        token.erase(remove(token.begin(), token.end(), '\''), token.end());
        char *newArg =
            new char[token.size() +
                     1]; // +1 for null terminator //need to delete these later.
        strcpy(newArg, token.c_str());
        newArgv.push_back(newArg); // Append the token to the new argv
      }
      params += 2;
      myArgs.insert(myArgs.begin() += params, newArgv.begin(), newArgv.end());
      break;
    }
    case 'k': { // k=list of known barcodes
      string file_name(optarg);
      string bc;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      cerr << "Setting known barcodes from " << file_name << "\n";
      if (!(file.good())) { // if the string given isn't a file
        stringstream bc_list(file_name);
        string s;
        while (getline(bc_list, bc, ',')) // tokenize
          known_barcodes.insert(bc);
      } else {
        // otherwise get the barcodes from the file..
        while (getline(file, line)) {
          istringstream line_stream(line);
          line_stream >> bc;
          known_barcodes.insert(bc);
        }
        file.close();
      }
      cerr << "Number of known barcodes: " << known_barcodes.size() << "\n";
      if (known_barcodes.size() == 0) {
        print_usage();
        exit(1); // case barcode file is empty
      }
      params += 2;
      break;
    }
    case 'i': {
      remove_barcodes = get_bool_opt_arg(optarg);
      cerr << "Setting read IDs to be replaced: " << remove_barcodes << "\n";
      params += 2;
      break;
    }
    case 'e': {
      edit_distance = atoi(optarg);
      cerr << "Setting max barcode edit distance to " << edit_distance << "\n";
      params += 2;
      break;
    }
    case 'f': {
      flank_edit_distance = atoi(optarg);
      cerr << "Setting max flanking sequence edit distance to "
           << flank_edit_distance << "\n";
      params += 2;
      break;
    }
    // x, u, b arguments inserts subpatterns to search_pattern
    case 'x': {
      search_pattern.push_back(std::make_pair("Unnamed Seq", optarg));
      cerr << "Adding flank sequence to search for: " << optarg << "\n";
      params += 2;
      break;
    }
    case 'u': {
      search_pattern.push_back(std::make_pair("UMI", optarg));
      cerr << "Setting UMI to search for: " << optarg << "\n";
      params += 2;
      break;
    }
    case 'b': {
      search_pattern.push_back(std::make_pair("BC", optarg));
      cerr << "Setting barcode to search for: " << optarg << "\n";
      params += 2;
      break;
    }
    case 'h': {
      print_usage();
      exit(0);
    }
    case 'n': {
      out_filename_prefix = optarg;
      cerr << "Setting output filename prefix to: " << out_filename_prefix
           << "\n";
      params += 2;
      break;
    }
    case 's': {
      split_file_by_barcode = get_bool_opt_arg(optarg);
      cerr << "Split read output into separate files by barcode: "
           << split_file_by_barcode << "\n";
      int max_split_bc = 50;
      if (known_barcodes.size() > max_split_bc) {
        cerr << "Too many barcodes to split into separate files: "
             << known_barcodes.size() << "> " << max_split_bc << "\n";
        split_file_by_barcode = false;
      }
      params += 2;
      break;
    }
    case 'p': {
      n_threads = atoi(optarg);
      cerr << "Setting number of threads to " << n_threads << endl;
      params += 2;
      break;
    }
    case '?': // other unknown options
      cerr << "Unknown option.. stopping" << endl;
      print_usage();
      exit(1);
    }
  }

  // default case when no x, u, b is speficied
  if (search_pattern.empty()) {
    cerr << "Using default search pattern: " << endl;
    search_pattern = {{"primer", "CTACACGACGCTCTTCCGATCT"},
                      {"BC", std::string(16, '?')},
                      {"UMI", std::string(12, '?')},
                      {"polyA", std::string(9, 'T')}};
    for (auto i : search_pattern)
      std::cerr << i.first << ": " << i.second << "\n";
  }

  cerr << "For usage information type: flexiplex -h" << endl;

  istream *in;
  ifstream reads_ifs;

  // check that a read file is given
  if (params >= myArgs.size()) {
    cerr << "No filename given... getting reads from stdin..." << endl;
    in = &std::cin;
  } else {
    // check that the reads fileis okay
    string reads_file = myArgs[params];
    reads_ifs.open(reads_file);
    if (!(reads_ifs.good())) {
      cerr << "Unable to open file " << reads_file << endl;
      print_usage();
      exit(1);
    }
    in = &reads_ifs;
  }

  /********* FIND BARCODE IN READS ********/
  string sequence;
  int bc_count = 0;
  int r_count = 0;
  int multi_bc_count = 0;

  ofstream out_stat_file;
  out_stat_filename = out_filename_prefix + "_" + out_stat_filename;
  out_bc_filename = out_filename_prefix + "_" + out_bc_filename;
  params += 2;

  if (known_barcodes.size() > 0) {
    out_stat_file.open(out_stat_filename);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI\n";
  }

  cerr << "Searching for barcodes...\n";

  bool is_fastq = true;
  unordered_map<string, int> barcode_counts;
  string read_id_line;

  if (getline(*in, read_id_line)) { // check the first line for file type
    if (read_id_line[0] == '>') {
      is_fastq = false;
    } else if (read_id_line[0] == '@') { // fasta - ignore!
    } else {
      cerr << "Unknown read format... exiting" << endl;
      exit(1);
    }
  }

  while (getline(*in, line)) {
    const int buffer_size = 2000; // number of reads to pass to one thread.

    vector<vector<SearchResult>> sr_v(n_threads);
    for (int i = 0; i < n_threads; i++) {
      sr_v[i] = vector<SearchResult>(buffer_size);
    }

    vector<thread> threads(n_threads);

    // get n_threads*buffer number of reads..
    for (int t = 0; t < n_threads; t++) {
      for (int b = 0; b < buffer_size; b++) {
        SearchResult &sr = sr_v[t][b];

        sr.line = line;
        string read_id;

        istringstream line_stream(read_id_line);
        line_stream >> sr.read_id;
        sr.read_id.erase(0, 1);

        if (!is_fastq) { // fasta (account for multi-lines per read)
          string buffer_string;
          while (getline(*in, buffer_string) && buffer_string[0] != '>')
            sr.line += buffer_string;
          read_id_line = buffer_string;
        } else { // fastq (get quality scores)
          for (int s = 0; s < 2; s++) {
            getline(*in, sr.qual_scores);
          }
          getline(*in, read_id_line);
        }

        r_count++; // progress counter
        if (// display for every 10,000 reads first
            (r_count < 100000 && r_count % 10000 == 0) ||
            // and then every 100,000 reads after
            (r_count % 100000 == 0)) {
          cerr << r_count / ((double)1000000) << " million reads processed.."
               << endl;
        }

        // this is quite ugly, must be a better way to do this..
        if (b == buffer_size - 1 && t == n_threads - 1) {
          break; // if it's the last in the chunk don't getline as this happens
                 // in the while statement
        } else if (!getline(*in, line)) { // case we are at the end of the
                                          // reads.
          sr_v[t].resize(b + 1);
          threads[t] = std::thread(search_read, ref(sr_v[t]),
                                   ref(known_barcodes), flank_edit_distance,
                                   edit_distance, ref(search_pattern));
          for (int t2 = t + 1; t2 < n_threads; t2++) {
            sr_v[t2].resize(0);
          }

          goto print_result; // advance the line
        }
      }
      // send reads to the thread
      threads[t] =
          std::thread(search_read, ref(sr_v[t]), ref(known_barcodes),
                      flank_edit_distance, edit_distance, ref(search_pattern));
    }

  print_result:
    // START print_result LABEL...

    // loop over the threads and print out ther results
    for (int t = 0; t < sr_v.size(); t++) {
      if (sr_v[t].size() > 0)
        threads[t].join(); // wait for the threads to finish before printing

      for (int r = 0; r < sr_v[t].size(); r++) { // loop over the reads

        for (int b = 0; b < sr_v[t][r].vec_bc_for.size(); b++)
          barcode_counts[sr_v[t][r].vec_bc_for.at(b).barcode]++;
        for (int b = 0; b < sr_v[t][r].vec_bc_rev.size(); b++)
          barcode_counts[sr_v[t][r].vec_bc_rev.at(b).barcode]++;

        if ((sr_v[t][r].vec_bc_for.size() + sr_v[t][r].vec_bc_rev.size()) > 0)
          bc_count++;
        if ((sr_v[t][r].vec_bc_for.size() + sr_v[t][r].vec_bc_rev.size()) > 1) {
          multi_bc_count++;
        }

        if (known_barcodes.size() !=
            0) { // if we are just looking for all possible barcodes don't
                 // output reads etc.

          print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_for, out_stat_file);
          print_stats(sr_v[t][r].read_id, sr_v[t][r].vec_bc_rev, out_stat_file);

          print_read(sr_v[t][r].read_id + "_+", sr_v[t][r].line,
                     sr_v[t][r].qual_scores, sr_v[t][r].vec_bc_for,
                     out_filename_prefix, split_file_by_barcode, found_barcodes,
                     remove_barcodes);
          reverse(sr_v[t][r].qual_scores.begin(), sr_v[t][r].qual_scores.end());
          if (remove_barcodes || sr_v[t][r].vec_bc_for.size() ==
                                     0) // case we just want to print read once
                                        // if multiple bc found.
            print_read(sr_v[t][r].read_id + "_-", sr_v[t][r].rev_line,
                       sr_v[t][r].qual_scores, sr_v[t][r].vec_bc_rev,
                       out_filename_prefix, split_file_by_barcode,
                       found_barcodes, remove_barcodes);
        }
      }
    }

    // END print_result LABEL
  }

  reads_ifs.close();

  // print summary statistics
  cerr << "Number of reads processed: " << r_count << "\n";
  cerr << "Number of reads where a barcode was found: " << bc_count << "\n";
  cerr << "Number of reads where more than one barcode was found: "
       << multi_bc_count << "\n";
  cerr << "All done!" << endl;

  if (known_barcodes.size() > 0) {
    out_stat_file.close();
    return (0);
  }

  if (barcode_counts.size() == 0)
    return (0);

  typedef std::pair<std::string, int> pair;
  vector<pair> bc_vec;

  copy(barcode_counts.begin(), barcode_counts.end(),
       back_inserter<vector<pair>>(bc_vec));

  sort(bc_vec.begin(), bc_vec.end(), [](const pair &l, const pair &r) {
    if (l.second != r.second)
      return l.second > r.second;
    return l.first < r.first;
  });

  vector<int> hist(bc_vec[0].second);
  ofstream out_bc_file;

  out_bc_file.open(out_bc_filename);

  for (auto const &bc_pair : bc_vec) {
    out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
    hist[bc_pair.second - 1]++;
  }

  out_bc_file.close();

  cout << "Reads\tBarcodes"
       << "\n";
  for (int i = hist.size() - 1; i >= 0; i--)
    cout << i + 1 << "\t" << hist[i] << "\n";
}
