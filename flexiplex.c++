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
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <queue>
#include <memory>

#include "edlib.h"
#include <cstring>

using namespace std;

// Append .1 to version for dev code, remove for release
// e.g. 1.00.1 (dev) goes to 1.01 (release)
const static string VERSION="1.02.5";

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
	   "-f 2 -k ? -b '' -u '' -i false"}}
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
  cerr << "     -c true/false   Add a _C suffix to the read identifier of any chimeric reads\n";
  cerr << "                     (default: false). For instance if,\n";
  cerr << "                       @BC_UMI#READID_+1of2\n";
  cerr << "                     is chimeric, it will become:\n";
  cerr << "                       @BC_UMI#READID_+1of2_C\n";
  cerr << "     -n prefix       Prefix for output filenames.\n";
  cerr << "     -e N            Maximum edit distance to barcode (default 2).\n";
  cerr << "     -f N            Maximum edit distance to primer+polyT (default 8).\n";
  cerr << "     -p N            Number of threads (default: 1).\n\n";

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
  cerr << "https://github.com/DavidsonGroup/flexiplex/issues\n\n" ;

  cerr << "If you use Flexiplex in your research, please cite our paper:\n" ;
  cerr << "O. Cheng et al., Flexiplex: a versatile demultiplexer and search tool for omics data, Bioinformatics, Volume 40, Issue 3, 2024 \n";

  cerr << endl;
}

// complement nucleotides - used to reverse complement string
char complement(char& c){
  switch(c){
  case 'A' : return 'T';
  case 'T' : return 'A';
  case 'G' : return 'C';
  case 'C' : return 'G';
  default: return 'N';
  }
}

//Inplace reverse complement
void reverse_complement(string & seq){
   reverse(seq.begin(),seq.end());
   transform(seq.begin(),seq.end(),seq.begin(),complement);
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
  int count;
  bool chimeric;
};

// A thread-safe queue for our producer-consumer model
template <typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cond_;

public:
    void push(T value) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(std::move(value));
        cond_.notify_one();
    }

    // Pop an element from the queue. Returns true on success, false if the queue is empty and likely to remain so.
    bool pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        // Wait until the queue is not empty.
        cond_.wait(lock, [this] { return !queue_.empty(); });
        
        // Although we waited, it's possible to be woken up spuriously.
        // Or, another thread could have popped the element.
        if (queue_.empty()) {
            return false;
        }
        
        value = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    void notify_all_finish() {
        std::lock_guard<std::mutex> lock(mutex_);
        cond_.notify_all();
    }
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

// extract UMI from the read after barcode matching
std::string get_umi(const std::string &seq,
                    const std::vector<std::pair<std::string, std::string>> &search_pattern,
                    const std::vector<int> &read_to_subpatterns,
                    const int umi_index, const int bc_index,
                    const bool sliding_window_match, // if true, use left_bound and endDistance
                    const int left_bound,
                    const int endDistance
                    ) {

  int umi_start, umi_length;
  std::string umi_pad = "";
  umi_length = search_pattern[umi_index].second.length();

  if (umi_index == -1) {
    return ""; // protocol does not have UMI

  } else if (umi_index == bc_index + 1) {
    // UMI right after BC
    if (sliding_window_match) {
      umi_start = left_bound + endDistance;
    } else {
      umi_start = read_to_subpatterns[bc_index] + search_pattern[bc_index].second.length();
    }
    if (seq.length() < umi_start + umi_length) {
      // read not long enough, pad N to the end
      umi_length = seq.length() - umi_start;
      umi_pad = std::string(search_pattern[umi_index].second.length() - umi_length, 'N');
    }
    return seq.substr(umi_start, umi_length) + umi_pad;

  } else if (umi_index == bc_index - 1) {
    // UMI right before BC
    int bc_start = sliding_window_match ? left_bound + endDistance : read_to_subpatterns[bc_index];
    // umi should start umi_offset bases before BC
    int umi_offset = search_pattern[bc_index].second.length() + search_pattern[umi_index].second.length();
    if (bc_start < umi_offset) {
      // not enough bases before BC
      umi_pad = std::string(umi_offset - bc_start, 'N');
      umi_start = 0;
      umi_length -= umi_offset - bc_start;
    } else {
      umi_start = bc_start - umi_offset;
    }
    return umi_pad + seq.substr(umi_start, umi_length);

  } else {
    // UMI not next to BC, no idea which side was truncated
    if (read_to_subpatterns.size() > umi_index + 1) {
      // UMI is not the last subpattern
      // use the start of the next subpattern as the end of UMI
      umi_start = read_to_subpatterns[umi_index];
      umi_length = min((int)read_to_subpatterns[umi_index + 1] - umi_start, umi_length);
      umi_pad = std::string(search_pattern[umi_index].second.length() - umi_length, 'N');
      return seq.substr(umi_start, umi_length) + umi_pad;
    } else {
      // UMI is the last subpattern
      umi_start = read_to_subpatterns[umi_index];
      umi_length = min((int) seq.length() - umi_start, umi_length);
      umi_pad = std::string(search_pattern[umi_index].second.length() - umi_length, 'N');
      return seq.substr(umi_start, umi_length) + umi_pad;
    }
  }
}

// Given a string 'seq' search for substring with primer and polyT sequence followed by
// a targeted search in the region for barcode
// Sequence seearch is performed using edlib

Barcode get_barcode(string & seq,
		    const unordered_set<string> *known_barcodes,
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

  // shorter than the pattern, skip search
  if (seq.length() < search_string.length()) {
    return barcode;
  }

  //search for the concatenated pattern
  EdlibAlignResult result = edlibAlign(search_string.c_str(), search_string.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK || result.numLocations==0 ){
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
    if (i_subpattern < subpattern_ends.size() && i_pattern >= subpattern_ends[i_subpattern]) {
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
    barcode.umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index, bc_index, false, 0, 0);
    return(barcode);
  }

  // otherwise widen our search space and the look for matches with errors

  int left_bound = max(
    (int)read_to_subpatterns[bc_index] - OFFSET, // widen the search by using an OFFSET
    0                                       // set a maximum starting character index of 0
  );
  int max_length = search_pattern[bc_index].second.length() + 2 * OFFSET;
  if (left_bound + max_length > seq.length()) {
      max_length = seq.length() - left_bound;
  }

  std::string barcode_seq = seq.substr(left_bound, max_length);

  //iterate over all the known barcodes, checking each sequentially
  auto known_barcodes_itr=known_barcodes->begin();
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
      barcode.umi = get_umi(seq, search_pattern, read_to_subpatterns, umi_index, bc_index, true, left_bound, endDistance);

      //if perfect match is found we're done.
      if (editDistance == 0) {
      	return(barcode);
      }
    }

  }
  return(barcode); //return the best matched barcode and associated information
}

//search a read for one or more barcodes (parent function that calls get_barcode)
vector<Barcode> big_barcode_search(const string & sequence, const unordered_set<string> & known_barcodes, int max_flank_editd, int max_editd, const std::vector<std::pair<std::string, std::string>> &search_pattern) {
    vector<Barcode> return_vec; //vector of all the barcodes found
    string masked_sequence = sequence; // Work on a copy

    while (true) {
      //search for barcode
      Barcode result=get_barcode(masked_sequence, &known_barcodes, max_flank_editd, max_editd, search_pattern);
      if(result.editd <= max_editd && result.unambiguous) { //add to return vector if edit distance small enough
        return_vec.push_back(result);
        // Mask the found region to prevent re-finding it
        int flank_length = result.flank_end - result.flank_start;
        if (flank_length > 0)
          masked_sequence.replace(result.flank_start, flank_length, string(flank_length, 'X'));
        else // if flank_length is 0, we risk an infinite loop
            break;
      } else {
        break; // No more barcodes found
      }
    }
    return(return_vec);
}


// utility function to check true/false input options
bool get_bool_opt_arg(string value){
  transform(value.begin(), value.end(), value.begin(), ::tolower);
  if (value.compare("true")==0 || value.compare("t")==0 || value.compare("1")==0){
    return true;
  } else if (value.compare("false")==0 || value.compare("f")==0 || value.compare("0")==0){
    return false;
  } else {
    cerr << "Unknown argument to boolean option\n";
    print_usage();
    exit(1);
  }
}

// print information about barcodes
void print_stats(const string& read_id, const vector<Barcode> & vec_bc, ostream & out_stream){
  for(int b=0; b<vec_bc.size() ; b++){
    out_stream << read_id << "\t"
	       << vec_bc.at(b).barcode << "\t"
	       << vec_bc.at(b).flank_editd << "\t"
	       << vec_bc.at(b).editd << "\t"
	       << vec_bc.at(b).umi << "\n";
  }
}

void print_line(const string& id, const string& read, const string& quals, ostream & out_stream){

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
void print_read(const string& read_id, const string& read, string qual,
		const vector<Barcode> & vec_bc, const string& prefix,
		bool split, map<string, ofstream> & barcode_file_streams,
        mutex & cout_mutex, mutex & file_mutex,
		bool trim_barcodes,
    bool chimeric){

    auto vec_size = vec_bc.size();

    //loop over the barcodes found... usually will just be one
    for (int b = 0; b < vec_size; b++) {

      // format the new read id. Using FLAMES format.
      stringstream ss;
      ss << (b + 1) << "of" << vec_size;
      if (chimeric) {
        ss << "_" << "C";
      }

      string barcode = vec_bc.at(b).barcode;
      // also add the proper FASTQ way: \tCB:Z:barcode\tUB:Z:umi
      string new_read_id =
          barcode + "_" + vec_bc.at(b).umi + "#" + read_id + ss.str() + "\tCB:Z:" + barcode + "\tUB:Z:" + vec_bc.at(b).umi;

      // work out the start and end base in case multiple barcodes
      // note: read_start+1, see issue #63
      int read_start = vec_bc.at(b).flank_end + 1;
      int read_length = read.length() - read_start;

      for (int f = 0; f < vec_size; f++) {
        int temp_read_length = vec_bc.at(f).flank_start - read_start;
        if (temp_read_length >= 0 && temp_read_length < read_length)
          read_length = temp_read_length;
      }
      string qual_new = ""; // don't trim the quality scores if it's a fasta file

      if (qual != "") {
        if((read_start+read_length)>(qual.length())){
            lock_guard<mutex> lock(cout_mutex);
            cerr << "WARNING: sequence and quality lengths diff for read: " << read_id << ". Ignoring read." << endl;
            return;
        }
        qual_new = qual.substr(read_start, read_length);
      }
      string read_new = read.substr(read_start, read_length);

      if (b == 0 && !trim_barcodes) { // override if read shouldn't be cut
        new_read_id = read_id;
        read_new = read;
        qual_new = qual;
        b = vec_size; // force loop to exit after this iteration
      }

      // Skip reads that become empty after trimming
      if (read_new.length() == 0) {
        continue;
      }

      if (split) { // to a file if spliting by barcode
        lock_guard<mutex> lock(file_mutex);
        // Check if the file stream is already open
        if (barcode_file_streams.find(barcode) == barcode_file_streams.end()) {
            // If not, open it for the first time and add it to the map
            string outname = prefix + "_" + barcode + (qual.empty() ? ".fasta" : ".fastq");
            barcode_file_streams[barcode].open(outname);
        }
        // Write to the already open stream
        print_line(new_read_id, read_new, qual_new, barcode_file_streams[barcode]);
      } else {
        lock_guard<mutex> lock(cout_mutex);
        print_line(new_read_id, read_new, qual_new, std::cout);
      }
    }
}

// The consumer thread function: searches for barcodes in a single read
void search_worker(
    ThreadSafeQueue<SearchResult*> &work_queue,
    ThreadSafeQueue<SearchResult*> &results_queue,
    const unordered_set<string> &known_barcodes,
    int flank_edit_distance,
    int edit_distance,
    const std::vector<std::pair<std::string, std::string>> &search_pattern
) {
    SearchResult* sr_ptr;
    while (work_queue.pop(sr_ptr) && sr_ptr) {
        //forward search
        sr_ptr->vec_bc_for = big_barcode_search(sr_ptr->line, known_barcodes, flank_edit_distance, edit_distance, search_pattern);

        // get reverse complement
        sr_ptr->rev_line = sr_ptr->line;
        reverse_complement(sr_ptr->rev_line);

        //Check the reverse complement of the read
        sr_ptr->vec_bc_rev = big_barcode_search(sr_ptr->rev_line, known_barcodes, flank_edit_distance, edit_distance, search_pattern);

        sr_ptr->count = sr_ptr->vec_bc_for.size() + sr_ptr->vec_bc_rev.size();
        sr_ptr->chimeric = !sr_ptr->vec_bc_for.empty() && !sr_ptr->vec_bc_rev.empty();

        results_queue.push(sr_ptr);
    }
    // Push a sentinel to the results queue to signal this worker is done
    results_queue.push(nullptr);
}


// check if file already exists
bool file_exists(const std::string &filename) {
  std::ifstream infile(filename);
  return infile.good();
}

// MAIN is here!!
int main(int argc, char **argv) {
  std::ios_base::sync_with_stdio(false);

  cerr << "FLEXIPLEX " << VERSION << "\n";

  // Variables to store user options
  // Set these to their defaults
  int edit_distance = 2;       //(e)
  int flank_edit_distance = 8; //(f)

  // set the output filenames
  string out_stat_filename = "reads_barcodes.txt";
  string out_bc_filename = "barcodes_counts.txt";
  string out_filename_prefix = "flexiplex"; //(n)

  bool split_file_by_barcode = false; //(s)
  bool remove_barcodes = true;        //(i)
  bool print_chimeric = false;        //(c)

  std::vector<std::pair<std::string, std::string>> search_pattern;

  // Set of known barcodes
  unordered_set<string> known_barcodes;

  // threads
  int n_threads = 1;

  /*** Pass command line option *****/\
  int c;
  int params = 1;
  ifstream file;
  string line;

  vector<char *> myArgs(argv, argv + argc);
  vector<char *> allocated_args;
  

  while ((c = getopt(myArgs.size(), myArgs.data(),
                     "d:k:i:b:u:x:e:f:n:s:hp:c:")) != EOF) {
    switch (c) {
    case 'd': { // d=predefined list of settings for various search/barcode schemes
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
        char *newArg = new char[token.size() + 1]; // +1 for null terminator
        strcpy(newArg, token.c_str());
        newArgv.push_back(newArg); // Append the token to the new argv
      }
      params += 2;
      myArgs.insert(myArgs.begin() + params, newArgv.begin(), newArgv.end());
      allocated_args.insert(allocated_args.end(), newArgv.begin(), newArgv.end());
      break;
    }
    case 'k': { // k=list of known barcodes
      string file_name(optarg);
      string bc;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      if (!(file.good())) { // if the string given isn't a file
        cerr << "Reading known barcodes from string: " << file_name << "\n";
        stringstream bc_list(file_name);
        while (getline(bc_list, bc, ',')) // tokenize
          known_barcodes.insert(bc);
      } else {
        // otherwise get the barcodes from the file..
        cerr << "Setting known barcodes from " << file_name << "\n";
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
    case 'c': {
      print_chimeric = get_bool_opt_arg(optarg);
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

  // free memory from -d option
  for (auto arg : allocated_args) {
      delete[] arg;
  }

  // default case when no x, u, b is speficied
  if (search_pattern.empty()) {
    cerr << "Using default search pattern: " << endl;
    search_pattern = {{"primer", "CTACACGACGCTCTTCCGATCT"},
                      {"BC", std::string(16, '?')},
                      {"UMI", std::string(12, '?')},
                      {"polyA", std::string(9, 'T')}};
    for (auto const& i : search_pattern)
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
  atomic<int> r_count(0);
  atomic<int> bc_count(0);
  atomic<int> multi_bc_count(0);

  ofstream out_stat_file;
  out_stat_filename = out_filename_prefix + "_" + out_stat_filename;
  out_bc_filename = out_filename_prefix + "_" + out_bc_filename;

  if (known_barcodes.size() > 0) {
    if (file_exists(out_stat_filename)) {
      cerr << "File " << out_stat_filename
           << " already exists, overwriting." << endl;
    }
    out_stat_file.open(out_stat_filename, std::ios::trunc);
    out_stat_file << "Read\tCellBarcode\tFlankEditDist\tBarcodeEditDist\tUMI\n";
  }

  cerr << "Searching for barcodes...\n";

  bool is_fastq = true;
  unordered_map<string, int> barcode_counts;
  string read_id_line;

  if (getline(*in, read_id_line)) { // check the first line for file type
    if (read_id_line[0] == '>') {
      is_fastq = false;
    } else if (read_id_line[0] != '@') {
      cerr << "Unknown read format... exiting" << endl;
      exit(1);
    }
  } else { // Empty file
      return 0;
  }

  // --- Threading and Pipeline Setup ---
  ThreadSafeQueue<SearchResult*> work_queue;
  ThreadSafeQueue<SearchResult*> results_queue;
  vector<thread> workers;

  // Mutexes for synchronized output
  mutex cout_mutex;
  mutex file_mutex;
  map<string, ofstream> barcode_file_streams;

  // Start worker threads
  for (int i = 0; i < n_threads; ++i) {
      workers.emplace_back(search_worker, ref(work_queue), ref(results_queue), ref(known_barcodes),
                           flank_edit_distance, edit_distance, ref(search_pattern));
  }

  // --- Producer Stage ---
  // The main thread reads the file and pushes work to the queue
  do {
    SearchResult *sr = new SearchResult();
    istringstream line_stream(read_id_line);
    line_stream >> sr->read_id;
    sr->read_id.erase(0, 1);

    if (!is_fastq) { // fasta (account for multi-lines per read)
        string buffer_string;
        while (getline(*in, buffer_string) && buffer_string[0] != '>')
            sr->line += buffer_string;
        read_id_line = buffer_string;
    } else { // fastq (get quality scores)
        getline(*in, sr->line);
        getline(*in, line); // skip '+' line
        getline(*in, sr->qual_scores);
        getline(*in, read_id_line);
    }
    work_queue.push(sr);
  } while (*in);

  // Push sentinels to the work queue to signal completion
  for (int i = 0; i < n_threads; ++i) {
      work_queue.push(nullptr);
  }

  // --- Consumer/Writer Stage ---
  // The main thread now processes results
  int finished_workers = 0;
  while(finished_workers < n_threads) {
      SearchResult* sr_ptr;
      results_queue.pop(sr_ptr);

      if (sr_ptr == nullptr) {
          finished_workers++;
          continue;
      }
      
      unique_ptr<SearchResult> sr(sr_ptr); // Manage memory

      r_count++;
      if (sr->count > 0) bc_count++;
      if (sr->chimeric) multi_bc_count++;

      for (const auto& bc : sr->vec_bc_for) barcode_counts[bc.barcode]++;
      for (const auto& bc : sr->vec_bc_rev) barcode_counts[bc.barcode]++;

      if (known_barcodes.size() != 0) {
          print_stats(sr->read_id, sr->vec_bc_for, out_stat_file);
          print_stats(sr->read_id, sr->vec_bc_rev, out_stat_file);

          print_read(sr->read_id + "_+", sr->line, sr->qual_scores, sr->vec_bc_for, out_filename_prefix,
                     split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex, remove_barcodes, print_chimeric && sr->chimeric);

          if (remove_barcodes || sr->vec_bc_for.empty()) {
              reverse(sr->qual_scores.begin(), sr->qual_scores.end());
              print_read(sr->read_id + "_-", sr->rev_line, sr->qual_scores, sr->vec_bc_rev, out_filename_prefix,
                         split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex, remove_barcodes, print_chimeric && sr->chimeric);
          }
      }
      
      if ((r_count < 100000 && r_count % 10000 == 0) || (r_count % 100000 == 0)) {
          cerr << r_count / 1000000.0 << " million reads processed.." << endl;
      }
  }

  // Join all worker threads
  for (auto& t : workers) {
      t.join();
  }

  reads_ifs.close();
  if (known_barcodes.size() > 0) {
    out_stat_file.close();
  }
  // Close all the split files
  for (auto& pair : barcode_file_streams) {
      if (pair.second.is_open()) {
          pair.second.close();
      }
  }

  // print summary statistics
  cerr << "Number of reads processed: " << r_count << "\n";
  cerr << "Number of reads where a barcode was found: " << bc_count << "\n";
  cerr << "Number of reads where more than one barcode was found: "
       << multi_bc_count << "\n";
  cerr << "All done!" << endl;
  cerr << "If you like Flexiplex, please cite us! https://doi.org/10.1093/bioinformatics/btae102" << endl;

  if (known_barcodes.size() > 0) {
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

  if (!bc_vec.empty()) {
    vector<int> hist(bc_vec[0].second);
    ofstream out_bc_file;

    out_bc_file.open(out_bc_filename);

    for (auto const &bc_pair : bc_vec) {
        out_bc_file << bc_pair.first << "\t" << bc_pair.second << "\n";
        if (bc_pair.second > 0) {
            hist[bc_pair.second - 1]++;
        }
    }
    out_bc_file.close();

    cout << "Reads\tBarcodes" << "\n";
    for (int i = hist.size() - 1; i >= 0; i--)
        if(hist[i] > 0)
            cout << i + 1 << "\t" << hist[i] << "\n";
  }
}
