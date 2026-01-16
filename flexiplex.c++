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

enum SegmentType { FIXED, MATCHED, RANDOM };

struct Segment {
    SegmentType type;
    std::string pattern;
    std::string name; // e.g., "BC1", "UMI", "Primer1"
    std::string bc_list_name; // Only for MATCHED, refers to a key in known_barcodes_map
    int buffer_size;            // Only for MATCHED
    int max_edit_distance;      // Only for MATCHED
};

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
  cerr << "                     Reads are considered chimeric if barcodes are found on both\n";
  cerr << "                     forward and reverse strands.\n";
  // what about when only one strand has multiple barcodes?
  cerr << "     -n prefix       Prefix for output filenames.\n";
  cerr << "     -e N            Maximum edit distance to barcode (default 2).\n";
  cerr << "     -f N            Maximum edit distance to primer+polyT (default 8).\n";
  cerr << "     -p N            Number of threads (default: 1).\n\n";

  cerr << "  Specifying adaptor / barcode structure : \n";
  cerr << "     -x sequence Append flanking sequence to search for\n";
  cerr << "     -b sequence Append the barcode pattern to search for (named CB by default)\n";
  cerr << "     -u sequence[,name:<name>] Append the UMI pattern to search for\n";
  cerr << "                 (named UB, UB2, UB3... by default if name not provided)\n";
  cerr << "     -B sequence[,name:<name>][,list:<file>][,buffer:<size>][,max_ed:<dist>]\n";
  cerr << "                 Append a barcode pattern with optional parameters\n";
  cerr << "                 - name: specify the name of the barcode segment (defaults to BC1, BC2...)\n";
  cerr << "                 - list: specify a file containing expected barcodes for this segment\n";
  cerr << "                 - buffer: number of bases to extend search region on either side (default 5)\n";
  cerr << "                 - max_ed: maximum edit distance for this segment (default is global -e value)\n";
  cerr << "     When no search pattern x,b,u option is provided, the following default pattern is used: \n";
  cerr << "          primer: CTACACGACGCTCTTCCGATCT\n";
  cerr << "          barcode (CB): ????????????????\n";
  cerr << "          UMI (UB): ????????????\n";
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
  std::map<std::string, std::string> features; // Map of segment name to extracted sequence
  int flank_editd;
  int flank_start;
  int flank_end;
  bool found_all_matched_segments;
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



// Given a string 'seq' search for substring with primer and polyT sequence followed by
// a targeted search in the region for barcode
// Sequence seearch is performed using edlib

Barcode extract_features(string & seq,
		    const std::vector<Segment> &segments,
		    const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
		    int global_flank_max_editd,
		    int global_barcode_edit_distance) {

  //initialise struct variables for return:
  Barcode barcode_result;
  barcode_result.found_all_matched_segments = true;
  barcode_result.flank_editd=100;

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
  EdlibAlignConfig edlibConf = {global_flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 28};

  // concatenate patterns in search_pattern in insertion order
  std::string search_string_template;
  std::vector<long unsigned int> segment_lengths;
  for (const auto &segment : segments) {
    search_string_template += segment.pattern;
    segment_lengths.push_back(segment.pattern.length());
  }

  // shorter than the pattern, skip search
  if (seq.length() < search_string_template.length()) {
    return barcode_result;
  }

  //search for the concatenated pattern
  EdlibAlignResult result = edlibAlign(search_string_template.c_str(), search_string_template.length(), seq.c_str(), seq.length(), edlibConf);
  if(result.status != EDLIB_STATUS_OK || result.numLocations==0 ){
    edlibFreeAlignResult(result);
    return(barcode_result); // no match found - return
  } //fill in info about found primer and polyT location
  barcode_result.flank_editd=result.editDistance;
  barcode_result.flank_start=result.startLocations[0];
  barcode_result.flank_end=result.endLocations[0];

  std::vector<long unsigned int> segment_template_ends;
  segment_template_ends.resize(segment_lengths.size());
  std::partial_sum(segment_lengths.begin(), segment_lengths.end(), segment_template_ends.begin());

  vector<int> read_to_segment_starts;
  read_to_segment_starts.reserve(segment_template_ends.size() + 1);
  read_to_segment_starts.emplace_back(barcode_result.flank_start);

  // initialise pointers
  int i_read = barcode_result.flank_start;
  int i_pattern = 0;
  int i_segment = 0;

  // walk through edlib aligment
  // 0 for match
  // 1 for insertion to target (read)
  // 2 for insertion to query (template)
  // 3 for mismatch
  std::vector<unsigned char> alignment_vector(result.alignment, result.alignment + result.alignmentLength);
  for (const auto& value : alignment_vector) {
    if (value != EDLIB_EDOP_INSERT) { // If not an insertion in the read
      i_read++;
    }
    if (value != EDLIB_EDOP_DELETE) { // If not an insertion in the template
      i_pattern++;
    }
    if (i_segment < segment_template_ends.size() && i_pattern >= segment_template_ends[i_segment]) {
      read_to_segment_starts.emplace_back(i_read);
      i_segment++;
    }
  }
  edlibFreeAlignResult(result);

  // Storage for refined positions from matched segments
  std::vector<int> refined_segment_starts(segments.size(), -1);
  std::vector<int> refined_segment_ends(segments.size(), -1);

  // Pass 1: Process MATCHED segments
  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment& s = segments[i];
    if (s.type != MATCHED) continue;

    int segment_read_start = read_to_segment_starts[i];
    int segment_read_end = (i + 1 < read_to_segment_starts.size()) ? read_to_segment_starts[i+1] : barcode_result.flank_end;

    const unordered_set<string>* current_bclist = nullptr;

    if (!s.bc_list_name.empty()) {
        auto it = known_barcodes_map->find(s.bc_list_name);
        if (it != known_barcodes_map->end()) {
            current_bclist = &(it->second);
        }
    } else if (known_barcodes_map->count("global")) {
        current_bclist = &(known_barcodes_map->at("global"));
    }

    if (current_bclist && !current_bclist->empty()) {
      // Apply buffer for search
      int search_start = max(0, segment_read_start - s.buffer_size);
      int search_end = min((int)seq.length(), segment_read_end + s.buffer_size);
      std::string search_region = seq.substr(search_start, search_end - search_start);

      unsigned int best_edit_distance = 100; // Higher than any reasonable max_ed
      unsigned int best_end_distance = 0;
      std::string best_barcode_match = "";
      bool current_segment_unambiguous = false;

      for (const auto& known_bc : *current_bclist) {
        unsigned int editDistance, endDistance;
        editDistance = edit_distance(search_region, known_bc, endDistance, s.max_edit_distance);

        if (editDistance == best_edit_distance) {
          current_segment_unambiguous = false;
        } else if (editDistance < best_edit_distance && editDistance <= s.max_edit_distance) {
          current_segment_unambiguous = true;
          best_edit_distance = editDistance;
          best_barcode_match = known_bc;
          best_end_distance = endDistance;

          // TODO: multiplex perfect match can exits due to the buffering around barcode segments
          if (editDistance == 0) {
            break; // perfect match found
          }
        }
      }
      if (best_edit_distance <= s.max_edit_distance && current_segment_unambiguous) {
          barcode_result.features[s.name] = best_barcode_match;
          
          // Refined positions
          refined_segment_ends[i] = search_start + best_end_distance - 1;
          // Approximate start based on match length (assuming not many indels)
          refined_segment_starts[i] = refined_segment_ends[i] - best_barcode_match.length() + 1;
      } else {
          barcode_result.found_all_matched_segments = false;
      }
    } else { // No barcode list found
      cerr << "Error: No barcode list found for segment " << s.name << ".\n";
      exit(1);
    }
  }

  if (!barcode_result.found_all_matched_segments) {
    return(barcode_result); // return early if any matched segments failed
  }

  // Pass 2: Process RANDOM and other segments
  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment& s = segments[i];
    if (s.type == MATCHED) continue;

    int extract_start = read_to_segment_starts[i];
    int extract_end = (i + 1 < read_to_segment_starts.size()) ? read_to_segment_starts[i+1] : barcode_result.flank_end;

    if (s.type == RANDOM) {
        // Check if there is a MATCHED segment immediately before
        if (i > 0 && segments[i-1].type == MATCHED && refined_segment_ends[i-1] != -1) {
            extract_start = refined_segment_ends[i-1] + 1;
            extract_end = extract_start + s.pattern.length();
        } 
        // Otherwise check if there is one immediately after
        else if (i + 1 < segments.size() && segments[i+1].type == MATCHED && refined_segment_starts[i+1] != -1) {
            extract_end = refined_segment_starts[i+1];
            extract_start = extract_end - s.pattern.length();
        }

        int actual_extract_start = max(0, extract_start);
        int actual_extract_end = min((int)seq.length(), extract_end);
        
        if (actual_extract_end < actual_extract_start) actual_extract_end = actual_extract_start;
        
        std::string extracted_seq = seq.substr(actual_extract_start, actual_extract_end - actual_extract_start);
        barcode_result.features[s.name] = extracted_seq;
    }
  }

  return(barcode_result); //return the best matched barcode and associated information
}

vector<Barcode> big_barcode_search(const string & sequence,
                                   const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                                   int global_flank_max_editd,
                                   int global_barcode_edit_distance,
                                   const std::vector<Segment> &segments) {
    vector<Barcode> return_vec;
    string masked_sequence = sequence; // Work on a copy

    while (true) {
        Barcode result = extract_features(masked_sequence, segments, known_barcodes_map, global_flank_max_editd, global_barcode_edit_distance);
        if (result.flank_editd <= global_flank_max_editd && result.found_all_matched_segments) {
            return_vec.push_back(result);

            // Mask the found region to prevent re-finding it
            // result.flank_end is inclusive, so length is end - start + 1
            int match_length = result.flank_end - result.flank_start + 1;

            if (match_length > 0) {
                masked_sequence.replace(result.flank_start, match_length, string(match_length, 'X'));
            } else {
                break; // Should not happen for valid match, but prevents infinite loop
            }
        } else {
            break;
        }
    }
    return return_vec;
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
void print_stats(const string& read_id, const vector<Barcode> & vec_bc, ostream & out_stream, const std::vector<Segment> &segments){
  for(int b=0; b<vec_bc.size() ; b++){
    out_stream << read_id;
    for (const auto& s : segments) {
        if (s.type != FIXED) {
            auto it = vec_bc.at(b).features.find(s.name);
            if (it != vec_bc.at(b).features.end()) {
                out_stream << "\t" << it->second;
            } else {
                out_stream << "\t";
            }
        }
    }
    out_stream << "\t" << vec_bc.at(b).flank_editd << "\n";
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
    bool chimeric,
    const std::vector<Segment> &segments){

    auto vec_size = vec_bc.size();

    for (int b = 0; b < vec_size; b++) {
      stringstream ss_id_suffix;
      ss_id_suffix << (b + 1) << "of" << vec_size;
      if (chimeric) {
        ss_id_suffix << "_" << "C";
      }

      stringstream new_read_id_ss;
      string primary_barcode_for_filename = "";

      for (const auto& s : segments) {
          if (s.type != FIXED) {
              auto it = vec_bc.at(b).features.find(s.name);
              if (it != vec_bc.at(b).features.end()) {
                  new_read_id_ss << it->second << "_";
                  if (s.type == MATCHED && primary_barcode_for_filename.empty()) {
                      primary_barcode_for_filename = it->second;
                  }
              } else {
                  new_read_id_ss << "NA_";
              }
          }
      }
      string new_read_id_prefix = new_read_id_ss.str();
      if (!new_read_id_prefix.empty() && new_read_id_prefix.back() == '_') {
          new_read_id_prefix.pop_back();
      }

      string new_read_id = new_read_id_prefix + "#" + read_id + ss_id_suffix.str();

      for (const auto& s : segments) {
          if (s.type != FIXED) {
              auto it = vec_bc.at(b).features.find(s.name);
              if (it != vec_bc.at(b).features.end()) {
                  new_read_id += "\t" + s.name + ":Z:" + it->second;
              }
          }
      }

      int read_start = vec_bc.at(b).flank_end + 1;
      // work out the start and end base in case multiple barcodes
      int next_barcode_start = read.length();
      for (int k = 0; k < vec_size; k++) {
          if (b == k) continue;
          if (vec_bc.at(k).flank_start >= read_start) {
              if (vec_bc.at(k).flank_start < next_barcode_start) {
                  next_barcode_start = vec_bc.at(k).flank_start;
              }
          }
      }

      int read_length = next_barcode_start - read_start;

      string qual_new = "";
      if (qual != "") {
        if((read_start+read_length)>(qual.length())){
            lock_guard<mutex> lock(cout_mutex);
            cerr << "WARNING: sequence and quality lengths diff for read: " << read_id << ". Ignoring read." << endl;
            return;
        }
        qual_new = qual.substr(read_start, read_length);
      }
      string read_new = read.substr(read_start, read_length);

      if (b == 0 && !trim_barcodes) {
        new_read_id = read_id;
        read_new = read;
        qual_new = qual;
      }

      if (read_new.length() == 0) {
        continue;
      }

      if (split) {
        lock_guard<mutex> lock(file_mutex);
        string barcode_for_split = primary_barcode_for_filename.empty() ? "UNKNOWN" : primary_barcode_for_filename;
        if (barcode_file_streams.find(barcode_for_split) == barcode_file_streams.end()) {
            string outname = prefix + "_" + barcode_for_split + (qual.empty() ? ".fasta" : ".fastq");
            barcode_file_streams[barcode_for_split].open(outname);
        }
        print_line(new_read_id, read_new, qual_new, barcode_file_streams[barcode_for_split]);
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
    const std::unordered_map<std::string, std::unordered_set<std::string>> &known_barcodes_map,
    int flank_edit_distance,
    int global_barcode_edit_distance,
    const std::vector<Segment> &search_pattern
) {
    SearchResult* sr_ptr;
    while (work_queue.pop(sr_ptr) && sr_ptr) {
        //forward search
        sr_ptr->vec_bc_for = big_barcode_search(sr_ptr->line, &known_barcodes_map, flank_edit_distance, global_barcode_edit_distance, search_pattern);

        // get reverse complement
        sr_ptr->rev_line = sr_ptr->line;
        reverse_complement(sr_ptr->rev_line);

        //Check the reverse complement of the read
        sr_ptr->vec_bc_rev = big_barcode_search(sr_ptr->rev_line, &known_barcodes_map, flank_edit_distance, global_barcode_edit_distance, search_pattern);

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

  std::vector<Segment> search_pattern;

  // Set of known barcodes (can be multiple, identified by name)
  std::unordered_map<std::string, std::unordered_set<std::string>> known_barcodes_map;
  // Flag to track if any -B segment specified its own barcode list
  bool individual_bclist_specified = false;

  // Counters for automatic naming
  int barcode_count = 0;
  int umi_count = 0;

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
                     "d:k:i:b:u:x:e:f:n:s:hp:c:B:")) != EOF) {
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
    case 'k': { // k=list of known barcodes (global barcode list)
      if (individual_bclist_specified) {
        cerr << "Error: Cannot use -k when -B segments have specified their own barcode files.\n";
        print_usage();
        exit(1);
      }
      string file_name(optarg);
      string bc;
      unordered_set<string> global_bclist;
      /**** READ BARCODE LIST FROM FILE ******/
      file.open(file_name);
      if (!(file.good())) { // if the string given isn't a file
        cerr << "Reading known barcodes from string: " << file_name << "\n";
        stringstream bc_list(file_name);
        while (getline(bc_list, bc, ',')) // tokenize
          global_bclist.insert(bc);
      } else {
        // otherwise get the barcodes from the file..
        cerr << "Setting known barcodes from " << file_name << "\n";
        while (getline(file, line)) {
          istringstream line_stream(line);
          line_stream >> bc;
          global_bclist.insert(bc);
        }
        file.close();
      }
      cerr << "Number of known barcodes in global list: " << global_bclist.size() << "\n";
      if (global_bclist.empty()) {
        print_usage();
        exit(1); // case barcode file is empty
      }
      known_barcodes_map["global"] = global_bclist;
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
      cerr << "Setting max barcode edit distance to " << edit_distance << " (global default).\n";
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
    case 'x': { // Fixed segment
      Segment s;
      s.type = FIXED;
      s.pattern = optarg;
      s.name = "Fixed_" + to_string(search_pattern.size());
      search_pattern.push_back(s);
      cerr << "Adding fixed sequence to search for: " << optarg << "\n";
      params += 2;
      break;
    }
    case 'u': { // Random (UMI) segment
      Segment s;
      s.type = RANDOM;
      s.name = ""; // Will be set later if not provided
      s.buffer_size = 0; // Default buffer for UMI

      string arg_str(optarg);
      size_t current_pos = 0;
      size_t next_comma = arg_str.find(',');

      // First part is always the pattern
      s.pattern = arg_str.substr(current_pos, next_comma - current_pos);
      cerr << "Setting UMI to search for: " << s.pattern;

      while (next_comma != string::npos) {
          current_pos = next_comma + 1;
          next_comma = arg_str.find(',', current_pos);
          string param = arg_str.substr(current_pos, next_comma - current_pos);

          size_t colon_pos = param.find(':');
          if (colon_pos == string::npos) {
              cerr << "Error: Malformed -u parameter: " << param << "\n";
              print_usage();
              exit(1);
          }
          string key = param.substr(0, colon_pos);
          string value = param.substr(colon_pos + 1);

          if (key == "name") {
              s.name = value;
              cerr << " (name: " << s.name << ")";
          } else if (key == "buffer") {
              s.buffer_size = stoi(value);
              cerr << " (buffer: " << s.buffer_size << ")";
          } else {
              cerr << "Error: Unknown -u parameter: " << key << "\n";
              print_usage();
              exit(1);
          }
      }
      
      // Assign default name if not provided
      if (s.name.empty()) {
          umi_count++;
          if (umi_count == 1) {
              s.name = "UB";
          } else {
              s.name = "UB" + to_string(umi_count);
          }
      } else {
          umi_count++; // Still increment counter for user-named UMIs
      }
      
      search_pattern.push_back(s);
      cerr << " -> named: " << s.name << "\n";
      params += 2;
      break;
    }
    case 'b': { // Old-style Barcode segment
      Segment s;
      s.type = MATCHED;
      s.pattern = optarg;
      barcode_count++;
      s.name = (barcode_count == 1) ? "CB" : "CB" + to_string(barcode_count);
      s.buffer_size = 5; // Default buffer for barcode
      s.max_edit_distance = edit_distance; // Will use global -e
      s.bc_list_name = "global"; // Will use global -k
      search_pattern.push_back(s);
      cerr << "Setting old-style barcode to search for: " << optarg 
           << " -> named: " << s.name 
           << ". Will use global barcode list (-k) and edit distance (-e).\n";
      params += 2;
      break;
    }
    case 'B': { // New-style Barcode segment
      Segment s;
      s.type = MATCHED;
      s.name = ""; // Will be set later if not provided
      s.buffer_size = 5; // Default buffer
      s.max_edit_distance = edit_distance; // Default to global -e

      string arg_str(optarg);
      size_t current_pos = 0;
      size_t next_comma = arg_str.find(',');

      // First part is always the pattern
      s.pattern = arg_str.substr(current_pos, next_comma - current_pos);
      cerr << "Adding new-style barcode segment: " << s.pattern;

      while (next_comma != string::npos) {
          current_pos = next_comma + 1;
          next_comma = arg_str.find(',', current_pos);
          string param = arg_str.substr(current_pos, next_comma - current_pos);

          size_t colon_pos = param.find(':');
          if (colon_pos == string::npos) {
              cerr << "Error: Malformed -B parameter: " << param << "\n";
              print_usage();
              exit(1);
          }
          string key = param.substr(0, colon_pos);
          string value = param.substr(colon_pos + 1);

          if (key == "name") {
              s.name = value;
              cerr << " (name: " << s.name << ")";
          } else if (key == "list") {
              s.bc_list_name = value;
              individual_bclist_specified = true;
              // Load barcode file
              unordered_set<string> current_bclist;
              file.open(value);
              if (!(file.good())) {
                  cerr << "Error: Unable to open barcode file: " << value << "\n";
                  print_usage();
                  exit(1);
              }
              while (getline(file, line)) {
                  istringstream line_stream(line);
                  string bc_entry;
                  line_stream >> bc_entry;
                  current_bclist.insert(bc_entry);
              }
              file.close();
              if (current_bclist.empty()) {
                  cerr << "Error: barcode file " << value << " is empty.\n";
                  print_usage();
                  exit(1);
              }
              known_barcodes_map[value] = current_bclist;
              cerr << " (barcode list: " << value << ")";
          } else if (key == "buffer") {
              s.buffer_size = stoi(value);
              cerr << " (buffer: " << s.buffer_size << ")";
          } else if (key == "max_ed") {
              s.max_edit_distance = stoi(value);
              cerr << " (max_ed: " << s.max_edit_distance << ")";
          } else {
              cerr << "Error: Unknown -B parameter: " << key << "\n";
              print_usage();
              exit(1);
          }
      }
      
      // Assign default name if not provided
      if (s.name.empty()) {
          barcode_count++;
          if (barcode_count == 1) {
              s.name = "CB";
          } else {
              s.name = "CB" + to_string(barcode_count);
          }
      } else {
          barcode_count++; // Still increment counter for user-named barcodes
      }
      
      search_pattern.push_back(s);
      cerr << " -> named: " << s.name << "\n";
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

  // default case when no x, u, b, B is speficied
  if (search_pattern.empty()) {
    cerr << "Using default search pattern: " << endl;
    Segment s_primer, s_bc, s_umi, s_polya;

    s_primer.type = FIXED;
    s_primer.pattern = "CTACACGACGCTCTTCCGATCT";
    s_primer.name = "primer";
    search_pattern.push_back(s_primer);

    s_bc.type = MATCHED;
    s_bc.pattern = std::string(16, '?');
    s_bc.name = "CB";
    s_bc.buffer_size = 5;
    s_bc.max_edit_distance = edit_distance;
    s_bc.bc_list_name = "global";
    search_pattern.push_back(s_bc);

    s_umi.type = RANDOM;
    s_umi.pattern = std::string(12, '?');
    s_umi.name = "UB";
    s_umi.buffer_size = 0;
    search_pattern.push_back(s_umi);

    s_polya.type = FIXED;
    s_polya.pattern = std::string(9, 'T');
    s_polya.name = "polyA";
    search_pattern.push_back(s_polya);

    for (auto const& s : search_pattern)
      std::cerr << s.name << ": " << s.pattern << "\n";
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

  if (!known_barcodes_map.empty()) {
    if (file_exists(out_stat_filename)) {
      cerr << "File " << out_stat_filename
           << " already exists, overwriting." << endl;
    }
    out_stat_file.open(out_stat_filename, std::ios::trunc);
    out_stat_file << "Read";
    for (const auto& s : search_pattern) {
        if (s.type != FIXED) {
            out_stat_file << "\t" << s.name;
        }
    }
    out_stat_file << "\tFlankEditDist\n";
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
      workers.emplace_back(search_worker, ref(work_queue), ref(results_queue), ref(known_barcodes_map),
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

      // Count barcodes from features map (first MATCHED segment)
      for (const auto& bc : sr->vec_bc_for) {
          for (const auto& s : search_pattern) {
              if (s.type == MATCHED) {
                  auto it = bc.features.find(s.name);
                  if (it != bc.features.end()) {
                      barcode_counts[it->second]++;
                      break; // Only count first barcode for statistics
                  }
              }
          }
      }
      for (const auto& bc : sr->vec_bc_rev) {
          for (const auto& s : search_pattern) {
              if (s.type == MATCHED) {
                  auto it = bc.features.find(s.name);
                  if (it != bc.features.end()) {
                      barcode_counts[it->second]++;
                      break; // Only count first barcode for statistics
                  }
              }
          }
      }

      if (known_barcodes_map.size() != 0) {
          print_stats(sr->read_id, sr->vec_bc_for, out_stat_file, search_pattern);
          print_stats(sr->read_id, sr->vec_bc_rev, out_stat_file, search_pattern);

          print_read(sr->read_id + "_+", sr->line, sr->qual_scores, sr->vec_bc_for, out_filename_prefix,
                     split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex, remove_barcodes, print_chimeric && sr->chimeric, search_pattern);

          if (remove_barcodes || sr->vec_bc_for.empty()) {
              reverse(sr->qual_scores.begin(), sr->qual_scores.end());
              print_read(sr->read_id + "_-", sr->rev_line, sr->qual_scores, sr->vec_bc_rev, out_filename_prefix,
                         split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex, remove_barcodes, print_chimeric && sr->chimeric, search_pattern);
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
  if (known_barcodes_map.size() > 0) {
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

  if (known_barcodes_map.size() > 0) {
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
