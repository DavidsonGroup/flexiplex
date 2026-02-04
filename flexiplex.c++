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
const static string VERSION="1.02.6";

enum SegmentType { FIXED, MATCHED, MATCHED_SPLIT, RANDOM };

struct Segment {
    SegmentType type;
    std::string pattern;
    std::string name; // e.g., "BC1", "UMI", "Primer1"
    std::string bc_list_name; // Only for MATCHED, refers to a key in known_barcodes_map
    int buffer_size;            // for MATCHED, MATCHED_SPLIT
    int max_edit_distance;      // for MATCHED, MATCHED_SPLIT
};

struct BarcodeGroup {
    std::string name;
    int max_edit_distance;
    std::vector<size_t> segment_indices;
};

// can have other properties
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
  cerr << "     -u sequence[,name:<name>] Append the UMI pattern to search for\n";
  cerr << "                 (named UB, UB2, UB3... by default if name not provided)\n";
  cerr << "     -b sequence[,name:<name>][,list:<file>][,group:<name>][,buffer:<size>][,max_ed:<dist>]\n";
  cerr << "                 Append a barcode pattern with optional parameters\n";
  cerr << "                 - name: specify the name of the barcode segment (defaults to BC1, BC2...)\n";
  cerr << "                 - list: specify a file containing expected barcodes for this segment\n";
  cerr << "                 - group: specify the group name for this segment\n";
  cerr << "                   Barcodes in the same group will be treated as parts of a single barcode\n";
  cerr << "                 - buffer: number of bases to extend search region on either side (default 5)\n";
  cerr << "                 - max_ed: maximum edit distance for this segment (default is global -e value)\n";
  cerr << "     -g name,list:<file>[,max_ed:<dist>]\n";
  cerr << "                 Define a barcode group (multiple barcode segments joined together)\n";
  cerr << "                 - name: unique name for the group\n";
  cerr << "                 - list: specify a file containing expected barcodes for this group\n";
  cerr << "                 - max_ed: maximum edit distance allowed for the entire group\n";
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

// Holds the found barcode and associated information
struct Barcode {
  std::map<std::string, std::string> features; // Map of segment name to extracted sequence
  int flank_editd = 100;
  int flank_start = -1;
  int flank_end = -1; // inclusive
  bool found_all_matched_segments = true;
};

// A thread-safe queue for our producer-consumer model
template <typename T>
class ThreadSafeQueue {
private:
    std::queue<T> queue_;
    mutable std::mutex mutex_;
    std::condition_variable cond_;
    std::condition_variable cond_not_full_;
    size_t capacity_ = 0; // 0 => unbounded

public:
    ThreadSafeQueue() = default;
    explicit ThreadSafeQueue(size_t capacity) : capacity_(capacity) {}

    void push(T value) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (capacity_ > 0) {
            cond_not_full_.wait(lock, [this] { return queue_.size() < capacity_; });
        }
        queue_.push(std::move(value));
        lock.unlock();
        cond_.notify_one();
    }

    // Pop an element from the queue. Blocks until an element is available.
    // Returns true if a value was popped.
    bool pop(T& value) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_.wait(lock, [this] { return !queue_.empty(); });

        value = std::move(queue_.front());
        queue_.pop();
        lock.unlock();

        if (capacity_ > 0) {
            cond_not_full_.notify_one();
        }
        return true;
    }

    void notify_all_finish() {
        std::lock_guard<std::mutex> lock(mutex_);
        cond_.notify_all();
        cond_not_full_.notify_all();
    }
};

struct SearchResult {
  string read_id;
  string qual_scores;
  string line;
  string rev_line;
  vector<Barcode> vec_bc_for;
  vector<Barcode> vec_bc_rev;
  int count = 0;
  bool chimeric = false;
};

// Forward declaration
vector<Barcode> big_barcode_search(
  const string &sequence,
  const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
  const std::unordered_map<std::string, BarcodeGroup> *group_map,
  int global_flank_max_editd,
  int global_barcode_edit_distance,
  const std::vector<Segment> &segments);

// Separated out from main so it can be run on batches in threads (deterministic output)
static void search_read_batch(
    std::vector<SearchResult> &reads,
    const std::unordered_map<std::string, std::unordered_set<std::string>> &known_barcodes_map,
    const std::unordered_map<std::string, BarcodeGroup> &group_map,
    int flank_edit_distance,
    int global_barcode_edit_distance,
    const std::vector<Segment> &search_pattern) {

  for (auto &sr : reads) {
    // forward search
    sr.vec_bc_for = big_barcode_search(sr.line, &known_barcodes_map, &group_map,
                                      flank_edit_distance,
                                      global_barcode_edit_distance,
                                      search_pattern);

    // reverse complement
    sr.rev_line = sr.line;
    reverse_complement(sr.rev_line);

    // reverse search
    sr.vec_bc_rev = big_barcode_search(sr.rev_line, &known_barcodes_map, &group_map,
                                      flank_edit_distance,
                                      global_barcode_edit_distance,
                                      search_pattern);

    sr.count = (int)sr.vec_bc_for.size() + (int)sr.vec_bc_rev.size();
    sr.chimeric = !sr.vec_bc_for.empty() && !sr.vec_bc_rev.empty();
  }
}

// Code for fast edit distance calculation for short sequences modified from
// https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#C++
// s2 is always assumned to be the shorter string (barcode)
unsigned int edit_distance(const std::string& s1, const std::string& s2, unsigned int &end, int max_editd){

  const std::size_t len1 = s1.size() + 1, len2 = s2.size() + 1;
  const char *s1_c = s1.c_str(); const char *s2_c = s2.c_str();

  // Reuse DP buffer per-thread to reduce allocator churn.
  thread_local std::vector<unsigned int> dist_holder;
  dist_holder.assign(len1 * len2, 0u);
  //initialise the edit distance matrix.
  //penalise for gaps at the start and end of the shorter sequence (j)
  //but not for shifting the start/end of the longer sequence (i,0)
  dist_holder[0]=0; //[0][0]
  for(unsigned int j = 1; j < len2; ++j) dist_holder[j] = j; //[0][j];
  for(unsigned int i = 1; i < len1; ++i) dist_holder[i*len2] = 0; //[i][0];

  int best = len2;
  end = len1 - 1;

  //loop over the distance matrix elements and calculate running distance
  for (unsigned int j = 1; j < len2; ++j) {
    bool any_below_threshold = false; // flag used for early exit
    for (unsigned int i = 1; i < len1; ++i) {
      const int sub = (s1_c[i - 1] == s2_c[j - 1]) ? 0 : 1; // are the bases the same?

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

// Forward declarations
bool align_read_to_pattern(const string &seq,
                           const std::vector<Segment> &segments,
                           int global_flank_max_editd,
                           Barcode &barcode_result,
                           std::vector<int> &read_to_segment_starts);
void refine_matched_segments(const string &seq,
                             const std::vector<Segment> &segments,
                             const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                             const std::vector<int> &read_to_segment_starts,
                             Barcode &barcode_result,
                             std::vector<int> &refined_segment_starts,
                             std::vector<int> &refined_segment_ends);
void refine_split_segments(const string &seq,
                           const std::vector<Segment> &segments,
                           const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                           const std::unordered_map<std::string, BarcodeGroup> *group_map,
                           const std::vector<int> &read_to_segment_starts,
                           Barcode &barcode_result,
                           std::vector<int> &refined_segment_starts,
                           std::vector<int> &refined_segment_ends);
void extract_random_segments(const string &seq,
                             const std::vector<Segment> &segments,
                             const std::vector<int> &read_to_segment_starts,
                             const std::vector<int> &refined_segment_starts,
                             const std::vector<int> &refined_segment_ends,
                             Barcode &barcode_result);

// Given a string 'seq' search for substring with primer and polyT sequence followed by
// a targeted search in the region for barcode
// Sequence seearch is performed using edlib

Barcode extract_features(string & seq,
                         const std::vector<Segment> &segments,
                         const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                         const std::unordered_map<std::string, BarcodeGroup> *group_map,
                         int global_flank_max_editd,
                         int global_barcode_edit_distance) {

  Barcode barcode_result;
  barcode_result.found_all_matched_segments = true;
  barcode_result.flank_editd = 100;

  vector<int> read_to_segment_starts;

  // 1. Align reads to get approximate positions
  bool aligned = align_read_to_pattern(seq, segments, global_flank_max_editd, barcode_result, read_to_segment_starts);

  if (!aligned) return barcode_result;

  // Storage for refined positions from matched segments
  std::vector<int> refined_segment_starts(segments.size(), -1);
  std::vector<int> refined_segment_ends(segments.size(), -1);

  // 2. Process MATCHED segments (Single barcodes)
  refine_matched_segments(seq, segments, known_barcodes_map, read_to_segment_starts, barcode_result, refined_segment_starts, refined_segment_ends);

  // 3. Process MATCHED_SPLIT segments (Grouped/Split barcodes)
  if (group_map && !group_map->empty()) {
    refine_split_segments(seq, segments, known_barcodes_map, group_map, read_to_segment_starts, barcode_result, refined_segment_starts, refined_segment_ends);
  }

  if (!barcode_result.found_all_matched_segments) {
    return barcode_result;
  }

  // 4. Process RANDOM segments using refined anchors
  extract_random_segments(seq, segments, read_to_segment_starts, refined_segment_starts, refined_segment_ends, barcode_result);

  return barcode_result;
}

// Helper: Align read to the concatenated pattern using Edlib
bool align_read_to_pattern(const string &seq,
                           const std::vector<Segment> &segments,
                           int global_flank_max_editd,
                           Barcode &barcode_result,
                           std::vector<int> &read_to_segment_starts) {

  // Edlib Configuration 
  EdlibEqualityPair additionalEqualities[28] = {
    {'R', 'A'}, {'R', 'G'}, {'K', 'G'}, {'K', 'T'},
    {'S', 'G'}, {'S', 'C'}, {'Y', 'C'}, {'Y', 'T'},
    {'M', 'A'}, {'M', 'C'}, {'W', 'A'}, {'W', 'T'},
    {'B', 'C'}, {'B', 'G'}, {'B', 'T'}, {'H', 'A'}, {'H', 'C'}, {'H', 'T'},
    {'?', 'A'}, {'?', 'C'}, {'?', 'G'}, {'?', 'T'},
    {'D', 'A'}, {'D', 'G'}, {'D', 'T'}, {'V', 'A'}, {'V', 'C'}, {'V', 'G'}
  };
  EdlibAlignConfig edlibConf = {global_flank_max_editd, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 28};

  // 1. Construct Search Template
  // Optimization Note: This construction happens every call. If segments are static, 
  // the template string and lengths could be pre-calculated.
  std::string search_string_template;
  std::vector<long unsigned int> segment_lengths;

  search_string_template.reserve(256);
  segment_lengths.reserve(segments.size());

  for (const auto &segment : segments) {
    search_string_template += segment.pattern;
    segment_lengths.push_back(segment.pattern.length());
  }

  if (seq.length() < search_string_template.length()) return false;

  // 2. Perform Alignment
  EdlibAlignResult result = edlibAlign(search_string_template.c_str(), search_string_template.length(), seq.c_str(), seq.length(), edlibConf);

  if (result.status != EDLIB_STATUS_OK || result.numLocations == 0) {
    edlibFreeAlignResult(result);
    return false;
  }

  barcode_result.flank_editd = result.editDistance;
  barcode_result.flank_start = result.startLocations[0];
  barcode_result.flank_end = result.endLocations[0];

  // 3. Map alignment to segment positions
  std::vector<long unsigned int> segment_template_ends;
  segment_template_ends.resize(segment_lengths.size());
  std::partial_sum(segment_lengths.begin(), segment_lengths.end(), segment_template_ends.begin());

  read_to_segment_starts.reserve(segment_template_ends.size() + 1);
  read_to_segment_starts.emplace_back(barcode_result.flank_start);

  int i_read = barcode_result.flank_start;
  int i_pattern = 0;
  size_t i_segment = 0;

  // Avoid copying the alignment to a temporary vector.
  for (int ai = 0; ai < result.alignmentLength; ++ai) {
    const unsigned char value = result.alignment[ai];
    if (value != EDLIB_EDOP_INSERT) i_read++;
    if (value != EDLIB_EDOP_DELETE) i_pattern++;

    if (i_segment < segment_template_ends.size() &&
        i_pattern >= (int)segment_template_ends[i_segment]) {
      read_to_segment_starts.emplace_back(i_read);
      i_segment++;
    }
  }

  edlibFreeAlignResult(result);
  return true;
}

// Helper: Process direct MATCHED segments (Pass 1)
void refine_matched_segments(const string &seq,
                             const std::vector<Segment> &segments,
                             const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                             const std::vector<int> &read_to_segment_starts,
                             Barcode &barcode_result,
                             std::vector<int> &refined_segment_starts,
                             std::vector<int> &refined_segment_ends) {

  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment& s = segments[i];
    if (s.type != MATCHED) continue;

    const int segment_read_start = read_to_segment_starts[i];
    const int segment_read_end =
      (i + 1 < read_to_segment_starts.size())
      ? read_to_segment_starts[i+1]
      : barcode_result.flank_end;

    const unordered_set<string>* current_bclist = nullptr;
    if (!s.bc_list_name.empty()) {
      auto it = known_barcodes_map->find(s.bc_list_name);
      if (it != known_barcodes_map->end()) current_bclist = &(it->second);
    } else if (known_barcodes_map->count("global")) {
      current_bclist = &(known_barcodes_map->at("global"));
    }

    if (current_bclist && !current_bclist->empty()) {
      const int search_start = max(0, segment_read_start - s.buffer_size);
      const int search_end = min((int)seq.length(), segment_read_end + s.buffer_size);
      const std::string search_region = seq.substr(search_start, search_end - search_start);

      unsigned int best_edit_distance = 100;
      unsigned int best_end_distance = 0;
      std::string best_barcode_match = "";
      int best_wiggle_offset = std::numeric_limits<int>::max(); // Initialize with a large value
      bool current_segment_unambiguous = false;

      for (const auto& known_bc : *current_bclist) {
        unsigned int editDistance, endDistance;
        editDistance = edit_distance(search_region, known_bc, endDistance, s.max_edit_distance);

        if (editDistance <= s.max_edit_distance) { // Only consider valid matches
            // Calculate the start of the barcode in the original sequence
            int current_match_start_in_seq = search_start + (endDistance - known_bc.length());
            // Calculate the wiggle offset from the ideal segment_read_start
            int current_wiggle_offset = std::abs(current_match_start_in_seq - segment_read_start);

            if (editDistance < best_edit_distance) {
                current_segment_unambiguous = true;
                best_edit_distance = editDistance;
                best_barcode_match = known_bc;
                best_end_distance = endDistance;
                best_wiggle_offset = current_wiggle_offset;
                if (editDistance == 0 && current_wiggle_offset == 0) break; // Exact match, no wiggle
            } else if (editDistance == best_edit_distance) {
                if (current_wiggle_offset < best_wiggle_offset) {
                    current_segment_unambiguous = true; // Found a better wiggle
                    best_edit_distance = editDistance; // Keep same edit distance
                    best_barcode_match = known_bc;
                    best_end_distance = endDistance;
                    best_wiggle_offset = current_wiggle_offset;
                } else if (current_wiggle_offset == best_wiggle_offset) {
                    current_segment_unambiguous = false; // Ambiguous if same edit distance and wiggle
                }
                // If current_wiggle_offset > best_wiggle_offset, do nothing, keep existing best
            }
        }
      }
      if (best_edit_distance <= s.max_edit_distance && current_segment_unambiguous) {
        barcode_result.features[s.name] = best_barcode_match;
        refined_segment_ends[i] = search_start + best_end_distance - 1;
        refined_segment_starts[i] = refined_segment_ends[i] - best_barcode_match.length() + 1;
      } else {
        barcode_result.found_all_matched_segments = false;
      }
    } else {
      // Discovery Mode (Fallback when no list provided)
      const int len = segment_read_end - segment_read_start;
      if (len > 0 && segment_read_start + len <= (int)seq.length()) {
        barcode_result.features[s.name] = seq.substr(segment_read_start, len);
        refined_segment_starts[i] = segment_read_start;
        refined_segment_ends[i] = segment_read_end - 1;
      } else {
        barcode_result.features[s.name] = "";
        refined_segment_starts[i] = segment_read_start;
        refined_segment_ends[i] = segment_read_start - 1;
      }
    }
  }
}

// Helper: Process MATCHED_SPLIT segments (Pass 2)
void refine_split_segments(const string &seq,
                           const std::vector<Segment> &segments,
                           const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
                           const std::unordered_map<std::string, BarcodeGroup> *group_map,
                           const std::vector<int> &read_to_segment_starts,
                           Barcode &barcode_result,
                           std::vector<int> &refined_segment_starts,
                           std::vector<int> &refined_segment_ends) {

  for (size_t i = 0; i < segments.size(); ++i) {
    if (segments[i].type != MATCHED_SPLIT) continue;
    if (refined_segment_starts[i] != -1) continue; // Already processed as part of a group

    const std::string group_name = segments[i].bc_list_name;
    if (group_map->find(group_name) == group_map->end()) {
      cerr << "Error: Undefined group " << group_name << ".\n"; exit(1);
    }

    const BarcodeGroup& bg = group_map->at(group_name);
    const std::vector<size_t>& split_group_indices = bg.segment_indices;

    std::vector<std::string> part_seqs;
    std::vector<int> part_starts;
    part_seqs.reserve(split_group_indices.size());
    part_starts.reserve(split_group_indices.size());
    bool possible_to_extract = true;

    for (size_t idx : split_group_indices) {
      const Segment& s = segments[idx];
      const int segment_read_start = read_to_segment_starts[idx];
      const int segment_read_end =
        (idx + 1 < read_to_segment_starts.size())
        ? read_to_segment_starts[idx + 1]
        : barcode_result.flank_end;

      const int search_start = max(0, segment_read_start - s.buffer_size);
      const int search_end = min((int)seq.length(), segment_read_end + s.buffer_size);

      if (search_end <= search_start) {
        possible_to_extract = false;
        break;
      }

      part_seqs.push_back(seq.substr(search_start, search_end - search_start));
      part_starts.push_back(search_start);
    }

    if (!possible_to_extract) {
      barcode_result.found_all_matched_segments = false;
      continue;
    }

    // Lookup group barcode list. If not present or empty, fall back to discovery mode.
    const unordered_set<string>* current_bclist = nullptr;
    auto it = known_barcodes_map->find(group_name);
    if (it != known_barcodes_map->end()) current_bclist = &(it->second);

    if (current_bclist && !current_bclist->empty()) {
      unsigned int best_edit_distance = 100;
      std::string best_barcode_match = "";
      std::vector<int> best_part_starts_in_read(split_group_indices.size(), 0);
      unsigned int best_total_wiggle_offset = std::numeric_limits<unsigned int>::max(); // Initialize with a large value
      bool current_group_unambiguous = false;

      std::vector<size_t> known_split_offsets;
      known_split_offsets.reserve(split_group_indices.size());
      size_t running_offset = 0;
      for (size_t idx : split_group_indices) {
        known_split_offsets.push_back(running_offset);
        running_offset += segments[idx].pattern.length();
      }

      for (const auto& known_bc : *current_bclist) {
        unsigned int total_edit_distance = 0;
        std::vector<int> current_part_starts;
        current_part_starts.reserve(split_group_indices.size());
        bool possible_match = true;

        for (size_t k = 0; k < split_group_indices.size(); ++k) {
          const size_t idx = split_group_indices[k];
          const int offset = known_split_offsets[k];
          int len = segments[idx].pattern.length();

          if (offset + len > (int)known_bc.length()) len = (int)known_bc.length() - offset;
          if (len <= 0) {
            possible_match = false;
            break;
          }

          const std::string known_part = known_bc.substr(offset, len);
          unsigned int endDist = 0;
          const unsigned int part_ed =
            edit_distance(part_seqs[k], known_part, endDist,
                          segments[idx].max_edit_distance);

          if (part_ed > (unsigned)segments[idx].max_edit_distance) {
            possible_match = false;
            break;
          }
          total_edit_distance += part_ed;
          current_part_starts.push_back(part_starts[k] + endDist - known_part.length());
        }

        if (!possible_match)
          continue;

        // Calculate current_total_wiggle_offset for the grouped barcode
        unsigned int current_total_wiggle_offset = 0;
        for (size_t k = 0; k < split_group_indices.size(); ++k) {
            const size_t idx = split_group_indices[k];
            current_total_wiggle_offset += std::abs(current_part_starts[k] - read_to_segment_starts[idx]);
        }

        if (total_edit_distance < best_edit_distance) {
          current_group_unambiguous = true;
          best_edit_distance = total_edit_distance;
          best_barcode_match = known_bc;
          best_part_starts_in_read = current_part_starts;
          best_total_wiggle_offset = current_total_wiggle_offset;
        } else if (total_edit_distance == best_edit_distance) {
            if (current_total_wiggle_offset < best_total_wiggle_offset) {
                current_group_unambiguous = true;
                best_edit_distance = total_edit_distance;
                best_barcode_match = known_bc;
                best_part_starts_in_read = current_part_starts;
                best_total_wiggle_offset = current_total_wiggle_offset;
            } else if (current_total_wiggle_offset == best_total_wiggle_offset) {
                current_group_unambiguous = false;
            }
        }
      }

      if (best_edit_distance <= (unsigned)bg.max_edit_distance && current_group_unambiguous) {
        size_t running_offset = 0;
        for (size_t k = 0; k < split_group_indices.size(); ++k) {
          const size_t idx = split_group_indices[k];
          int len = segments[idx].pattern.length();
          if (running_offset + len > best_barcode_match.length())
            len = best_barcode_match.length() - running_offset;

          barcode_result.features[segments[idx].name] =
            best_barcode_match.substr(running_offset, len);
          refined_segment_starts[idx] = best_part_starts_in_read[k];
          refined_segment_ends[idx] = best_part_starts_in_read[k] + len - 1;
          running_offset += segments[idx].pattern.length();
        }
        barcode_result.features[group_name] = best_barcode_match;
      } else {
        barcode_result.found_all_matched_segments = false;
      }

    } else {
      // Discovery mode for MATCHED_SPLIT:
      // No known list for the group, so extract each part based on the approximate
      // alignment-derived boundaries and concatenate to a group-level feature.
      std::string group_concat;
      group_concat.reserve(256);

      for (size_t k = 0; k < split_group_indices.size(); ++k) {
        const size_t idx = split_group_indices[k];
        const Segment& s = segments[idx];

        int segment_read_start = read_to_segment_starts[idx];
        int segment_read_end =
          (idx + 1 < read_to_segment_starts.size())
          ? read_to_segment_starts[idx + 1]
          : barcode_result.flank_end;

        int extract_start = max(0, segment_read_start);
        int extract_end = min((int)seq.length(), segment_read_end);

        // Anchor left boundary to previous refined split-part if we have it.
        if (k > 0) {
          const size_t prev_idx = split_group_indices[k - 1];
          if (refined_segment_ends[prev_idx] != -1) {
            extract_start = max(extract_start, refined_segment_ends[prev_idx] + 1);
          }
        }
        // Anchor right boundary to next refined split-part if we have it.
        if (k + 1 < split_group_indices.size()) {
          const size_t next_idx = split_group_indices[k + 1];
          if (refined_segment_starts[next_idx] != -1) {
            extract_end = min(extract_end, refined_segment_starts[next_idx]);
          }
        }

        // If still unresolved, fall back to expected pattern length, starting at the approximate position.
        if (extract_end <= extract_start) {
          extract_start = max(0, segment_read_start);
          extract_end = min((int)seq.length(), extract_start + (int)s.pattern.length());
        }

        if (extract_end < extract_start) extract_end = extract_start;

        const std::string part = seq.substr(extract_start, extract_end - extract_start);
        barcode_result.features[s.name] = part;
        refined_segment_starts[idx] = extract_start;
        refined_segment_ends[idx] = extract_start + (int)part.length() - 1;

        group_concat += part;
      }

      barcode_result.features[group_name] = group_concat;
    }
  }
}

// Helper: Process RANDOM segments (Pass 3)
void extract_random_segments(const string &seq,
                             const std::vector<Segment> &segments,
                             const std::vector<int> &read_to_segment_starts,
                             const std::vector<int> &refined_segment_starts,
                             const std::vector<int> &refined_segment_ends,
                             Barcode &barcode_result) {

  for (size_t i = 0; i < segments.size(); ++i) {
    const Segment &s = segments[i];
    if (s.type == MATCHED || s.type == MATCHED_SPLIT) continue;

    int extract_start = read_to_segment_starts[i];
    int extract_end = (i + 1 < read_to_segment_starts.size())
      ? read_to_segment_starts[i + 1]
      : barcode_result.flank_end;

    if (s.type == RANDOM) {
      // Anchor to previous matched segment if available
      if (i > 0 &&
        (segments[i - 1].type == MATCHED || segments[i - 1].type == MATCHED_SPLIT) &&
          refined_segment_ends[i - 1] != -1) {
        extract_start = refined_segment_ends[i - 1] + 1;
        extract_end = extract_start + (int)s.pattern.length();
      }
      // Anchor to next matched segment if available
      else if (i + 1 < segments.size() &&
        (segments[i + 1].type == MATCHED || segments[i + 1].type == MATCHED_SPLIT) &&
          refined_segment_starts[i + 1] != -1) {
        extract_end = refined_segment_starts[i + 1];
        extract_start = extract_end - (int)s.pattern.length();
      }

      int actual_extract_start = max(0, extract_start);
      int actual_extract_end = min((int)seq.length(), extract_end);

      if (actual_extract_end < actual_extract_start) actual_extract_end = actual_extract_start;

      barcode_result.features[s.name] =
        seq.substr(actual_extract_start, actual_extract_end - actual_extract_start);
    }
  }
}

vector<Barcode> big_barcode_search(
  const string &sequence,
  const std::unordered_map<std::string, std::unordered_set<std::string>> *known_barcodes_map,
  const std::unordered_map<std::string, BarcodeGroup> *group_map,
  int global_flank_max_editd,
  int global_barcode_edit_distance,
  const std::vector<Segment> &segments) {

  vector<Barcode> return_vec;
  string masked_sequence = sequence; // Work on a copy

  while (true) {
    Barcode result = extract_features(masked_sequence, segments, known_barcodes_map,
                                      group_map, global_flank_max_editd,
                                      global_barcode_edit_distance);

    if (result.flank_editd <= global_flank_max_editd && result.found_all_matched_segments) {
      return_vec.push_back(result);

      // Mask the found region to prevent re-finding it
      // result.flank_end is inclusive, so length is end - start + 1
      const int match_length = result.flank_end - result.flank_start + 1;

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
      out_stream << "@" << id << endl;
    else
      out_stream << ">" << id << endl;
    out_stream << read << "\n";
    if(is_fastq){
      out_stream << "+" << id << endl;
      out_stream << quals << "\n";
    }
}

std::string compose_new_id(
    const std::string &read_id, const Barcode &bc, int which, int total,
    bool chimeric, const std::vector<Segment> &segments,
    const std::unordered_map<std::string, BarcodeGroup> &group_map,
    std::string &primary_barcode) {

  std::ostringstream ss_suffix;
  ss_suffix << which << "of" << total;
  if (chimeric)
    ss_suffix << "_C";

  std::ostringstream prefix;

  // print the grouped barcodes as a whole first
  for (const auto &g : group_map) {
    auto it = bc.features.find(g.first);
    if (it != bc.features.end()) {
      prefix << it->second << "_";
      if (primary_barcode.empty())
        primary_barcode = it->second;
    } else {
      prefix << "NA_";
    }
  }

  for (const auto &s : segments) {
    if (s.type == FIXED)
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      prefix << it->second << "_";
      if (s.type == MATCHED && primary_barcode.empty())
        primary_barcode = it->second;
    } else {
      prefix << "NA_";
    }
  }

  std::string id_prefix = prefix.str();
  if (!id_prefix.empty() && id_prefix.back() == '_')
    id_prefix.pop_back();

  // Add CB tag (assume always present)
  std::string new_id = id_prefix + "#" + read_id + ss_suffix.str() + "\t" + "CB:Z:";
  std::string delim = "";
  for (const auto &g : group_map) {
    auto it = bc.features.find(g.first);
    if (it != bc.features.end()) {
      new_id += delim + g.first + ":" + it->second;
      delim = ",";
    }
  }
  for (const auto &s : segments) {
    if (s.type != MATCHED) 
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      new_id += delim + s.name + (s.name == "" ? "" : ":") + it->second;
      delim = ",";
    }
  }

  // Add UB tag (if present, concatonate all UB segments if multiple present)
  std::string ub;
  for (const auto &s : segments) {
    if (s.type != RANDOM)
      continue;
    auto it = bc.features.find(s.name);
    if (it != bc.features.end()) {
      ub += it->second;
    }
  }
  if (!ub.empty()) {
    new_id += "\tUB:Z:" + ub;
  }

  return new_id;
}

//print fastq or fasta lines..
void print_read(const string& read_id, const string& read, 
                const string qual,
                const vector<Barcode> & vec_bc, const string& prefix,
                bool split, map<string, ofstream> & barcode_file_streams,
                mutex & cout_mutex, mutex & file_mutex,
                bool trim_barcodes,
                bool chimeric,
                const std::vector<Segment> &segments,
                const std::unordered_map<std::string, BarcodeGroup> &group_map){

  const size_t vec_size = vec_bc.size();

  for (int b = 0; b < vec_size; b++) {
    const Barcode &bc = vec_bc.at(b);

    if (bc.flank_end < 0) {
      continue;
    }

    string primary_barcode_for_filename;
    string new_read_id =
      compose_new_id(read_id, bc, b + 1, vec_size, chimeric, segments, group_map, primary_barcode_for_filename);

    int read_start = bc.flank_end + 1;
    // work out the start and end base in case multiple barcodes
    int next_barcode_start = read.length();
    for (int k = 0; k < vec_size; k++) {
      if (b == k) continue;
      if (vec_bc.at(k).flank_start >= read_start) {
        next_barcode_start = std::min(next_barcode_start, vec_bc.at(k).flank_start);
      }
    }
    int read_length = next_barcode_start - read_start;
    read_length = std::max(0, read_length);

    string qual_new;
    if (!qual.empty()) {
      if (read_start + read_length > (int)qual.length()){
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

    if (read_new.empty()) {
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



// check if file already exists
bool file_exists(const std::string &filename) {
  std::ifstream infile(filename);
  return infile.good();
}

// Helper function to parse sub-options
void parse_sub_options(const string& opt_string, string& primary, map<string, string>& options) {
    size_t comma_pos = opt_string.find(',');
    if (comma_pos == string::npos) {
        primary = opt_string;
        return;
    }
    primary = opt_string.substr(0, comma_pos);
    
    size_t start = comma_pos + 1;
    while (start < opt_string.length()) {
        comma_pos = opt_string.find(',', start);
        string kv_pair;
        if (comma_pos == string::npos) {
            kv_pair = opt_string.substr(start);
            start = opt_string.length();
        } else {
            kv_pair = opt_string.substr(start, comma_pos - start);
            start = comma_pos + 1;
        }

        size_t colon_pos = kv_pair.find(':');
        if (colon_pos != string::npos) {
            options[kv_pair.substr(0, colon_pos)] = kv_pair.substr(colon_pos + 1);
        }
    }
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
  // Flag to track if any -b segment specified its own barcode list

  // Map to store group settings
  std::unordered_map<std::string, BarcodeGroup> group_map;
  
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
                     "d:k:i:b:u:x:e:f:n:s:hp:c:g:")) != EOF) {
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
      // If -g is used or -B with list is used, -k might be ambiguous or disallowed.
      // Allow -k only if NO groups or explicit lists are ever defined to keep it simple.
      
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
      // Move to avoid copying potentially large sets.
      known_barcodes_map["global"] = std::move(global_bclist);
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
      s.name = ""; 
      s.buffer_size = 0; 
      
      map<string, string> options;
      parse_sub_options(optarg, s.pattern, options);
      cerr << "Setting UMI to search for: " << s.pattern;

      if (options.count("name")) {
         s.name = options["name"];
         cerr << " (name: " << s.name << ")";
      }
      if (options.count("buffer")) {
         s.buffer_size = stoi(options["buffer"]);
         cerr << " (buffer: " << s.buffer_size << ")";
      }
      
      if (s.name.empty()) {
          umi_count++;
          s.name = (umi_count == 1) ? "UB" : "UB" + to_string(umi_count);
      } else {
          umi_count++; 
      }
      
      search_pattern.push_back(s);
      cerr << " -> named: " << s.name << "\n";
      params += 2;
      break;
    }
    case 'b': {
      Segment s;
      s.type = MATCHED; 
      s.name = ""; 
      s.buffer_size = 5; 
      s.max_edit_distance = edit_distance; 
      s.bc_list_name = "global"; 
      
      string group_name = "";
      map<string, string> options;
      parse_sub_options(optarg, s.pattern, options);
      
      cerr << "Adding barcode segment: " << s.pattern;

      if (options.count("name")) {
          s.name = options["name"];
          cerr << " (name: " << s.name << ")";
      }
      if (options.count("buffer")) {
          s.buffer_size = stoi(options["buffer"]);
          cerr << " (buffer: " << s.buffer_size << ")";
      }
      if (options.count("max_ed")) {
          s.max_edit_distance = stoi(options["max_ed"]);
          cerr << " (max_ed: " << s.max_edit_distance << ")";
      }
      if (options.count("group")) {
          group_name = options["group"];
          cerr << " (group: " << group_name << ")";
      }
      if (options.count("list")) {
          string file_name = options["list"];
          // Loading list logic
          if (known_barcodes_map.find(file_name) == known_barcodes_map.end()) {
             unordered_set<string> bclist;
             ifstream bc_file(file_name);
             if (bc_file.good()) {
                  string bc_line, bc;
                  while (getline(bc_file, bc_line)) {
                     istringstream line_stream(bc_line);
                     line_stream >> bc;
                     bclist.insert(bc);
                  }
                  known_barcodes_map[file_name] = std::move(bclist);
                  cerr << " (loaded " << known_barcodes_map[file_name].size() << " barcodes from " << file_name << ")";
             } else {
                 cerr << "\nError: Could not open barcode list file: " << file_name << "\n";
                 exit(1);
             }
          } else {
              cerr << " (using loaded list " << file_name << ")";
          }
          s.bc_list_name = file_name;
      }

      if (s.name.empty()) {
          barcode_count++;
          s.name = (barcode_count == 1) ? "" : "BC" + to_string(barcode_count);
      } else {
         barcode_count++; 
      }
      
      if (!group_name.empty()) {
          if (options.count("list")) {
               cerr << "\nError: Please specify the list file in the corresponding -g group definition if the segment is part of a group.\n";
               exit(1);
          }
          s.bc_list_name = group_name;
          s.type = MATCHED_SPLIT;
      }
      
      search_pattern.push_back(s);
      cerr << " -> named: " << s.name << "\n";
      params += 2;
      break;
    }
    case 'g': { // Define a group
        BarcodeGroup bg;
        bg.max_edit_distance = 2; 

        string group_name_dummy;
        map<string, string> options;
        parse_sub_options(optarg, bg.name, options); // The primary value is the group name

        cerr << "Defining group: " << bg.name;
        
        if (options.count("list")) {
            string list_file = options["list"];
            cerr << " (list: " << list_file << ")";
            
            // Re-using logic to load list, but here we store it under bg.name
            if (known_barcodes_map.find(bg.name) != known_barcodes_map.end()) {
                 cerr << "\nError: Group name " << bg.name << " already defined or used.\n";
                 exit(1);
            }
            
            unordered_set<string> bclist;
            ifstream bc_file(list_file);
            if (bc_file.good()) {
                 string bc_line, bc;
                 while (getline(bc_file, bc_line)) {
                    istringstream line_stream(bc_line);
                    line_stream >> bc;
                    bclist.insert(bc);
                 }
                 known_barcodes_map[bg.name] = std::move(bclist);
                 cerr << " (loaded " << known_barcodes_map[bg.name].size() << " barcodes)";
            } else {
                cerr << "\nError: Could not open barcode list file: " << list_file << "\n";
                exit(1);
            }
        } else {
            cerr << "\nError: Group definition must include 'list:filename'.\n";
            exit(1);
        }

        if (options.count("max_ed")) {
            bg.max_edit_distance = stoi(options["max_ed"]);
            cerr << " (max_ed: " << bg.max_edit_distance << ")";
        }
        
        group_map[bg.name] = bg;
        cerr << "\n";
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
  
  // Validation and sanity checks
  bool uses_global_list = false;
  bool uses_explicit_list = false;
  
  if (known_barcodes_map.find("global") != known_barcodes_map.end()) {
      uses_global_list = true;
  }
  
  // Check what search patterns rely on
  for (const auto& s : search_pattern) {
      if (s.type == MATCHED || s.type == MATCHED_SPLIT) {
          if (s.bc_list_name == "global") {
               if (!uses_global_list) {
                   // This segment needs global list but it's not provided
                   // TODO: What to do with multiple -b segments and descovery mode?
               }
          } else {
              uses_explicit_list = true;
          }
      }
  }

  if (uses_global_list && uses_explicit_list) {
      cerr << "Error: You cannot use -k (global list) when -b or -g specifies explicit barcode lists.\n";
      print_usage();
      exit(1);
  }

  // Populate group segment indices
  for (size_t i = 0; i < search_pattern.size(); ++i) {
      if (search_pattern[i].type == MATCHED_SPLIT) {
          string grp_name = search_pattern[i].bc_list_name;
          if (group_map.find(grp_name) != group_map.end()) {
              group_map[grp_name].segment_indices.push_back(i);
          } else {
              cerr << "Error: Barcode segment " << search_pattern[i].name << " refers to undefined group '" << grp_name << "'.\n";
              exit(1);
          }
      }
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
    cerr << "Reading reads from file " << reads_file << endl;
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
            out_stat_file << "\t" << (s.name.empty() ? "CB" : s.name);
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

  // --- Batch Threading ---
  // Read input in fixed-size chunks, distribute chunks to threads,
  // join threads, then write results in the same order as input.

  const int buffer_size = 2000; // reads per thread batch

  // Mutexes for synchronized output (only needed for split files)
  mutex cout_mutex;
  mutex file_mutex;
  map<string, ofstream> barcode_file_streams;

  while (*in) {
    std::vector<std::vector<SearchResult>> sr_v(n_threads);
    for (int t = 0; t < n_threads; ++t) {
      sr_v[t] = std::vector<SearchResult>();
      sr_v[t].reserve(buffer_size);
    }

    // Fill per-thread batches in order
    bool any_reads = false;
    for (int t = 0; t < n_threads; ++t) {
      for (int b = 0; b < buffer_size; ++b) {
        if (!*in) break;
        if (read_id_line.empty()) break;
        if (read_id_line[0] != '@' && read_id_line[0] != '>') {
          // tolerate trailing blank lines
          break;
        }

        SearchResult sr;
        istringstream line_stream(read_id_line);
        line_stream >> sr.read_id;
        if (!sr.read_id.empty()) sr.read_id.erase(0, 1);

        if (!is_fastq) { // fasta (account for multi-lines per read)
          string buffer_string;
          while (getline(*in, buffer_string) && !buffer_string.empty() && buffer_string[0] != '>') {
            sr.line += buffer_string;
          }
          read_id_line = buffer_string;
        } else { // fastq
          if (!getline(*in, sr.line)) {
            read_id_line.clear();
            break;
          }
          getline(*in, line); // skip '+' line
          getline(*in, sr.qual_scores);
          getline(*in, read_id_line);
        }

        sr_v[t].push_back(std::move(sr));
        any_reads = true;
      }

      // Stop allocating further thread batches if input exhausted
      if (!*in || read_id_line.empty()) {
        for (int t2 = t + 1; t2 < n_threads; ++t2) sr_v[t2].clear();
        break;
      }
    }

    if (!any_reads) break;

    // Launch threads to process each batch
    std::vector<std::thread> threads;
    threads.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
      if (sr_v[t].empty()) continue;
      threads.emplace_back(search_read_batch, ref(sr_v[t]), ref(known_barcodes_map), ref(group_map),
                           flank_edit_distance, edit_distance, ref(search_pattern));
    }

    // Join all threads
    for (auto &th : threads) th.join();

    // Write results in deterministic order: t then r (matches fill order above)
    for (int t = 0; t < (int)sr_v.size(); ++t) {
      for (int r = 0; r < (int)sr_v[t].size(); ++r) {
        SearchResult &sr = sr_v[t][r];

        r_count++;
        if (sr.count > 0) bc_count++;
        if (sr.chimeric) multi_bc_count++;

        // Count barcodes from features map (first MATCHED segment)
        for (const auto& bc : sr.vec_bc_for) {
          for (const auto& s : search_pattern) {
            if (s.type == MATCHED) {
              auto it = bc.features.find(s.name);
              if (it != bc.features.end()) {
                barcode_counts[it->second]++;
                break;
              }
            }
          }
        }
        for (const auto& bc : sr.vec_bc_rev) {
          for (const auto& s : search_pattern) {
            if (s.type == MATCHED) {
              auto it = bc.features.find(s.name);
              if (it != bc.features.end()) {
                barcode_counts[it->second]++;
                break;
              }
            }
          }
        }

        if (!known_barcodes_map.empty()) {
          print_stats(sr.read_id, sr.vec_bc_for, out_stat_file, search_pattern);
          print_stats(sr.read_id, sr.vec_bc_rev, out_stat_file, search_pattern);

          print_read(sr.read_id + "_+", sr.line, sr.qual_scores, sr.vec_bc_for, out_filename_prefix,
                     split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex,
                     remove_barcodes, print_chimeric && sr.chimeric, search_pattern, group_map);

          if (remove_barcodes || sr.vec_bc_for.empty()) {
            string rev_qual = sr.qual_scores;
            reverse(rev_qual.begin(), rev_qual.end());
            print_read(sr.read_id + "_-", sr.rev_line, rev_qual, sr.vec_bc_rev, out_filename_prefix,
                       split_file_by_barcode, barcode_file_streams, cout_mutex, file_mutex,
                       remove_barcodes, print_chimeric && sr.chimeric, search_pattern, group_map);
          }
        }

        if ((r_count < 100000 && r_count % 10000 == 0) || (r_count % 100000 == 0)) {
          cerr << r_count / 1000000.0 << " million reads processed.." << endl;
        }
      }
    }
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
