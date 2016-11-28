/*
 * extract.cc
 *
 *  Created on: Nov 27, 2016
 *      Author: isovic
 */

#include "consensus/consensus.h"
#include "log_system/log_system.h"
#include "utility/utility_general.h"
#include <stdint.h>
#include <algorithm>
#include <sstream>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>

void ExtractWindowFromSampledSeqs(const SingleSequence *contig, IntervalTreeSampled &interval_tree,
                                  int64_t window_start, int64_t window_end, double qv_threshold, bool use_contig_qvs,
                                  std::vector<std::string> &window_seqs, std::vector<std::string> &window_qv, std::vector<const SingleSequence *> &window_refs,
                                  std::vector<uint32_t> &window_starts, std::vector<uint32_t> &window_ends,
                                  std::vector<uint32_t> &starts_on_read, std::vector<uint32_t> &ends_on_read, FILE *fp_window) {
  if (window_start > window_end) {
    return;
  }

  int64_t temp_window_end = std::min((int64_t) window_end, (int64_t) (contig->get_sequence_length()-1));
  window_refs.push_back(contig);
  window_seqs.push_back(GetSubstring((char *) (contig->get_data() + window_start), (temp_window_end - window_start + 1)));
  if (use_contig_qvs == false || contig->get_quality() == NULL || contig->get_quality_length() == 0) {
    std::string dummy_quals((temp_window_end - window_start + 1), '!');
    window_qv.push_back(dummy_quals);
  } else {
    window_qv.push_back(GetSubstring((char *) (contig->get_quality() + window_start), (temp_window_end - window_start + 1)));
  }
  window_starts.push_back(0);
  window_ends.push_back(temp_window_end - window_start);
  starts_on_read.push_back(window_start);
  ends_on_read.push_back(window_end - 1);

  // Find seqs which fall into the window region.
  std::vector<IntervalSampled> intervals;
  interval_tree.findOverlapping(window_start, temp_window_end, intervals);

  // For each seq, extract its segment which falls into the window.
  for (int64_t i=0; i<intervals.size(); i++) {
    std::shared_ptr<SampledAlignment> sal = intervals[i].value;
    std::shared_ptr<SingleSequence> seq = sal->seq;

    int32_t start_seq = sal->find(window_start);
    int32_t end_seq = sal->find(temp_window_end);

    auto aln = sal->seq->get_aln();

    uint32_t seq_start_in_window = 0;
    uint32_t seq_end_in_window = temp_window_end - window_start;

    // If the window starts within the interval (after the window_start location).
    if (start_seq == -1 && intervals[i].start > window_start) {
      start_seq = intervals[i].value->mhap.Astart;
      seq_start_in_window = intervals[i].start - window_start;
      seq_start_in_window = std::max((uint32_t) 0, (uint32_t) ((int32_t) seq_start_in_window - 0));
    } else if (start_seq < 0) {
      assert(start_seq >= 0 && "ERROR: start_seq is < 0");
    }

    if (end_seq == -1 && intervals[i].stop < window_end) {
      end_seq = intervals[i].value->mhap.Aend;
      seq_end_in_window = intervals[i].stop - window_start;
      seq_end_in_window = std::min((uint32_t) (temp_window_end - window_start), (uint32_t) ((int32_t) seq_end_in_window + 0));
    } else if (start_seq < 0) {
      assert(end_seq >= 0 && ("ERROR: start_seq is < 0"));
    }

    std::string seq_data = GetSubstring((char *) (seq->get_data() + start_seq), end_seq - start_seq + 1);
    std::string seq_qual = (seq->get_quality() != NULL) ? (GetSubstring((char *) (seq->get_quality() + start_seq), end_seq - start_seq + 1)) : (std::string((end_seq - start_seq + 1), '!' + 0));

    // Safety percaution.
    if (seq_data.size() < 2) { continue; }

    double avg_qual;
    for (int64_t j=0; j<seq_qual.size(); j++) {
      avg_qual += (double) (seq_qual[j] - '!');
    }
    avg_qual /= std::max((double) seq_qual.size(), 1.0);

    if (avg_qual >= qv_threshold) {
      window_refs.push_back(&(*seq));
      window_seqs.push_back(seq_data);
      window_starts.push_back(seq_start_in_window);
      window_ends.push_back(seq_end_in_window);
      window_qv.push_back(seq_qual);

      starts_on_read.push_back(start_seq);
      ends_on_read.push_back(end_seq);
    }

    if (fp_window) {
      #ifndef WINDOW_OUTPUT_IN_FASTQ
        fprintf (fp_window, ">%s Window_%d_to_%d\n%s\n", seq->get_header(), window_start, temp_window_end, window_seqs.back().c_str());
      #else
        fprintf (fp_window, "@%s Window_%d_to_%d\n%s\n", seq->get_header(), window_start, temp_window_end, window_seqs.back().c_str());
        fprintf (fp_window, "+\n%s\n", window_qv.back().c_str());
      #endif
    }
  }
}




