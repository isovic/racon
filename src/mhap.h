/*
 * mhap.h
 *
 *  Created on: May 29, 2016
 *      Author: isovic
 */

#ifndef SRC_MHAP_H_
#define SRC_MHAP_H_

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <stdio.h>

#include "log_system/log_system.h"
#include "sequences/sequence_file.h"

class MHAPLine {
 public:
  MHAPLine() :
    Aid(0), Bid(0), perc_err(0.0), shared_minmers(0), Arev(0), Astart(0), Aend(0), Alen(0), Brev(0), Bstart(0), Bend(0), Blen(0) {
  }

  int Parse(std::string &line) {
    std::istringstream iss(line);
    if (!(iss >> Aid >> Bid >> perc_err >> shared_minmers >> Arev >> Astart >> Aend >> Alen >> Brev >> Bstart >> Bend >> Blen)) {
      return 1;
    }
    return 0;
  }

  std::string Verbose() const {
    std::stringstream ss;
    ss << Aid << " " << Bid << " " << perc_err << " " << shared_minmers << " " << Arev << " " << Astart << " " << Aend << " " << Alen << " " << Brev << " " << Bstart << " " << Bend << " " << Blen;
    return ss.str();
  }

  int CheckConstraints(double max_dist_ratio) const {
    double Adist = Aend - Astart;
    double Bdist = Bend - Bstart;
    double ratio = (Adist > Bdist) ? (1.0f - Bdist / Adist) : (1.0f - Adist/Bdist);
    if (ratio > max_dist_ratio) { return 1; }
    return 0;
  }

  int64_t Aid, Bid;
  double perc_err;
  int64_t shared_minmers;
  int64_t Arev, Astart, Aend, Alen;     // start is zero-based, end points to a position right after the last inclusive base.
  int64_t Brev, Bstart, Bend, Blen;
};

int ParseMHAP(const std::string &mhap_path, std::vector<MHAPLine> &ret_overlaps);
int FilterMHAP(const std::vector<MHAPLine> &overlaps_in, std::vector<MHAPLine> &overlaps_out, float error_rate);
int FilterMHAPErc(const std::vector<MHAPLine> &overlaps_in, std::vector<MHAPLine> &overlaps_out, float error_rate);
int EdlibNWWrapper(const int8_t *read_data, int64_t read_length,
                   const int8_t *reference_data, int64_t reference_length,
                   int64_t* ret_alignment_position_start, int64_t *ret_alignment_position_end,
                   int64_t *ret_edit_distance, std::vector<unsigned char> &ret_alignment);
int AlignMHAP(const SequenceFile &refs, const SequenceFile &reads, const std::vector<MHAPLine> &overlaps, int32_t num_threads, SequenceFile &aligned);

#endif /* SRC_MHAP_H_ */
