/*
 * parameters.h
 *
 *  Created on: Feb 23, 2016
 *      Author: isovic
 */

#ifndef SRC_RACON_PARAMETERR_H_
#define SRC_RACON_PARAMETERR_H_

#include <stdint.h>
#include <string>
#include "argparser.h"

namespace is {

class OverlapFormat {
 public:
  static OverlapFormat Paf() { return OverlapFormat(kPaf); }
  static OverlapFormat Mhap() { return OverlapFormat(kMhap); }
  static OverlapFormat Sam() { return OverlapFormat(kSam); }

  bool isPaf() { return x_ == kPaf; }
  bool isMhap() { return x_ == kMhap; }
  bool isSam() { return x_ == kSam; }

 private:
  explicit OverlapFormat(int32_t x): x_(x) { }
  int32_t x_;
  enum OvlFormatTypes {
    kPaf = 0,
    kMhap = 1,
    kSam = 2
  };
};

class Parameters {
 public:
  Parameters(int argc, char* argv[]);

 private:
  // Input output options.
  std::string raw_contigs_path_ = "";
  std::string aln_path_ = "";
  std::string consensus_path_ = "";
  std::string out_fmt_ = "fasta";
  std::string reads_path_ = "";

  // Input format type. Default is PAF.
  OverlapFormat overlap_format_;
  // This should be deprecated.
  bool is_paf_ = true;
  bool is_mhap_ = false;
  bool is_sam_ = false;

  // Algorithmic/debug options.
  int64_t window_len_ = 500;
  int64_t batch_of_windows_ = 20000;
  int64_t num_batches_ = -1;
  int64_t start_window_ = 0;
  double qv_threshold_ = 10.0;
  int32_t num_threads_ = 4;
  double win_ovl_margin_ = 0.00;     // Window overlapping is disabled by default.
  bool use_contig_qvs_ = false;      // If true and contigs are in FASTQ format, the QVs which are in the file will be used for window consensuses. If false, dummy QVs of 0 will be used instead.
  bool do_erc_ = false;

  double error_rate_ = 0.30;

  int32_t verbose_level_ = 1;

  // Information about the program being run.
  std::string program_bin_;
  std::string program_folder_;
  std::vector<std::string> cmd_arguments_;

  // Alignment options for SPOA.
  int32_t match_ = 5;
  int32_t mismatch_ = -4;
  int32_t gap_open_ = -8;
  int32_t gap_ext_ = -6;
  int32_t aln_type_ = 1;     // SW 0, NW 1, OV 2
};

}

#endif /* SRC_PARAMETERS_H_ */
