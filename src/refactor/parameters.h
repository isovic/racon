/*
 * parameters.h
 *
 *  Created on: January 18, 2017
 *      Author: Ivan Sovic
 */

#ifndef SRC_RACON_PARAMETERR_H_
#define SRC_RACON_PARAMETERR_H_

#include <stdint.h>
#include <memory>
#include <string>
#include "argparser.h"
#include "types.h"

namespace is {

class Parameters;

std::unique_ptr<Parameters> createParameters(int argc, char* argv[]);

class Parameters {
 public:
  friend std::unique_ptr<Parameters> createParameters(int argc, char* argv[]);

  const std::string& aln_path() const {
    return aln_path_;
  }

  int32_t aln_type() const {
    return aln_type_;
  }

  int64_t batch_of_windows() const {
    return batch_of_windows_;
  }

  const std::string& consensus_path() const {
    return consensus_path_;
  }

  const std::string& contigs_path() const {
    return contigs_path_;
  }

  bool do_erc() const {
    return do_erc_;
  }

  double error_rate() const {
    return error_rate_;
  }

  int32_t gap_ext() const {
    return gap_ext_;
  }

  int32_t gap_open() const {
    return gap_open_;
  }

  int32_t match() const {
    return match_;
  }

  int32_t mismatch() const {
    return mismatch_;
  }

  int64_t num_batches() const {
    return num_batches_;
  }

  int32_t num_threads() const {
    return num_threads_;
  }

  const std::string& out_fmt() const {
    return out_fmt_;
  }

  const OverlapFormat& overlap_format() const {
    return overlap_format_;
  }

  double qv_threshold() const {
    return qv_threshold_;
  }

  const std::string& reads_path() const {
    return reads_path_;
  }

  int64_t start_window() const {
    return start_window_;
  }

  bool use_contig_qvs() const {
    return use_contig_qvs_;
  }

  int32_t verbose_level() const {
    return verbose_level_;
  }

  double win_ovl_margin() const {
    return win_ovl_margin_;
  }

  int64_t window_len() const {
    return window_len_;
  }

 private:
  Parameters(int argc, char* argv[]);

  // Input output options.
  std::string contigs_path_ = "";
  std::string aln_path_ = "";
  std::string consensus_path_ = "";
  std::string out_fmt_ = "fasta";
  std::string reads_path_ = "";

  // Input format type. Default is PAF.
  OverlapFormat overlap_format_;
  // This should be deprecated.
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
