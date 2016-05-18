/*
 * parameters.h
 *
 *  Created on: Feb 23, 2016
 *      Author: isovic
 */

#ifndef SRC_RACON_PARAMETERR_H_
#define SRC_RACON_PARAMETERR_H_

#include <string>


struct ProgramParameters {
  std::string raw_contigs_path = "";
  std::string aln_path = "";
  std::string temp_alt_contig_path = "";
  std::string consensus_path = "";
  std::string out_fmt = "fasta";
  int64_t window_len = 500;
  int64_t batch_of_windows = 20000;
  int64_t num_batches = -1;
  int64_t start_window = 0;
  double qv_threshold = 10.0;
  int32_t num_threads = 4;

  std::string program_bin;
  std::string program_folder;
  std::vector<std::string> cmd_arguments;

  int32_t match = 5;
  int32_t mismatch = -4;
  int32_t gap_open = -8;
  int32_t gap_ext = -6;
  int32_t aln_type = 1;     // SW 0, NW 1, OV 2

  double win_ovl_margin = 0.00; // 0.05;

  bool do_realign = false;
  std::string realigned_aln_path = "";
  bool do_pileup = false;

  int32_t verbose_level = 1;
};

#endif /* SRC_PARAMETERS_H_ */
