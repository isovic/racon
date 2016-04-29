/*
 * parameters.h
 *
 *  Created on: Feb 23, 2016
 *      Author: isovic
 */

#ifndef SRC_CONSISE_PARAMETERR_H_
#define SRC_CONSISE_PARAMETERR_H_

#include <string>


struct ProgramParameters {
  std::string raw_contigs_path = "";
  std::string aln_path = "";
  std::string temp_alt_contig_path = "";
  std::string consensus_path = "";
  std::string temp_window_path = "";
  std::string out_fmt = "fasta";
  int64_t window_len = 1000;
  int64_t batch_of_windows = 200;
  int64_t num_batches = -1;
  int64_t start_window = 0;
  double qv_threshold = 10.0;
  std::string msa = "mafft";
  std::string mafft_folder = "../tools/mafft-7.273-with-extensions/";
  std::string poav2_folder = "../tools/poaV2";
  int32_t num_threads = 4;

  std::string program_bin;
  std::string program_folder;
  std::vector<std::string> cmd_arguments;

  int32_t match = 1;
  int32_t mismatch = -1;
  int32_t gap_open = -1;
  int32_t gap_ext = -1;
  int32_t aln_type = 1;     // SW 0, NW 1, OV 2

  double win_ovl_margin = 0.00; // 0.05;

  std::string realigned_aln_path = "";
};

#endif /* SRC_PARAMETERS_H_ */
