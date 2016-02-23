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
  std::string alt_contig_path = "";
  std::string out_fmt = "fasta";
  int64_t window_len = 5000;
  double qv_threshold = 10.0;
};



#endif /* SRC_PARAMETERS_H_ */
