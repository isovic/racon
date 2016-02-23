/*
 * main.cc
 *
 *  Created on: Jan 24, 2016
 *      Author: isovic
 */

#include <stdio.h>
#include "log_system/log_system.h"
#include <sstream>
#include "sequences/sequence_file.h"
#include "sequences/sequence_alignment_test.h"
#include "consensus/consensus.h"
#include "argparser.h"
#include "parameters.h"

void RunTests() {
  TEST_CLASS_SEQUENCE_ALIGNMENT();
  exit(0);
}

//void PrintUsage() {
//  fprintf (stderr, "Expected 3 parameters!\n\n");
//  fprintf (stderr, "Usage:\n");
//  fprintf (stderr, "\t%s <layout.gfa> <alignments.sam> <out_alt_contigs.sam>\n", argv[0]);
//}

int main(int argc, char* argv[]) {
  //  RunTests();

  bool help = false;
  ProgramParameters parameters;
  ArgumentParser argparser;
  argparser.AddArgument(&(parameters.raw_contigs_path), VALUE_TYPE_STRING, "", "raw", "", "Path to the raw contig sequences (output from the layout step). GFA, FASTA/FASTQ or SAM formats allowed.", -3, "Input/Output options");
  argparser.AddArgument(&(parameters.aln_path), VALUE_TYPE_STRING, "", "aln", "", "Path to a SAM file with read-to-raw contig alignments.", -2, "Input/Output options");
  argparser.AddArgument(&(parameters.alt_contig_path), VALUE_TYPE_STRING, "", "temp", "", "Extracted alternate contigs. Output is in SAM format.", -1, "Input/Output options");
  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other options");

  argparser.ProcessArguments(argc, argv);

  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
//    PrintUsage();
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  std::string gfa = argv[1];
  SequenceFile seqs_gfa(SEQ_FORMAT_GFA, parameters.raw_contigs_path);
  seqs_gfa.Verbose(stdout);

  std::string sam = argv[2];
  SequenceFile seqs_sam(SEQ_FORMAT_SAM, parameters.aln_path);
  seqs_sam.Verbose(stdout);

  std::string alt_contig_path = argv[3];

  Consensus(seqs_gfa, seqs_sam, alt_contig_path);

	return 0;
}
