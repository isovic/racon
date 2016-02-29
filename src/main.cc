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
  argparser.AddArgument(&(parameters.consensus_path), VALUE_TYPE_STRING, "", "out", "", "Output consensus sequence.", -1, "Input/Output options");
  argparser.AddArgument(&(parameters.temp_window_path), VALUE_TYPE_STRING, "", "winpath", "temp.window.fasta", "A window of alternate contigs, used to feed MSA tools from disk.", 0, "Input/Output options");
//  argparser.AddArgument(&(parameters.temp_window_path), VALUE_TYPE_STRING, "", "winpath", "", "A window of alternate contigs.", 0, "Input/Output options");

  argparser.AddArgument(&(parameters.num_threads), VALUE_TYPE_INT32, "t", "threads", "4", "Number of threads to use.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.qv_threshold), VALUE_TYPE_DOUBLE, "", "bq", "10.0", "Threshold for the average base quality of the input reads. If a read has average BQ < specified, the read will be skipped.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.new_seq_percent), VALUE_TYPE_DOUBLE, "", "minnewseq", "0.80", "Minimum percentage of new sequence included in an alternate contig. E.g. 0.80 would mean there can be at least 80% of bases obtained from alternate reads aligned to the contig, and 20% used from the input raw draft contig.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.percent_overlap), VALUE_TYPE_DOUBLE, "", "maxovl", "0.01", "When generating alternate contigs, take alignments of sequences which overlap at most this much.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.window_len), VALUE_TYPE_INT64, "w", "winlen", "500", "Length of the window to perform POA on.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.batch_of_windows), VALUE_TYPE_INT64, "b", "winbatch", "200", "Size of the batch in which to process windows. After a batch is finished, consensus of the windows is joined and output to file.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.msa), VALUE_TYPE_STRING, "", "msa", "poa", "Choose the MSA algorithm to use. Available options: 'maff', 'poa', 'poav2'.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.mafft_folder), VALUE_TYPE_STRING, "", "mafft", "../tools/mafft-7.273-with-extensions/", "Relative path to the folder containing MAFFT. Relative to the folder this binary resides in.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.poav2_folder), VALUE_TYPE_STRING, "", "poav2", "../tools/poaV2/", "Relative path to the folder containing POAv2. Relative to the folder this binary resides in.", 0, "Algorithm");

  argparser.AddArgument(&(parameters.temp_alt_contig_path), VALUE_TYPE_STRING, "", "altctgs", "", "Extracted alternate contigs. Output is in SAM format.", 0, "Debug");

  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other options");

  if (argc == 1) {
    fprintf (stderr, "  %s [options] <raw> <aln> <temp>\n\n", argv[0]);
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  // Process the command line arguments.
  argparser.ProcessArguments(argc, argv);

  // Store the command line arguments for later use.
  for (int32_t i=0; i<argc; i++) { parameters.cmd_arguments.push_back(argv[i]); }
  parameters.program_folder = parameters.cmd_arguments[0].substr(0, parameters.cmd_arguments[0].find_last_of("\\/"));
  parameters.program_bin = parameters.cmd_arguments[0].substr(parameters.cmd_arguments[0].find_last_of("\\/") + 1);

  // Verbose the current state of the parameters after.
  fprintf (stderr, "%s\n\n", argparser.VerboseArguments().c_str());

  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
    fprintf (stderr, "  %s [options] <raw> <aln> <temp>\n\n", argv[0]);
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  std::string gfa = argv[1];
  SequenceFile seqs_gfa(SEQ_FORMAT_AUTO, parameters.raw_contigs_path);
  seqs_gfa.Verbose(stdout);

  std::string sam = argv[2];
  SequenceFile seqs_sam(SEQ_FORMAT_SAM, parameters.aln_path);
  seqs_sam.Verbose(stdout);

  std::string alt_contig_path = argv[3];

  Consensus(parameters, seqs_gfa, seqs_sam);

	return 0;
}
