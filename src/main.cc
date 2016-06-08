/*
 * main.cc
 *
 *  Created on: Jan 24, 2016
 *      Author: isovic
 */

#include <sequences/sequence_test.h>
#include <stdio.h>
#include "log_system/log_system.h"
#include <sstream>
#include "sequences/sequence_file.h"
#include "consensus/consensus.h"
#include "argparser.h"
#include "parameters.h"
#include "mhap.h"

void RunTests() {
  TEST_CLASS_SEQUENCE_ALIGNMENT();
  exit(0);
}

int main(int argc, char* argv[]) {
  //  RunTests();

  bool help = false;
  ProgramParameters parameters;
  ArgumentParser argparser;
  argparser.AddArgument(&(parameters.reads_path), VALUE_TYPE_STRING, "", "reads", "", "Reads in FASTQ format.", -4, "Input/Output options");
  argparser.AddArgument(&(parameters.aln_path), VALUE_TYPE_STRING, "", "alnpath", "", "Path to a MHAP file with read-to-target overlaps.", -3, "Input/Output options");
  argparser.AddArgument(&(parameters.raw_contigs_path), VALUE_TYPE_STRING, "", "raw", "", "Path to the raw contig/read sequences (output from the layout step). GFA, FASTA/FASTQ or SAM formats allowed.", -2, "Input/Output options");
  argparser.AddArgument(&(parameters.consensus_path), VALUE_TYPE_STRING, "", "out", "", "Output consensus sequence.", -1, "Input/Output options");

  argparser.AddArgument(&(parameters.is_sam), VALUE_TYPE_BOOL, "", "sam", "0", "SAM is provided instead of MHAP. The reads file will be ignored, and seq and qual fields from the SAM file will be used.", 0, "Input/Output options");
  argparser.AddArgument(&(parameters.is_paf), VALUE_TYPE_BOOL, "", "paf", "0", "Overlaps are in PAF format instead of MHAP.", 0, "Input/Output options");
//  argparser.AddArgument(&(parameters.is_mhap), VALUE_TYPE_BOOL, "", "mhap", "0", "Overlaps are in PAF format instead of MHAP.", 0, "Input/Output options");

  argparser.AddArgument(&(parameters.qv_threshold), VALUE_TYPE_DOUBLE, "", "bq", "10.0", "Threshold for the average base quality of the input reads. If a read has average BQ < specified, the read will be skipped. If value is < 0.0, filtering is disabled.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.window_len), VALUE_TYPE_INT64, "w", "winlen", "500", "Length of the window to perform POA on.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.do_pileup), VALUE_TYPE_BOOL, "", "pileup", "0", "Simple pileup + majority vote consensus will be performed instead of using Spoa. Much faster, but less accurate.", 0, "Algorithm");
  argparser.AddArgument(&(parameters.num_threads), VALUE_TYPE_INT32, "t", "threads", "4", "Number of threads to use.", 0, "Control");
  argparser.AddArgument(&(parameters.batch_of_windows), VALUE_TYPE_INT64, "b", "winbatch", "20000", "Size of the batch in which to process windows. After a batch is finished, consensus of the windows is joined and output to file.", 0, "Control");
  argparser.AddArgument(&(parameters.num_batches), VALUE_TYPE_INT64, "", "num-batches", "-1", "The number of batches which to process", 0, "Control");
  argparser.AddArgument(&(parameters.start_window), VALUE_TYPE_INT64, "", "start-window", "0", "ID of the window to start processing from.", 0, "Control");
  argparser.AddArgument(&(parameters.do_erc), VALUE_TYPE_BOOL, "", "erc", "0", "Perform error-correction instead of contig consensus. The only difference is in the type of parallelization to achieve better performance.", 0, "Control");
  argparser.AddArgument(&(parameters.error_rate), VALUE_TYPE_DOUBLE, "e", "error-rate", "0.30", "Maximum allowed error rate. Used for filtering faulty overlaps.", 0, "Algorithm");

  argparser.AddArgument(&(parameters.match), VALUE_TYPE_INT32, "M", "match", "5", "Match score (positive value).", 0, "Alignment");
  argparser.AddArgument(&(parameters.mismatch), VALUE_TYPE_INT32, "X", "mismatch", "-4", "Mismatch penalty (negative value expected).", 0, "Alignment");
  argparser.AddArgument(&(parameters.gap_open), VALUE_TYPE_INT32, "G", "gapopen", "-8", "Gap open penalty (negative value expected).", 0, "Alignment");
  argparser.AddArgument(&(parameters.gap_ext), VALUE_TYPE_INT32, "E", "gapext", "-6", "Gap extend penalty (negative value expected).", 0, "Alignment");
//  argparser.AddArgument(&(parameters.aln_type), VALUE_TYPE_INT32, "a", "align", "1", "Alignment algorithm: 0 for local, 1 for global and 2 for overlap.", 0, "Alignment");
//  parameters.aln_type = 0;
  argparser.AddArgument(&(parameters.verbose_level), VALUE_TYPE_INT32, "v", "verbose", "5", "Verbose level. 0 off, 1 low, 2 medium, 3 high, 4 and 5 all levels, 6-9 debug.", 0, "Other");

  // TODO: Deprecated feature. Consider removing permanently.
//  argparser.AddArgument(&(parameters.win_ovl_margin), VALUE_TYPE_DOUBLE, "", "ovl-margin", "0.0", "Fraction of the window size to overlap the windows by.", 0, "Algorithm");
//  argparser.AddArgument(&(parameters.temp_alt_contig_path), VALUE_TYPE_STRING, "", "altctgs", "", "Extracted alternate contigs. Output is in SAM format.", 0, "Debug");

  // TODO: These two options are currently in dev.
//  argparser.AddArgument(&(parameters.do_realign), VALUE_TYPE_BOOL, "", "realign", "0", "If enabled, input SAM file will be realigned. In this mode, BQ filtering cannot be used. Realigned SAM will be output to stdout unless --rsam is specified.", 0, "Control");
//  argparser.AddArgument(&(parameters.realigned_aln_path), VALUE_TYPE_STRING, "", "rsam", "", "If specified, the input SAM file will be realigned and output to the specified path.", 0, "Control");


  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other options");

  if (argc == 1) {
    fprintf (stderr, "  %s [options] <reads.fastq> <overlaps.mhap> <raw_contigs.fasta> <out_consensus.fasta>\n\n", argv[0]);
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    exit(1);
  }

  // Process the command line arguments.
  argparser.ProcessArguments(argc, argv);
  // Store the command line arguments for later use.
  for (int32_t i=0; i<argc; i++) { parameters.cmd_arguments.push_back(argv[i]); }
  parameters.program_folder = parameters.cmd_arguments[0].substr(0, parameters.cmd_arguments[0].find_last_of("\\/"));
  parameters.program_bin = parameters.cmd_arguments[0].substr(parameters.cmd_arguments[0].find_last_of("\\/") + 1);
  //  // Verbose the current state of the parameters after.
  //  fprintf (stderr, "%s\n\n", argparser.VerboseArguments().c_str());

  // Sanity check on parameter values.
  if (parameters.is_sam == true && parameters.is_paf == true) {
    fprintf (stderr, "ERROR: More than one input overlap/alignment format specified. Exiting.\n");
    exit(1);
  }
  if (parameters.is_sam == false && parameters.reads_path == "") {
    fprintf (stderr, "ERROR: Reads file not specified! Exiting.\n");
    exit(1);
  }



  if (parameters.verbose_level == 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
  } else if (parameters.verbose_level > 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
  }

  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(parameters.verbose_level);

  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
    fprintf (stderr, "  %s [options] <raw> <aln> <temp>\n\n", argv[0]);
    fprintf (stderr, "%s\n", argparser.VerboseUsage().c_str());
    fflush(stderr);
    exit(1);
  }

  std::string gfa = argv[1];
  SequenceFile seqs_gfa(SEQ_FORMAT_AUTO, parameters.raw_contigs_path);

  SequenceFile *seqs_sam = NULL;

  if (parameters.is_sam == true) {
    std::string sam = parameters.aln_path;
    LOG_ALL("Using SAM for input alignments. (%s)\n", sam.c_str());
    LOG_ALL("Parsing the SAM file.\n");
    seqs_sam = new SequenceFile(SEQ_FORMAT_SAM, parameters.aln_path);
  } else {
    std::string mhap = parameters.aln_path;
    LOG_ALL("Using MHAP for input alignments. (%s)\n", mhap.c_str());
    std::vector<MHAPLine> overlaps, overlaps_filtered, overlaps_final;

    LOG_ALL("Parsing the MHAP file.\n");
    ParseMHAP(mhap, overlaps);
    LOG_ALL("Filtering MHAP overlaps.\n");
    FilterMHAP(overlaps, overlaps_final, parameters.error_rate);

//    if (parameters.do_erc == false) {
//      FilterMHAP(overlaps, overlaps_final, parameters.error_rate);
//    } else {
//      FilterMHAPErc(overlaps, overlaps_filtered, parameters.error_rate);
//      DuplicateAndSwitch(overlaps_filtered, overlaps_final);
//    }

    LOG_ALL("Loading reads.\n");
    SequenceFile seqs_reads(SEQ_FORMAT_AUTO, parameters.reads_path);
    seqs_sam = new SequenceFile();
    LOG_ALL("Aligning overlaps.\n");
    AlignMHAP(seqs_gfa, seqs_reads, overlaps_final, parameters.num_threads, *seqs_sam);
  }

  // Sanity check to see if the reads have quality values.
  for (int64_t i=0; i<seqs_sam->get_sequences().size(); i++) {
    if (seqs_sam->get_sequences()[i]->get_quality() == NULL || seqs_sam->get_sequences()[i]->get_quality_length() == 0) {
      fprintf (stderr, "ERROR: Reads are not specified in a format which contains quality information. Exiting.\n");
      exit(1);
    }
  }

  ConsensusDirectFromAln(parameters, seqs_gfa, *seqs_sam);
  seqs_sam->Clear();
  if (seqs_sam) { delete seqs_sam; }

//  std::string alt_contig_path = argv[3];

//  std::string mhap_path = "results/temp/consensus-lambda_30x_ont-mhap-iter1.mhap";
//  std::string reads_path = "test-data/lambda/reads.fastq";
//  std::vector<MHAPLine> overlaps, overlaps_filtered;
//  ParseMHAP(mhap_path, overlaps);
//  FilterMHAP(overlaps, overlaps_filtered);
//  SequenceFile seqs_reads(SEQ_FORMAT_AUTO, reads_path);
//  SequenceFile seqs_generated_sam;
//  AlignMHAP(seqs_gfa, seqs_reads, overlaps_filtered, seqs_generated_sam);
//
//  //  Consensus(parameters, seqs_gfa, seqs_sam);
//  ConsensusDirectFromAln(parameters, seqs_gfa, seqs_generated_sam);

	return 0;
}
