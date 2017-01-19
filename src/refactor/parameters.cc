#include "parameters.h"
#include <iostream>
#include "log_system/log_system.h"
#include "argparser.h"

namespace is{

Parameters::Parameters(int argc, char* argv[]) :	overlap_format_(OverlapFormat::Paf()),
							window_len_(500)
							  {

  bool help = false;

  ArgumentParser argparser;

  // Input/output options.
  argparser.AddArgument(&(reads_path_), VALUE_TYPE_STRING, "", "reads", "", "Reads in FASTQ format.", -4, "Input/Output options");
  argparser.AddArgument(&(aln_path_), VALUE_TYPE_STRING, "", "alnpath", "", "Path to a MHAP file with read-to-target overlaps.", -3, "Input/Output options");
  argparser.AddArgument(&(raw_contigs_path_), VALUE_TYPE_STRING, "", "raw", "", "Path to the raw contig/read sequences (output from the layout step). GFA, FASTA/FASTQ or SAM formats allowed.", -2, "Input/Output options");
  argparser.AddArgument(&(consensus_path_), VALUE_TYPE_STRING, "", "out", "", "Output consensus sequence.", -1, "Input/Output options");
  argparser.AddArgument(&(is_sam_), VALUE_TYPE_BOOL, "", "sam", "0", "SAM is provided instead of MHAP. The reads file will be ignored, and seq and qual fields from the SAM file will be used.", 0, "Input/Output options");
  argparser.AddArgument(&(is_mhap_), VALUE_TYPE_BOOL, "", "mhap", "0", "Overlaps are in PAF format instead of MHAP.", 0, "Input/Output options");

  // Algorithmic options.
  argparser.AddArgument(&(qv_threshold_), VALUE_TYPE_DOUBLE, "", "bq", "10.0", "Threshold for the average base quality of the input reads. If a read has average BQ < specified, the read will be skipped. If value is < 0.0, filtering is disabled.", 0, "Algorithm");
  argparser.AddArgument(&(use_contig_qvs_), VALUE_TYPE_BOOL, "", "use-contig-qv", "0", "If false, dummy QVs equal to '!' will be assigned to each contig base during window consensus. Otherwise, QVs will be loaded from the contigs file if the file is in FASTQ format.", 0, "Algorithm");
  argparser.AddArgument(&(window_len_), VALUE_TYPE_INT64, "w", "winlen", "500", "Length of the window to perform POA on.", 0, "Algorithm");
  argparser.AddArgument(&(num_threads_), VALUE_TYPE_INT32, "t", "threads", "4", "Number of threads to use.", 0, "Control");
  argparser.AddArgument(&(batch_of_windows_), VALUE_TYPE_INT64, "b", "winbatch", "20000", "Size of the batch in which to process windows. After a batch is finished, consensus of the windows is joined and output to file.", 0, "Control");
  argparser.AddArgument(&(num_batches_), VALUE_TYPE_INT64, "", "num-batches", "-1", "The number of batches which to process", 0, "Control");
  argparser.AddArgument(&(start_window_), VALUE_TYPE_INT64, "", "start-window", "0", "ID of the window to start processing from.", 0, "Control");
  argparser.AddArgument(&(win_ovl_margin_), VALUE_TYPE_DOUBLE, "", "ovl-margin", "0.0", "Fraction of the window size to overlap the windows by.", 0, "Algorithm");
  argparser.AddArgument(&(do_erc_), VALUE_TYPE_BOOL, "", "erc", "0", "Perform error-correction instead of contig consensus. The only difference is in the type of parallelization to achieve better performance.", 0, "Control");
  argparser.AddArgument(&(error_rate_), VALUE_TYPE_DOUBLE, "e", "error-rate", "0.30", "Maximum allowed error rate. Used for filtering faulty overlaps.", 0, "Algorithm");

  // Alignment options.
  argparser.AddArgument(&(match_), VALUE_TYPE_INT32, "M", "match", "5", "Match score (positive value).", 0, "Alignment");
  argparser.AddArgument(&(mismatch_), VALUE_TYPE_INT32, "X", "mismatch", "-4", "Mismatch penalty (negative value expected).", 0, "Alignment");
  argparser.AddArgument(&(gap_open_), VALUE_TYPE_INT32, "G", "gapopen", "-8", "Gap open penalty (negative value expected).", 0, "Alignment");
  argparser.AddArgument(&(gap_ext_), VALUE_TYPE_INT32, "E", "gapext", "-6", "Gap extend penalty (negative value expected).", 0, "Alignment");

  // Other options.
  argparser.AddArgument(&(verbose_level_), VALUE_TYPE_INT32, "v", "verbose", "5", "Verbose level. 0 off, 1 low, 2 medium, 3 high, 4 and 5 all levels, 6-9 debug.", 0, "Other");
  argparser.AddArgument(&help, VALUE_TYPE_BOOL, "h", "help", "0", "View this help.", 0, "Other");

  // Process the command line arguments.
  argparser.ProcessArguments(argc, argv);


  std::string usage_cmd = std::string(argv[0]) +
                std::string(" [options] <reads.fastq> <overlaps.paf> <raw_contigs.fasta> <out_consensus.fasta>");

  // Program was run with no parameters. Verbose usage and exit.
  if (argc == 1) {
    std::cerr << "  " << usage_cmd << std::endl << std::endl;
    std::cerr << argparser.VerboseUsage() << std::endl;
    exit(1);
  }

  /// Check if help was triggered.
  if (argparser.GetArgumentByLongName("help")->is_set == true) {
  	std::cerr << "  " << usage_cmd << std::endl << std::endl;
  	std::cerr << argparser.VerboseUsage() << std::endl;
    exit(1);
  }  

  // Store the command line arguments for later use.
  for (int32_t i=0; i<argc; i++) {
  	cmd_arguments_.push_back(argv[i]);
  }
  program_folder_ = cmd_arguments_[0].substr(0, cmd_arguments_[0].find_last_of("\\/"));
  program_bin_ = cmd_arguments_[0].substr(cmd_arguments_[0].find_last_of("\\/") + 1);

  if (is_mhap_ == false && is_sam_ == false) {
    is_paf_ = true;
  }

  // Sanity check for input format specification.
  if (is_sam_ == true && is_paf_ == true) {
  	std::cerr << "ERROR: More than one input overlap/alignment format specified. Exiting." << std::endl;
    exit(1);
  }
  if (is_sam_ == false && reads_path_ == "") {
  	std::cerr << "ERROR: Reads file not specified! Exiting." << std::endl;
    exit(1);
  }

  // Adjust the verbose level for logging.
  if (verbose_level_ == 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_STD;
  } else if (verbose_level_ > 1) {
    LogSystem::GetInstance().LOG_VERBOSE_TYPE = LOG_VERBOSE_FULL | LOG_VERBOSE_STD;
  }

  // Set the verbose level for the execution of this program.
  LogSystem::GetInstance().SetProgramVerboseLevelFromInt(verbose_level_);
}

}
