#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "cigar.hpp"
#include "polisher.hpp"
#include "sequence.hpp"
#include "util.hpp"
#include "vcf.hpp"
#ifdef CUDA_ENABLED
#include "cuda/cudapolisher.hpp"
#endif

#ifndef RACON_VERSION
#error "Undefined version for Racon. Please pass version using -DRACON_VERSION macro."
#endif

static const char* version = RACON_VERSION;
static const int32_t CUDAALIGNER_INPUT_CODE = 10000;
static const int32_t CUDAALIGNER_BAND_WIDTH_INPUT_CODE = 10001;

static struct option options[] = {
    {"include-unpolished", no_argument, 0, 'u'},
    {"fragment-correction", no_argument, 0, 'f'},
    {"window-length", required_argument, 0, 'w'},
    {"quality-threshold", required_argument, 0, 'q'},
    {"error-threshold", required_argument, 0, 'e'},
    {"no-trimming", no_argument, 0, 'T'},
    {"liftover", required_argument, 0, 'L'},
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 'x'},
    {"gap", required_argument, 0, 'g'},
    {"threads", required_argument, 0, 't'},
    {"bed", required_argument, 0, 'B'},
    {"version", no_argument, 0, 'v'},
    {"help", no_argument, 0, 'h'},
#ifdef CUDA_ENABLED
    {"cudapoa-batches", optional_argument, 0, 'c'},
    {"cuda-banded-alignment", no_argument, 0, 'b'},
    {"cudaaligner-batches", required_argument, 0, CUDAALIGNER_INPUT_CODE},
    {"cudaaligner-band-width", required_argument, 0, CUDAALIGNER_BAND_WIDTH_INPUT_CODE},
#endif
    {0, 0, 0, 0}
};

enum class LiftoverOutFormat {
    PAF,
    SAM,
    VCF,
    Unknown,
};

LiftoverOutFormat ParseLiftoverFormatFromExt(const std::string& out_file) {
    std::vector<std::string> tokens = racon::Tokenize(out_file, '.');
    if (tokens.size() < 2) {
        return LiftoverOutFormat::Unknown;
    }
    std::string ext = tokens.back();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                        [](unsigned char c){ return std::tolower(c); });
    if (ext == "paf") {
        return LiftoverOutFormat::PAF;
    } else if (ext == "sam") {
        return LiftoverOutFormat::SAM;
    } else if (ext == "vcf") {
        return LiftoverOutFormat::VCF;
    }
    return LiftoverOutFormat::Unknown;
}

void WriteLiftoverFile(FILE *fp_out, LiftoverOutFormat liftover_fmt, const std::unique_ptr<racon::Polisher>& polisher,
                        const std::vector<std::unique_ptr<racon::Sequence>>& polished_sequences) {
    // Header.
    if (liftover_fmt == LiftoverOutFormat::SAM) {
        fprintf(fp_out, "@HD\tVN:1.5\n");
        for (const auto& it: polished_sequences) {
            const std::string cons_header = racon::TokenizeToWhitespaces(it->name())[0];
            fprintf(fp_out, "@SQ\tSN:%s\tLN:%lu\n", cons_header.c_str(), it->data().size());
        }
    } else if (liftover_fmt == LiftoverOutFormat::VCF) {
        fprintf(fp_out, "##fileformat=VCFv4.2\n");
        fprintf(fp_out, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
        fprintf(fp_out, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
        // Draft contig headers.
        for (const auto& it: polished_sequences) {
            // Get the input draft sequence, needed for length.
            const auto& draft_seq = polisher->sequences()[it->id()];
            const std::string draft_header = racon::TokenizeToWhitespaces(draft_seq->name())[0];
            fprintf(fp_out, "##contig=<ID=%s,length=%lu>\n", draft_header.c_str(), draft_seq->data().size());
        }
        fprintf(fp_out, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFORMAT_VALS\n");
    }

    for (const auto& it: polished_sequences) {
        // Sanity check.
        if (it->id() < 0) {
            std::ostringstream oss;
            oss << "Invalid target id: " << it->id() << "\n";
            throw std::runtime_error(oss.str());
        }

        // Get the input draft sequence, needed for length.
        const auto& draft_seq = polisher->sequences()[it->id()];

        // Parse the headers only up to the first whitespace.
        const std::string cons_header = racon::TokenizeToWhitespaces(it->name())[0];
        const std::string draft_header = racon::TokenizeToWhitespaces(draft_seq->name())[0];

        // PAF output.
        if (liftover_fmt == LiftoverOutFormat::PAF) {
            fprintf(fp_out, "%s\t%lu\t%lu\t%lu\t+\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\t60\tcg:Z:%s\n",
                                cons_header.c_str(), it->data().size(), 0, it->data().size(),
                                draft_header.c_str(), draft_seq->data().size(), 0, draft_seq->data().size(),
                                it->data().size(), draft_seq->data().size(),
                                it->cigar().c_str());

        } else if (liftover_fmt == LiftoverOutFormat::SAM) {
            fprintf(fp_out, "%s\t0\t%s\t1\t60\t%s\t*\t0\t0\t%s\t*\n",
                                cons_header.c_str(), draft_header.c_str(),it->cigar().c_str(), it->data().c_str());

        } else if (liftover_fmt == LiftoverOutFormat::VCF) {
            // Get the VCF events.
            auto cigar = racon::ParseCigarString(it->cigar());
            const auto& qseq = it->data();
            const auto& tseq = draft_seq->data();
            auto vcf_diffs = ExtractVCFEventsFromCigarString(cigar, qseq, tseq);

            // Write out the VCF events.
            for (const auto& v: vcf_diffs) {
                fprintf(fp_out, "%s\t%d\t.\t%s\t%s\t%d\tPASS\t.\tGT\t%s\n", draft_header.c_str(), v.pos + 1, v.ref.c_str(), v.alt.c_str(), 60, "1/1");
            }

        } else {
            throw std::runtime_error("Currently unsupported output liftover format.");
        }
    }
}

void help();

int main(int argc, char** argv) {

    std::vector<std::string> input_paths;

    uint32_t window_length = 500;
    double quality_threshold = 10.0;
    double error_threshold = 0.3;
    bool trim = true;

    int8_t match = 3;
    int8_t mismatch = -5;
    int8_t gap = -4;
    uint32_t type = 0;

    bool drop_unpolished_sequences = true;
    uint32_t num_threads = 1;

    uint32_t cudapoa_batches = 0;
    uint32_t cudaaligner_batches = 0;
    uint32_t cudaaligner_band_width = 0;
    bool cuda_banded_alignment = false;

    std::string bed_file;
    std::string out_liftover_file;

    std::string optstring = "ufw:q:e:m:x:g:t:B:L:h";
#ifdef CUDA_ENABLED
    optstring += "bc::";
#endif

    int32_t argument;
    while ((argument = getopt_long(argc, argv, optstring.c_str(), options, nullptr)) != -1) {
        switch (argument) {
            case 'u':
                drop_unpolished_sequences = false;
                break;
            case 'f':
                type = 1;
                break;
            case 'w':
                window_length = atoi(optarg);
                break;
            case 'q':
                quality_threshold = atof(optarg);
                break;
            case 'e':
                error_threshold = atof(optarg);
                break;
            case 'T':
                trim = false;
                break;
            case 'L':
                out_liftover_file = optarg;
                break;
            case 'm':
                match = atoi(optarg);
                break;
            case 'x':
                mismatch = atoi(optarg);
                break;
            case 'g':
                gap = atoi(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'B':
                bed_file = std::string(optarg);
                break;
            case 'v':
                printf("%s\n", version);
                exit(0);
            case 'h':
                help();
                exit(0);
#ifdef CUDA_ENABLED
            case 'c':
                //if option c encountered, cudapoa_batches initialized with a default value of 1.
                cudapoa_batches = 1;
                // next text entry is not an option, assuming it's the arg for option 'c'
                if (optarg == NULL && argv[optind] != NULL
                    && argv[optind][0] != '-') {
                    cudapoa_batches = atoi(argv[optind++]);
                }
                // optional argument provided in the ususal way
                if (optarg != NULL) {
                    cudapoa_batches = atoi(optarg);
                }
                break;
            case 'b':
                cuda_banded_alignment = true;
                break;
            case CUDAALIGNER_INPUT_CODE: // cudaaligner-batches
                cudaaligner_batches = atoi(optarg);
                break;
            case CUDAALIGNER_BAND_WIDTH_INPUT_CODE: // cudaaligner-band-width
                cudaaligner_band_width = atoi(optarg);
                break;
#endif
            default:
                exit(1);
        }
    }

    for (int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 3) {
        fprintf(stderr, "[racon::] error: missing input file(s)!\n");
        help();
        exit(1);
    }

    std::cerr << "BED file: '" << bed_file << "'\n";

    // Prepare output for the liftover if required.
    const bool produce_liftover = (out_liftover_file.empty() ? false : true);
    LiftoverOutFormat liftover_fmt = LiftoverOutFormat::Unknown;
    FILE* fp_out_liftover = NULL;
    if (produce_liftover) {
        liftover_fmt = ParseLiftoverFormatFromExt(out_liftover_file);
        // Parse the output format, and make sure it's good.
        if (liftover_fmt == LiftoverOutFormat::Unknown) {
            throw std::runtime_error("Unknown output liftover format for file: '" + out_liftover_file + "'.");
        }
        // Open the output liftover file if required, and sanity check.
        fp_out_liftover = fopen(out_liftover_file.c_str(), "w");
        if (fp_out_liftover == NULL) {
            throw std::runtime_error("Cannot open file '" + out_liftover_file + "' for writing!.");
        }
    }

    auto polisher = racon::createPolisher(input_paths[0], input_paths[1],
        input_paths[2], bed_file, type == 0 ? racon::PolisherType::kC :
        racon::PolisherType::kF, window_length, quality_threshold,
        error_threshold, trim, produce_liftover, match, mismatch, gap, num_threads,
        cudapoa_batches, cuda_banded_alignment, cudaaligner_batches,
        cudaaligner_band_width);

    polisher->initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    polisher->polish(polished_sequences, drop_unpolished_sequences);

    for (const auto& it: polished_sequences) {
        fprintf(stdout, ">%s\n%s\n", it->name().c_str(), it->data().c_str());
    }

    // Write the liftover file if required.
    if (produce_liftover) {
        fprintf(stderr, "[racon::] Writing the liftover file.\n");
        WriteLiftoverFile(fp_out_liftover, liftover_fmt, polisher, polished_sequences);
        if (fp_out_liftover) {
            fclose(fp_out_liftover);
        }
    }

    return 0;
}

void help() {
    printf(
        "usage: racon [options ...] <sequences> <overlaps> <target sequences>\n"
        "\n"
        "    #default output is stdout\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences used for correction\n"
        "    <overlaps>\n"
        "        input file in MHAP/PAF/SAM format (can be compressed with gzip)\n"
        "        containing overlaps between sequences and target sequences\n"
        "    <target sequences>\n"
        "        input file in FASTA/FASTQ format (can be compressed with gzip)\n"
        "        containing sequences which will be corrected\n"
        "\n"
        "    options:\n"
        "        -u, --include-unpolished\n"
        "            output unpolished target sequences\n"
        "        -f, --fragment-correction\n"
        "            perform fragment correction instead of contig polishing\n"
        "            (overlaps file should contain dual/self overlaps!)\n"
        "        -w, --window-length <int>\n"
        "            default: 500\n"
        "            size of window on which POA is performed\n"
        "        -q, --quality-threshold <float>\n"
        "            default: 10.0\n"
        "            threshold for average base quality of windows used in POA\n"
        "        -e, --error-threshold <float>\n"
        "            default: 0.3\n"
        "            maximum allowed error rate used for filtering overlaps\n"
        "        --no-trimming\n"
        "            disables consensus trimming at window ends\n"
        "        -L, --liftover <string>\n"
        "            default: ''\n"
        "            optional output liftover file which converts the draft\n"
        "            sequence to the output consensus\n"
        "        -m, --match <int>\n"
        "            default: 3\n"
        "            score for matching bases\n"
        "        -x, --mismatch <int>\n"
        "            default: -5\n"
        "            score for mismatching bases\n"
        "        -g, --gap <int>\n"
        "            default: -4\n"
        "            gap penalty (must be negative)\n"
        "        -B, --bed <str>\n"
        "            default: ''\n"
        "            path to a BED file with regions to polish\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        --version\n"
        "            prints the version number\n"
        "        -h, --help\n"
        "            prints the usage\n"
#ifdef CUDA_ENABLED
        "        -c, --cudapoa-batches <int>\n"
        "            default: 0\n"
        "            number of batches for CUDA accelerated polishing per GPU\n"
        "        -b, --cuda-banded-alignment\n"
        "            use banding approximation for alignment on GPU\n"
        "        --cudaaligner-batches <int>\n"
        "            default: 0\n"
        "            number of batches for CUDA accelerated alignment per GPU\n"
        "        --cudaaligner-band-width <int>\n"
        "            default: 0\n"
        "            Band width for cuda alignment. Must be >= 0. Non-zero allows user defined \n"
        "            band width, whereas 0 implies auto band width determination.\n"
#endif
    );
}
