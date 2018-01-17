#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <string>
#include <vector>

#include "polisher.hpp"

static struct option options[] = {
    {"fragment-correction", no_argument, 0, 'f'},
    {"window-length", required_argument, 0, 'w'},
    {"quality-threshold", required_argument, 0, 'q'},
    {"error-threshold", required_argument, 0, 'e'},
    {"match", required_argument, 0, 'm'},
    {"mismatch", required_argument, 0, 'x'},
    {"gap", required_argument, 0, 'g'},
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::vector<std::string> input_paths;

    uint32_t window_length = 500;
    double quality_threshold = 10.0;
    double error_threshold = 0.3;

    int8_t match = 5;
    int8_t mismatch = -4;
    int8_t gap = -8;
    uint32_t type = 0;

    uint32_t num_threads = 1;

    char argument;
    while ((argument = getopt_long(argc, argv, "fw:q:e:m:x:g:t:h", options, nullptr)) != -1) {
        switch (argument) {
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
            case 'h':
            default:
                help();
                exit(1);
        }
    }

    for (int32_t i = optind; i < argc; ++i) {
        input_paths.emplace_back(argv[i]);
    }

    if (input_paths.size() < 3) {
        fprintf(stderr, "racon:: error: missing input file(s)!\n");
        help();
        exit(1);
    }

    auto polisher = racon::createPolisher(input_paths[0], input_paths[1],
        input_paths[2], type == 0 ? racon::PolisherType::kC :
        racon::PolisherType::kF, window_length, quality_threshold,
        error_threshold, match, mismatch, gap, num_threads);

    return 0;
}

void help() {
    printf(
        "usage: racon [options ...] <sequences> <overlaps> <target sequences>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format containing seqeunces used for\n"
        "        correction\n"
        "    <overlaps>\n"
        "        input file in MHAP/PAF/SAM format containing overlaps between\n"
        "        sequences and target sequences\n"
        "    <target seqeunces>\n"
        "        input file in FASTA/FASTQ format containing sequences which will\n"
        "        be corrected\n"
        "\n"
        "    options:\n"
        "        -f, --fragment-correction\n"
        "            perform fragment correction instead of contig polishing\n"
        "        -w, --window-length <int>\n"
        "            default: 500\n"
        "            size of window on which POA is performed\n"
        "        -m, --match <int>\n"
        "            default: 5\n"
        "            score for matching bases\n"
        "        -x, --mismatch <int>\n"
        "            default: -4\n"
        "            score for mismatching bases\n"
        "        -g, --gap <int>\n"
        "            default: -8\n"
        "            gap penalty (must be negative)\n"
        "        -t, --threads <int>\n"
        "            default: 1\n"
        "            number of threads\n"
        "        -h, --help\n"
        "            prints out the help\n");
}
