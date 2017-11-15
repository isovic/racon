#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <getopt.h>

#include <string>
#include <vector>

#include "polisher.hpp"

static struct option options[] = {
    {"quality-threshold", required_argument, 0, 'q'},
    {"threads", required_argument, 0, 't'},
    {"help", no_argument, 0, 'h'},
    {0, 0, 0, 0}
};

void help();

int main(int argc, char** argv) {

    std::vector<std::string> input_paths;

    double quality_threshold = 10.0;
    uint32_t num_threads = 1;

    char argument;
    while ((argument = getopt_long(argc, argv, "q:t:h", options, nullptr)) != -1) {
        switch (argument) {
            case 'q':
                quality_threshold = atof(optarg);
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
        help();
        exit(1);
    }

    auto polisher = racon::createPolisher();

    return 0;
}

void help() {
    printf(
        "usage: racon [options ...] <reads> <mappings> <contigs/target reads>\n"
        "\n"
        "    <reads>\n"
        "        input file in FASTA/FASTQ format containing reads\n"
        "    <mappings>\n"
        "        input file in MHAP/PAF format containing mappings to contigs or\n"
        "        targer reads\n"
        "    <contigs/target reads>\n"
        "        input file in FASTA/FASTQ format containing contigs or reads\n"
        "        which will be corrected\n"
        "\n"
        "    options:\n");
}
