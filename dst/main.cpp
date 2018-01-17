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
        fprintf(stderr, "racon:: error: missing input file(s)!\n");
        help();
        exit(1);
    }

    //auto polisher = racon::createPolisher(input_paths[0], input_paths[1],
        //input_paths[2]);

    return 0;
}

void help() {
    printf(
        "usage: racon [options ...] <sequences> <mappings> <target sequences>\n"
        "\n"
        "    <sequences>\n"
        "        input file in FASTA/FASTQ format containing seqeunces used for\n"
        "        correction\n"
        "    <mappings>\n"
        "        input file in MHAP/PAF/SAM format containing mappings between\n"
        "        sequences and target sequences\n"
        "    <target seqeunces>\n"
        "        input file in FASTA/FASTQ format containing sequences which will\n"
        "        be corrected\n"
        "\n"
        "    options:\n");
}
