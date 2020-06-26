#!/usr/bin/env python

from __future__ import print_function
import os, sys, time, shutil, argparse, subprocess

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#*******************************************************************************

class RaconWrapper:

    __racon = '@racon_path@'
    __rampler = '@rampler_path@'

    def __init__(self, sequences, overlaps, target_sequences, split, subsample,
        include_unpolished, fragment_correction, window_length, quality_threshold,
        error_threshold, match, mismatch, gap, threads):

        self.sequences = os.path.abspath(sequences)
        self.subsampled_sequences = None
        self.overlaps = os.path.abspath(overlaps)
        self.target_sequences = os.path.abspath(target_sequences)
        self.split_target_sequences = []
        self.chunk_size = split
        self.reference_length, self.coverage = subsample if subsample is not None\
            else (None, None)
        self.include_unpolished = include_unpolished
        self.fragment_correction = fragment_correction
        self.window_length = window_length
        self.quality_threshold = quality_threshold
        self.error_threshold = error_threshold
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.threads = threads
        self.work_directory = os.getcwd() + '/racon_work_directory_' + str(time.time())

    def __enter__(self):
        try:
            os.makedirs(self.work_directory)
        except OSError:
            if (not os.path.isdir(self.work_directory)):
                eprint('[RaconWrapper::__enter__] error: unable to create work directory!')
                sys.exit(1)

    def __exit__(self, exception_type, exception_value, traceback):
        try:
            shutil.rmtree(self.work_directory)
        except OSError:
            eprint('[RaconWrapper::__exit__] warning: unable to clean work directory!')

    def run(self):
        # run preprocess
        eprint('[RaconWrapper::run] preparing data with rampler')
        if (self.reference_length is not None and self.coverage is not None):
            try:
                p = subprocess.Popen([RaconWrapper.__rampler, '-o', self.work_directory,
                    'subsample', self.sequences, self.reference_length, self.coverage])
            except OSError:
                eprint('[RaconWrapper::run] error: unable to run rampler!')
                sys.exit(1)
            p.communicate()
            if (p.returncode != 0):
                sys.exit(1)

            base_name = os.path.basename(self.sequences).split('.')[0]
            extension = '.fasta' if (self.sequences.endswith('.fasta') or\
                self.sequences.endswith('.fasta.gz') or\
                self.sequences.endswith('.fa') or\
                self.sequences.endswith('.fa.gz')) else\
                '.fastq'
            self.subsampled_sequences = os.path.join(self.work_directory, base_name) +\
                '_' + self.coverage + 'x' + extension
            if (not os.path.isfile(self.subsampled_sequences)):
                eprint('[RaconWrapper::run] error: unable to find subsampled sequences!')
                sys.exit(1)
        else:
            self.subsampled_sequences = self.sequences

        if (self.chunk_size is not None):
            try:
                p = subprocess.Popen([RaconWrapper.__rampler, '-o', self.work_directory,
                    'split', self.target_sequences, self.chunk_size])
            except OSError:
                eprint('[RaconWrapper::run] error: unable to run rampler!')
                sys.exit(1)
            p.communicate()
            if (p.returncode != 0):
                sys.exit(1)

            base_name = os.path.basename(self.target_sequences).split('.')[0]
            extension = '.fasta' if (self.target_sequences.endswith('.fasta') or\
                self.target_sequences.endswith('.fasta.gz') or\
                self.target_sequences.endswith('.fa') or\
                self.target_sequences.endswith('.fa.gz')) else\
                '.fastq'

            i = 0
            while (True):
                target_sequences_part = os.path.join(self.work_directory, base_name) +\
                    '_' + str(i) + extension
                if (not os.path.isfile(target_sequences_part)):
                    break
                self.split_target_sequences.append(target_sequences_part)
                i += 1
            if (len(self.split_target_sequences) == 0):
                eprint('[RaconWrapper::run] error: unable to find split target sequences!')
                sys.exit(1)
        else:
            self.split_target_sequences.append(self.target_sequences)

        racon_params = [RaconWrapper.__racon]
        if (self.include_unpolished == True): racon_params.append('-u')
        if (self.fragment_correction == True): racon_params.append('-f')
        racon_params.extend(['-w', str(self.window_length),
            '-q', str(self.quality_threshold),
            '-e', str(self.error_threshold),
            '-m', str(self.match),
            '-x', str(self.mismatch),
            '-g', str(self.gap),
            '-t', str(self.threads),
            self.subsampled_sequences, self.overlaps, ""])

        for target_sequences_part in self.split_target_sequences:
            eprint('[RaconWrapper::run] processing data with racon')
            racon_params[-1] = target_sequences_part
            try:
                p = subprocess.Popen(racon_params)
            except OSError:
                eprint('[RaconWrapper::run] error: unable to run racon!')
                sys.exit(1)
            p.communicate()
            if (p.returncode != 0):
                sys.exit(1)

        self.subsampled_sequences = None
        self.split_target_sequences = []

#*******************************************************************************

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''Racon_wrapper encapsulates
        racon and adds two additional features from the outside to enable easier
        usage to the end-user. Sequences can now be subsampled to decrease the
        total execution time (accuracy might be lower) while target
        sequences can be split into smaller chunks and run sequentially to
        decrease memory consumption. Both features can be run at the same time
        as well! The usage equals the one of racon.''',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sequences', help='''input file in FASTA/FASTQ format
        (can be compressed with gzip) containing sequences used for correction''')
    parser.add_argument('overlaps', help='''input file in MHAP/PAF/SAM format
        (can be compressed with gzip) containing overlaps between sequences and
        target sequences''')
    parser.add_argument('target_sequences', help='''input file in FASTA/FASTQ
        format (can be compressed with gzip) containing sequences which will be
        corrected''')
    parser.add_argument('--split', help='''split target sequences into chunks of
        desired size in bytes''')
    parser.add_argument('--subsample', nargs=2, help='''subsample sequences to
        desired coverage (2nd argument) given the reference length (1st argument)''')
    parser.add_argument('-u', '--include-unpolished', action='store_true',
        help='''output unpolished target sequences''')
    parser.add_argument('-f', '--fragment-correction', action='store_true',
        help='''perform fragment correction instead of contig polishing
        (overlaps file should contain dual/self overlaps!)''')
    parser.add_argument('-w', '--window-length', default=500, help='''size of
        window on which POA is performed''')
    parser.add_argument('-q', '--quality-threshold', default=10.0,
        help='''threshold for average base quality of windows used in POA''')
    parser.add_argument('-e', '--error-threshold', default=0.3, help='''maximum
        allowed error rate used for filtering overlaps''')
    parser.add_argument('-m', '--match', default=5, help='''score for matching
        bases''')
    parser.add_argument('-x', '--mismatch', default=-4, help='''score for
        mismatching bases''')
    parser.add_argument('-g', '--gap', default=-8, help='''gap penalty (must be
        negative)''')
    parser.add_argument('-t', '--threads', default=1, help='''number of threads''')

    args = parser.parse_args()

    racon = RaconWrapper(args.sequences, args.overlaps, args.target_sequences,
        args.split, args.subsample, args.include_unpolished,
        args.fragment_correction, args.window_length, args.quality_threshold,
        args.error_threshold, args.match, args.mismatch, args.gap, args.threads)

    with racon:
        racon.run()
