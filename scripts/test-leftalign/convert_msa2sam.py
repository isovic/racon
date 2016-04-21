#! /usr/bin/python

import os;
import sys;
import subprocess;

import fastqparser;
import utility_sam;

import traceback;

from msa2sam import *

def load_seqs(infile):
	[headers, seqs, quals] = fastqparser.read_fastq(infile);
	# try:
	# 	fp = open(infile, 'r');
	# except:
	# 	sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n');
	# 	exit(1);
	# alns = [line.strip() for line in fp.readlines() if len(line.strip()) != 0];
	# fp.close();
	# return alns;
	return seqs;

def main(msa_file, leftaligned_msa_file):
#	[folder, filename] = os.path.splitext(msa_file);
	folder = os.path.dirname(msa_file);
	filename = os.path.splitext(os.path.basename(msa_file))[0];

	temp_folder = '%s/temp-msa-leftalign' % (folder);

	if (not os.path.exists(temp_folder)):
		os.makedirs(temp_folder);

	alns = load_seqs(msa_file);
	sam_file = '%s/%s.sam' % (temp_folder, filename);
	bam_file = '%s/%s.bam' % (folder, filename);
	leftaligned_bam_file = '%s/%s.leftaligned.bam' % (temp_folder, filename);
	leftaligned_sam_file = '%s/%s.leftaligned.sam' % (temp_folder, filename);
#	leftaligned_msa_file = '%s/temp/%s.leftaligned.txt' % (temp_folder, filename);
	ref_file = '%s/%s.fasta' % (temp_folder, filename);

	fp_sam = open(sam_file, 'w');
	fp_ref = open(ref_file, 'w');
	fp_sam.write('@HD\tVN:1.0\tSO:unknown\n');

	sam_lines = [];
	ref_headers = [];
	ref_seqs = [];
	# for aln_id in xrange(0, len(alns)):
	# 	fp_sam.write('@SQ	SN:gi|48994873|gb|U00096.2|_Escherichia_coli_str._K-12_substr._MG1655,_complete_genome	LN:4639675');

	if (len(alns) < 2):
		sys.stderr.write('Not enough alignments in the MSA file! len(alns) = %d\n' % (len(alns)));
		exit(1);

	for aln_id in xrange(0, len(alns)):
		alns[aln_id] = alns[aln_id].upper();

	for aln_id in xrange(1, len(alns)):
		[ref_header, ref, sam_line] = convert_to_sam(alns[0], alns[aln_id], 'read_0', ('read_%d' % aln_id));

		ref_headers.append(ref_header);
		ref_seqs.append(ref);
		sam_lines.append(sam_line);

	fp_ref.write('>%s\n%s\n' % (ref_headers[0], ref_seqs[0]));
	fp_sam.write('@SQ\tSN:%s\tLN:%d\n' % (ref_headers[0], len(ref_seqs[0])));
	for i in xrange(len(sam_lines)):
		fp_sam.write('%s\n' % (sam_lines[i]));
	fp_sam.close();
	fp_ref.close();

	command = 'samtools faidx %s' % (ref_file);
	execute_command(command, sys.stderr, dry_run=False);

	command = 'samtools view -Sb -h %s > %s' % (sam_file, bam_file);
	execute_command(command, sys.stderr, dry_run=False);

	command = 'cat %s | bamleftalign/bamleftalign -d -f %s > %s' % (bam_file, ref_file, leftaligned_bam_file);
	execute_command(command, sys.stderr, dry_run=False);

	command = 'samtools view -h %s > %s' % (leftaligned_bam_file, leftaligned_sam_file);
	execute_command(command, sys.stderr, dry_run=False);

#	convert_sam_to_msa_fasta(ref_file, leftaligned_sam_file, leftaligned_msa_file);
	convert_rname_sams_to_msa(ref_file, leftaligned_sam_file, 'read_0', leftaligned_msa_file);

if __name__ == "__main__":
#	msa_file = 'data/msa1.txt';
#	msa_file = 'data/msa2.txt';
	# msa_file = 'data/window.fasta.1.msa';
	# msa_file = 'data/window.fasta.3.msa';
#	msa_file = 'data/window.fasta.5.msa';
	if (len(sys.argv) != 3):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('  %s <msa_file> <out_leftaligned_msa>\n' % (sys.argv[0]));
		exit(1);

	msa_file = sys.argv[1];
	out_msa_file = sys.argv[2];

	main(msa_file, out_msa_file);

# ../../../samscripts/src/fastqfilter.py touppercase data/window.fasta.5.msa data/window.fasta.5.uppercase.msa
# dnadiff -p data/temp/dnadiff/normal ../../temp/NC_001416.fa data/temp/window.fasta.5.consensus.fasta
# dnadiff -p data/temp/dnadiff/leftaligned ../../temp/NC_001416.fa data/temp/window.fasta.5.leftaligned.consensus.fasta
