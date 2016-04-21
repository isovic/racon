#! /usr/bin/python

import os;
import sys;
import subprocess;

import fastqparser;
import utility_sam;

import traceback;
def execute_command(command, fp_log, dry_run=True):
    if (dry_run == True):
        fp_log.write('Executing (dryrun): "%s".\n' % (command));
        fp_log.write('\n');
        return 0;
    if (dry_run == False):
        fp_log.write('Executing: "%s".\n' % (command));
        rc = subprocess.call(command, shell=True);
        if (rc != 0):
            fp_log.write('ERROR: subprocess call returned error code: %d.\n' % (rc));
            fp_log.write('Traceback:\n');
            # traceback.print_stack(fp_log);
            # exit(1);
        return rc;

def load_seqs(infile):
	try:
		fp = open(infile, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % (infile));
		exit(1);

	alns = [[]];
	for line in fp:
		line = line.strip();
		if (line == ''):
			alns.append([]);
		else:
			alns[-1].append(line);
	fp.close();

	return alns;


def convert_to_sam(alnseq1, alnseq2, rname, qname):
	if (len(alnseq1) != len(alnseq2)):
		print ''
		sys.stderr.write('len(alnseq1) != len(alnseq2)\n');
		sys.stderr.write('len(alnseq1) = %d\n' % (len(alnseq1)));
		sys.stderr.write('len(alnseq2) = %d\n' % (len(alnseq2)));
		return None;

	OP_MATCH = 'M';
	OP_MISMATCH = 'M';
	OP_INS = 'I';
	OP_DEL = 'D';

	ref = '';
	seq = '';
	cigar = [];

	alignment = [];

	for i in xrange(0, len(alnseq1)):
		# sys.stdout.write('[i = %d] alnseq1[%d] = %s, alnseq2[%d] = %s' % (i, i, alnseq1[i], i, alnseq2[i]));
		op = None;
		if (alnseq1[i] == alnseq2[i] and alnseq1[i] != '-'):
			alignment.append(OP_MATCH);
			ref += (alnseq1[i]);
			seq += (alnseq2[i]);
			op = OP_MATCH;
		elif (alnseq1[i] == '-' and alnseq2[i] != '-'):
			alignment.append(OP_INS);
			seq += (alnseq2[i]);
			op = OP_INS;
		elif (alnseq1[i] != '-' and alnseq2[i] == '-'):
			alignment.append(OP_DEL);
			ref += (alnseq1[i]);
			op = OP_DEL;
		elif (alnseq1[i] != alnseq2[i]):
			alignment.append(OP_MISMATCH);
			ref += (alnseq1[i]);
			seq += (alnseq2[i]);
			op = OP_MISMATCH;
		if (op == None):
			continue;
		elif  (len(cigar) == 0):
			cigar.append([op, 1]);
		elif (op == cigar[-1][0]):
			cigar[-1][1] += 1;
		elif (op != cigar[-1][0]):
			cigar.append([op, 1]);

	cigar_string = ''.join(['%d%s' % (val[1], val[0]) for val in cigar]);

	sam_line = '';
#	sam_line += 'read_%d_q\t' % (id);			# qname
	sam_line += '%s\t' % (qname);
	sam_line += '0\t';				# flag
#	sam_line += 'read_%d_r\t' % (id);			# rname
	sam_line += '%s\t' % (rname);
	sam_line += '1\t';				# pos
	sam_line += '40\t';				# mapq
	sam_line += '%s\t' % cigar_string;	#cigar
	sam_line += '*\t';				# mrname
	sam_line += '0\t';				# mpos
	sam_line += '0\t';				# insert size
	sam_line += '%s\t' % (seq);	# seq
	sam_line += '*';				# qual

#	ref_header = 'read_%d_r' % (id);
	ref_header = '%s' % (rname);

	return [ref_header, ref, sam_line];	

### Returns: alnseq1 (reference), alnseq2 (query) and their alignment.
def sam_line_to_msa(ref_seq, sam_line):
	alnseq1 = '';
	alnseq2 = '';
	alignment = '';

	split_ops = sam_line.CalcCigarStartingPositions(separate_matches_in_individual_bases=True);
	for split_op in split_ops:
		[cigar_count, cigar_op, pos_on_ref, pos_on_read] = split_op;
		# print split_op;
		if (cigar_op in 'M=X'):
			if (ref_seq[pos_on_ref] == sam_line.seq[pos_on_read]):
				alignment += '|';
			else:
				alignment += 'X';
			alnseq1 += ref_seq[pos_on_ref];
			alnseq2 += sam_line.seq[pos_on_read];
		elif (cigar_op == 'I'):
			alnseq1 += '-' * cigar_count;
			alnseq2 += sam_line.seq[pos_on_read:(pos_on_read+cigar_count)];
			alignment += ' ' * cigar_count;
		elif (cigar_op == 'D'):
			alnseq1 += ref_seq[pos_on_ref:(pos_on_ref+cigar_count)];
			alnseq2 += '-' * cigar_count;
			alignment += ' ' * cigar_count;

	return [alnseq1, alnseq2, alignment];

def convert_sam_to_msa_pairwise(ref_path, sam_path, out_msa_path):
	[headers_sam, sam_lines] = utility_sam.LoadSAM(sam_path);
	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	fp = open(out_msa_path, 'w');

	seqs_ref_hash = {};
	for i in xrange(0, len(seqs_ref)):
		seqs_ref_hash[headers_ref[i]] = seqs_ref[i];
		seqs_ref_hash[headers_ref[i].split()[0]] = seqs_ref[i];

	for sam_line in sam_lines:
		[alnseq1, alnseq2, alignment] = sam_line_to_msa(seqs_ref_hash[sam_line.rname], sam_line);
		fp.write('%s\n' % (alnseq1));
		fp.write('%s\n' % (alnseq2));
		fp.write('\n');

	fp.close();

def main(msa_file):
	alns = load_seqs(msa_file);
	sam_file = '%s.sam' % (os.path.splitext(msa_file)[0]);
	bam_file = '%s.bam' % (os.path.splitext(msa_file)[0]);
	leftaligned_bam_file = '%s.leftaligned.bam' % (os.path.splitext(msa_file)[0]);
	leftaligned_sam_file = '%s.leftaligned.sam' % (os.path.splitext(msa_file)[0]);
	leftaligned_msa_file = '%s.leftaligned.txt' % (os.path.splitext(msa_file)[0]);
	ref_file = '%s.fasta' % (os.path.splitext(msa_file)[0]);

	fp_sam = open(sam_file, 'w');
	fp_ref = open(ref_file, 'w');
	fp_sam.write('@HD\tVN:1.0\tSO:unknown\n');

	sam_lines = [];
	ref_headers = [];
	ref_seqs = [];
	# for aln_id in xrange(0, len(alns)):
	# 	fp_sam.write('@SQ	SN:gi|48994873|gb|U00096.2|_Escherichia_coli_str._K-12_substr._MG1655,_complete_genome	LN:4639675');

	for aln_id in xrange(0, len(alns)):
		aln = alns[aln_id];
		if (len(aln) != 2):
			sys.stderr.write('len(aln) = %d, should be 2. Continuing.\n' % (len(aln)));
			continue;
		# print aln[0]
		# print aln[1]
		rname = 'read_%d_r' % aln_id;
		qname = 'read_%d_q' % aln_id;
		# print aln_id;
		converted = convert_to_sam(aln[0], aln[1], rname, qname);
		if (converted == None):
			continue;
		[ref_header, ref, sam_line] = converted;
		ref_headers.append(ref_header);
		ref_seqs.append(ref);
		sam_lines.append(sam_line);

	for i in xrange(len(ref_seqs)):
		fp_ref.write('>%s\n%s\n' % (ref_headers[i], ref_seqs[i]));
		fp_sam.write('@SQ\tSN:%s\tLN:%d\n' % (ref_headers[i], len(ref_seqs[i])));
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

	convert_sam_to_msa_pairwise(ref_file, leftaligned_sam_file, leftaligned_msa_file);

if __name__ == "__main__":
	# msa_file = 'data/msa1.txt';
	# msa_file = 'data/msa2.txt';
	# msa_file = 'data/msa3.txt';
	# msa_file = 'data2/lambda_500_pairwise.txt';

	if (len(sys.argv) != 2):
		sys.stderr.write('Usage:\n');
		sys.stderr.write('  %s <in_msa.txt>\n' % (sys.argv[0]));
		sys.stderr.write('\n');
		sys.stderr.write('Output will be stored in two files in the same path as the input file, with extensions changed to: <in_msa.sam> and <in_msa.leftaligned.sam>. The first one is the direct conversion from the input MSA file to SAM format, and the other is its leftaligned version.\n');

		sys.stderr.write('This script loads MSA alignments in the format:\n');
		sys.stderr.write('aligned_seqs_1_a\n');
		sys.stderr.write('aligned_seqs_1_b\n');
		sys.stderr.write('\n');
		sys.stderr.write('aligned_seqs_2_a\n');
		sys.stderr.write('aligned_seqs_2_b\n');
		sys.stderr.write('...\n\n');
		sys.stderr.write('There are no headers, and all MSA should be pairwise (only two rows per alignment). The first row is seq1 in a pair, and the second is seq2.\n');
		sys.stderr.write('\n');

		exit(1);

	msa_file = sys.argv[1];

	main(msa_file);
