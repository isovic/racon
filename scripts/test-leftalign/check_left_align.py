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

def convert_sam_to_msa(ref_path, sam_path, out_msa_path):
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
		# print alnseq1;
		# print alnseq2;
		# # print alignment;
		# print '';

	fp.close();




def convert_to_sam(alnseq1, alnseq2, id):
	if (len(alnseq1) != len(alnseq2)): return -1;

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
	sam_line += 'read_%d_q\t' % (id);			# qname
	sam_line += '0\t';				# flag
	sam_line += 'read_%d_r\t' % (id);			# rname
	sam_line += '1\t';				# pos
	sam_line += '40\t';				# mapq
	sam_line += '%s\t' % cigar_string;	#cigar
	sam_line += '*\t';				# mrname
	sam_line += '0\t';				# mpos
	sam_line += '0\t';				# insert size
	sam_line += '%s\t' % (seq);	# seq
	sam_line += '*';				# qual

	ref_header = 'read_%d_r' % (id);

	return [ref_header, ref, sam_line];



# def check_left_alignment(aln, aln_id):
# 	if (len(aln) < 2):
# 		return;

# 	[ref_header, ref, sam_line] = convert_to_sam(aln[0], aln[1], aln_id);
# 	print sam_line;

# 	last_base_pos = [];
# 	for i in xrange(len(aln)):
# 		last_base_pos.append(-1);
# 	for i in xrange(0, len(aln[0])):
# 		pass;

def load_seqs(infile):
	try:
		fp = open(infile, 'r');
	except:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n');
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
		[ref_header, ref, sam_line] = convert_to_sam(aln[0], aln[1], aln_id);
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

	convert_sam_to_msa(ref_file, leftaligned_sam_file, leftaligned_msa_file);

if __name__ == "__main__":
	msa_file = 'data/msa1.txt';
	msa_file = 'data/msa2.txt';

	main(msa_file);
