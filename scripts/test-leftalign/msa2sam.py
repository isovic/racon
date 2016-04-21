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

def convert_sam_to_msa_fasta(ref_path, sam_path, out_msa_path):
	[headers_sam, sam_lines] = utility_sam.LoadSAM(sam_path);
	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	fp = open(out_msa_path, 'w');

	seqs_ref_hash = {};
	for i in xrange(0, len(seqs_ref)):
		seqs_ref_hash[headers_ref[i]] = seqs_ref[i];
		seqs_ref_hash[headers_ref[i].split()[0]] = seqs_ref[i];

	aln_id = 0;
	for sam_line in sam_lines:
		[alnseq1, alnseq2, alignment] = sam_line_to_msa(seqs_ref_hash[sam_line.rname], sam_line);
		if (aln_id == 0):
			fp.write('>Ref\n%s\n' % (alnseq1));
			is_ref_output = True;
		fp.write('>Read_%d\n%s\n' % (aln_id, alnseq2));
		aln_id += 1;

	fp.close();


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

class PileupBase:
	def __init__(self, op, count, qid, depth_id, seq):
		self.op = op;
		self.count = count;
		self.qid = qid;			# Query ID of the alignment in the SAM file.
		self.seq = seq;
		self.id = depth_id;		# Ordinal number of the alignment in the MSA. If all alignments from the SAM file are taken into account, then self.id == self.qid. Otherwise, qid might be larger.
	# def verbose(self):
	# 	ret = '%'

def convert_rname_sams_to_msa(ref_path, sam_path, rname_to_convert, out_msa_file):
	[headers_sam, sam_lines] = utility_sam.LoadSAM(sam_path);
	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	ref_sam = utility_sam.SAMLine();
	for i in xrange(0, len(seqs_ref)):
		if (rname_to_convert == headers_ref[i]):
			ref_sam.qname = rname_to_convert;
			ref_sam.pos = 1;
			ref_sam.clipped_pos = 1;
			ref_sam.rname = rname_to_convert;
			ref_sam.cigar = '%dM' % (len(seqs_ref[i]));
			ref_sam.flag = 0;
			ref_sam.mapq = 40;
			ref_sam.seq = seqs_ref[i];
			ref_sam.qual = '*';
	if (ref_sam.qname == ''):
		sys.stderr.write('ERROR: rname not found in reference FASTQ file!\n');
		return 1;

	pileup = [];
	all_ref_cigar_pos = ref_sam.CalcCigarStartingPositions(separate_matches_in_individual_bases=True);
	for i in xrange(0, len(all_ref_cigar_pos)):
		[cigar_count, cigar_op, pos_on_ref, pos_on_read] = all_ref_cigar_pos[i];
		pileup.append([PileupBase('M', cigar_count, -1, 0, ref_sam.seq[pos_on_read])]);
		if (pos_on_ref != i):
			sys.stderr.write('ERROR: Something went wrong, position on reference is not consecutive to the M CIGAR operations! i = %d, pos_on_ref = %d\n' % (i, pos_on_ref));
			sys.stderr.write('%s\n' % str(all_ref_cigar_pos[i]));
			exit(1);
	
	pileup.append([PileupBase('M', 1, -1, 0, '-')]);
	# for i in xrange(0, len(pileup)):
	# 	print '%d\t%s' % (i, str(pileup[i][0].seq));

#	qid = 0;

#	for sam_line in sam_lines:
	cov_depth = 0;
	msa_headers = [rname_to_convert];
	for qid in xrange(0, len(sam_lines)):
		sam_line = sam_lines[qid];

		if (sam_line.rname == rname_to_convert):
			cov_depth += 1;
			msa_headers.append(sam_line.qname);

			all_cigar_pos = sam_line.CalcCigarStartingPositions(separate_matches_in_individual_bases=True);

			### Check if the alignment starts at the beginning of the 'reference'. If not, fill the gap with '-'.
			if (len(all_cigar_pos) > 0):
				[cigar_count, cigar_op, pos_on_ref, pos_on_read] = all_cigar_pos[0];
				for i in xrange(0, pos_on_ref):
					pileup[i].append(PileupBase('M', 1, qid, (cov_depth), '-'));

			i = 0;
			for cigar_pos in all_cigar_pos:
				i += 1;	
				# sys.stderr.write('%s, i = %d, len(all_cigar_pos) = %d\n' % (str(cigar_pos), i, len(all_cigar_pos)));
				[cigar_count, cigar_op, pos_on_ref, pos_on_read] = cigar_pos;
				if (cigar_op == 'M' or cigar_op == '=' or cigar_op == 'X'):
					pileup[pos_on_ref].append(PileupBase('M', cigar_count, qid, (cov_depth), sam_line.seq[pos_on_read:(pos_on_read+cigar_count)]));
				elif (cigar_op == 'I' or cigar_op == 'S'):
					pileup[pos_on_ref].append(PileupBase('I', cigar_count, qid, (cov_depth), sam_line.seq[pos_on_read:(pos_on_read+cigar_count)]));
				elif (cigar_op == 'D'):
					for di in xrange(0, cigar_count):
						pileup[pos_on_ref+di].append(PileupBase('M', 1, qid, (cov_depth), '-'));

			### Check if the alignment ends at the end of the 'reference'. If not, fill the gap with '-'.
			if (len(all_cigar_pos) > 0):
				[cigar_count, cigar_op, pos_on_ref, pos_on_read] = all_cigar_pos[-1];
				for i in xrange((pos_on_ref + cigar_count), len(ref_sam.seq)):
					pileup[i].append(PileupBase('M', 1, qid, (cov_depth+500), '-'));

			### Add an extra base to every query line.
			# for i in xrange(1, len(pileup)):
			pileup[-1].append(PileupBase('M', 1, qid, (cov_depth+500), '-'));
	cov_depth += 1;

#		qid += 1;

	# for i in xrange(0, len(pileup)):
	# 	pileup[i] = sorted(pileup[i], key=lambda x: (x.op, -x.id, -len(x.seq)), reverse=True)
	# 	sys.stderr.write('%d\t' % i);
	# 	for k in xrange(0, len(pileup[i])):
	# 		if (k > 0): sys.stderr.write(',');
	# 		if (pileup[i][k].op == 'I'): sys.stderr.write('[%s, %s, %d]' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
	# 		else: sys.stderr.write('(%s, %s, %d)' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
	# 	sys.stderr.write('\n');

	### Initialize the MSA matrix for zero elements.
	msa_matrix = [];
	for i in xrange(0, cov_depth): msa_matrix.append('');

	for i in xrange(0, len(pileup)):
		### Sort first by the op so that matches/mismatches/deletion bases come first, and then inserions.
		### Sort the bases in ascending order of qid, and sort insertions in the ascending order of insertion lengths.
		pileup[i] = sorted(pileup[i], key=lambda x: (x.op, -len(x.seq), -x.id), reverse=True)
		### The number of non-insertion bases should be constant (this includes the '-' bases) and equal to cov_depth.
#		for j in xrange(0, cov_depth):
#			msa_matrix[j] = pileup[i][j].seq;
		# Expand the insertions to the length of the largest insertion, and add it to the corresponding M base.
		max_ins_len = len(pileup[i][-1].seq) if (pileup[i][-1].op == 'I') else 0;
		for j in xrange(cov_depth, len(pileup[i])):
			ordid = pileup[i][j].id;
			pileup[i][j].seq = ('-' * (max_ins_len - len(pileup[i][j].seq))) + pileup[i][j].seq;
			pileup[i][ordid].seq = pileup[i][j].seq + pileup[i][ordid].seq;
#			pileup[i][ordid].seq += pileup[i][j].seq + ('-' * (max_ins_len - len(pileup[i][j].seq) - 1));
			# print 'ordid = %d' % (ordid);

		# sys.stderr.write('%d\t' % i);
		# for k in xrange(0, len(pileup[i])):
		# 	if (k > 0): sys.stderr.write(',');
		# 	if (pileup[i][k].op == 'I'): sys.stderr.write('[%s, %s, %d]' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
		# 	else: sys.stderr.write('(%s, %s, %d)' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
		# sys.stderr.write('\n');

		for j in xrange(0, cov_depth):
			# try:
			if (len(pileup[i][j].seq) == 1):
				pileup[i][j].seq = ('-' * (max_ins_len)) + pileup[i][j].seq;
			msa_matrix[j] += pileup[i][j].seq;

#		if (i >= 0 and i < 30):
		# sys.stderr.write('%d\t' % i);
		# for k in xrange(0, len(pileup[i])):
		# 	if (k > 0): sys.stderr.write(',');
		# 	if (pileup[i][k].op == 'I'): sys.stderr.write('[%s, %s, %d]' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
		# 	else: sys.stderr.write('(%s, %s, %d)' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id));
		# sys.stderr.write('\n');

			# except:
			# 	print j, len(pileup[i]), cov_depth;
			# 	sys.stderr.write('%d\t' % i);
			# 	for k in xrange(0, len(pileup[i])):
			# 		if (k > 0): sys.stderr.write(',');
			# 		if (pileup[i][k].op == 'I'): sys.stderr.write('[%s]' % (pileup[i][k].seq));
			# 		else: sys.stderr.write('%s' % (pileup[i][k].seq));
			# 	sys.stderr.write('\n');

	fp = open(out_msa_file, 'w');
	# for i in xrange(0, len(pileup)):
	# 	fp.write('%d\t' % i);
	# 	for k in xrange(0, len(pileup[i])):
	# 		if (k > 0): fp.write(',');
	# 		# if (pileup[i][k].op == 'I'): fp.write('[%s, %s, %d, %d]' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id, pileup[i][k].qid));
	# 		# else: fp.write('(%s, %s, %d, %d)' % (pileup[i][k].seq, pileup[i][k].op, pileup[i][k].id, pileup[i][k].qid));

	# 		# if (pileup[i][k].op == 'I'): fp.write('[%s, %d]' % (pileup[i][k].seq, cov_depth-k));
	# 		# else: fp.write('(%s, %d)' % (pileup[i][k].seq, cov_depth-k));

	# 		if (pileup[i][k].op == 'I'): fp.write('[%s]' % (pileup[i][k].seq));
	# 		else: fp.write('%s' % (pileup[i][k].seq));
	# 	fp.write('\n');
	for i in xrange(0, cov_depth):
		fp.write('>%s\n%s\n' % (msa_headers[i], msa_matrix[i]));



	fp.close();
	return 0;



if __name__ == "__main__":
	pass;
