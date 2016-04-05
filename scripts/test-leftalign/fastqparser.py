#! /usr/bin/python

# Script written by Ivan Sovic.
# Copyright Ivan Sovic, 2014. All rights reserved.
#               www.sovic.org
#
# Module for parsing FASTA/FASTQ files.

import sys;
import numpy as np;

def peek(fp, num_chars):
	data = fp.read(num_chars);
	if len(data) == 0:
		return '';
	fp.seek(num_chars * -1, 1);
	return data;

# Returns a single read from the given FASTA/FASTQ file.
# Parameter header contains only the header of the read.
# Parameter lines contains all lines of the read, which include:
# - header
# - seq
# - '+' if FASTQ
# - quals if FASTQ
# Parameter lines is an array of strings, each for one component.
# Please note that multiline FASTA/FASTQ entries (e.g. sequence line)
# will be truncated into one single line.
# def get_single_read_faulty(fp):
# 	lines = [];
	
# 	line = fp.readline();
# 	header = line.rstrip();
# 	header_leading_char = '';
# 	if (len(header) > 0):
# 		sequence_separator = header[0];
# 		header_leading_char = header[0];
# 		header = header[1:];			# Strip the '>' or '@' sign from the beginning.
# 	else:
# 		return ['', []];
	
# 	next_char = peek(fp, 1);
	
# 	line_string = '';
# 	lines.append(header_leading_char + header);
	
# 	num_lines = 1;
# 	#while len(next_char) > 0 and next_char != sequence_separator or (next_char == '@' and num_lines < 4):
# 	while (len(next_char) > 0 and (next_char != sequence_separator or (next_char == '@' and num_lines < 4))):
# 		line = fp.readline();
# 		if (line.rstrip() == '+' or line.rstrip() == ('+' + header)):
# 		#if (line.rstrip()[0] == '+'):
# 			lines.append(line_string);
# 			lines.append(line.rstrip());
# 			line_string = '';
# 		else:
# 			line_string += line.rstrip();
# 		next_char = peek(fp, 1);
# 		num_lines += 1;
		
# 	lines.append(line_string);
	
# 	return [header, lines];

def get_single_read(fp):
	lines = [];
	
	line = fp.readline();
	header = line.rstrip();
	header_leading_char = '';
	if (len(header) > 0):
		sequence_separator = header[0];
		header_leading_char = header[0];
		header = header[1:];			# Strip the '>' or '@' sign from the beginning.
	else:
		return ['', []];
	
	if (header_leading_char == '>'):
		num_expected_lines = 2;
	elif (header_leading_char == '@'):
		num_expected_lines = 4;
	else:
		sys.stderr.write('ERROR: Sequence is not in FASTA nor FASTQ format! Leading character in the header is "%s". Exiting.\n' % (header_leading_char));
		exit(1);

	next_char = peek(fp, 1);
	num_loaded_lines = 1;
	
	line_string = '';
	lines.append(header_leading_char + header);

	num_lines = 1;
	while (num_loaded_lines < num_expected_lines and len(next_char) > 0):
		line = fp.readline().strip();
		next_char = peek(fp, 1);

		if (header_leading_char == '@'):
			if ((num_loaded_lines == 1 and next_char == '+') or
				(num_loaded_lines == 2 and line[0] == '+') or
				(num_loaded_lines == 3 and next_char == sequence_separator)):
					line_string += line.rstrip();
					lines.append(line_string);
					num_loaded_lines += 1;
					line_string = '';
			elif (num_loaded_lines == 2 and line[0] != '+'):
				sys.stderr.write('ERROR: File is not in FASTQ format! Separator between seq and qual is not "+" but "%s". Exiting.\n' % (line[0]));
				exit(1);
			else:
				line_string += line.rstrip();

		elif (header_leading_char == '>'):
			if (num_loaded_lines == 1 and next_char == sequence_separator):
				line_string += line.rstrip();
				lines.append(line_string);
				num_loaded_lines = 2;
				line_string = '';
			else:
				line_string += line.rstrip();

		num_lines += 1;
	if (line_string != ''):
		lines.append(line_string);
		num_loaded_lines += 1;

	return [header, lines];


def read_fastq(fastq_path):
	headers = [];
	seqs = [];
	quals = [];
	
	fp_in = None;

	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % fastq_path);
		return;
	
	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		qual = '';
		if (len(read) == 4):
			qual = read[3];
		headers.append(header);
		seqs.append(seq);
		quals.append(qual);
		
	fp_in.close();
	
	return [headers, seqs, quals];

def convert_to_fasta(fastq_path, out_fasta_path):
	headers = [];
	seqs = [];
	quals = [];
	
	fp_in = None;
	fp_out = None;
	
	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		return;
	
	try:
		fp_out = open(out_fasta_path, 'w');
	except IOError:
		print 'ERROR: Could not open file "%s" for writing!' % out_fasta_path;
		fp_in.close();
		return;
	
	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		
		fp_out.write('>' + header + '\n');
		fp_out.write(seq + '\n');

	fp_in.close();
	fp_out.close();

def count_seq_length(fastq_path):
	fp_in = None;
	
	# total_seq_len = 0;
	num_seqs = 0;
	average_seq_len = 0;
	
	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		return;
	
	all_read_lengths = [];
	while True:
		[header, read] = get_single_read(fp_in);
		
		if (len(header) == 0):
			break;
		
		seq = read[1];
		# total_seq_len += len(seq);
		all_read_lengths.append(len(seq));
		num_seqs += 1;

	fp_in.close();
	
	# average_seq_len = (float(total_seq_len) / float(num_seqs)) if (num_seqs != 0) else 0.0;

	# read_length_stats = [np.mean(all_read_lengths), np.std(all_read_lengths), np.median(all_read_lengths), np.min(all_read_lengths), np.max(all_read_lengths)];

	if (len(all_read_lengths) > 0):
		total_seq_len = sum(all_read_lengths);
		average_seq_len = float(np.mean(all_read_lengths));
		min_seq_len = np.min(all_read_lengths);
		max_seq_len = np.max(all_read_lengths);
		std_seq_len = float(np.std(all_read_lengths));
		median_seq_len = float(np.median(all_read_lengths));
	else:
		total_seq_len = 0;
		average_seq_len = 0.0;
		min_seq_len = 0;
		max_seq_len = 0;
		std_seq_len = 0.0;
		median_seq_len = 0;
		sys.stderr.write('Warning: No sequences were loaded from file "%s" (function: fastqparser.count_seq_length).\n' % (fastq_path));
		sys.stderr.flush();

	ret_string = 'Number of sequences: %d\n' % num_seqs;
	ret_string += 'Total sequence length: %d\n' % total_seq_len;
	ret_string += 'Average sequence length: %.2f\n' % average_seq_len;
	ret_string += 'STD of sequence lengths: %.2f\n' % std_seq_len;
	ret_string += 'Median of sequence lengths: %.2f\n' % median_seq_len;
	ret_string += 'Min sequence length: %d\n' % min_seq_len;
	ret_string += 'Max sequence length: %d\n' % max_seq_len;
	
	return [ret_string, num_seqs, total_seq_len, average_seq_len, max_seq_len];

def get_headers_and_lengths(fastq_path):
	headers = [];
	lengths = [];
	fp_in = None;
	try:
		fp_in = open(fastq_path, 'r');
	except IOError:
		print 'ERROR: Could not open file "%s" for reading!' % fastq_path;
		exit(1);
	seq_id = 0;
	while True:
		[header, read] = get_single_read(fp_in);
		if (len(header) == 0):
			break;
		headers.append(header);
		lengths.append(len(read[1]));
		seq_id += 1;
	fp_in.close();
	return [headers, lengths];

def complement_base(base):
	if (base == 'A'):
		return 'T';
	if (base == 'C'):
		return 'G';
	if (base == 'T'):
		return 'A';
	if (base == 'G'):
		return 'C';
	return 'N';

def revcomp_seq(sequence):
	ret_seq = '';
	i = 0;
	while (i < len(sequence)):
		ret_seq += complement_base(sequence[len(sequence)-i-1]);
		i += 1;
	return ret_seq;

if __name__ == "__main__":
	pass;
