#! /usr/bin/python

# Copyright Ivan Sovic, 2015. www.sovic.org
#
# Module for parsing and processing of SAM files.

import os;
import re;
import sys;

CIGAR_OPERATIONS_ALL = ['M', 'I', 'D', 'N', 'S', 'H', 'P', '=', 'X'];
CIGAR_OPERATIONS_BASIC = ['M', 'I', 'D', 'S', 'H'];
CIGAR_OPERATIONS_EXTENDED = ['M', 'I', 'D', 'S', 'H', '=', 'X'];



class SAMLine:
	def __init__(self, line='', sam_basename=''):
		if line == '':
			self.Clear();
		else:
			self.ParseLine(line, sam_basename);

	def Clear(self):
		# (I) These parameters are parsed from the SAM file.
		self.qname = '';
		self.flag = 0;
		self.rname = '';
		self.pos = 0;
		self.mapq = 0;
		self.cigar = '';
		self.mrnm = '';
		self.mpos = 0;
		self.isize = 0;
		self.seq = '';
		self.qual = '';
		self.optional = {};
		self.original_line = '';
		self.evalue = -1.0;
		self.alignment_score = -1;
		self.edit_distance = -1;
		self.fp_filter = -1.0;
		self.optional_graphmap = {};

		# (II) These parameters are parsed from the parameters parsed from the SAM file.
		self.clipped_pos = 0;			# If the read has clipping operation at the beginning, then self.clipped_pos = self.pos - num_clipped_bases . This is to allow easy comparison to the reference.
		self.clip_op_front = '';
		self.clip_count_front= 0;
		self.clip_op_back = '';
		self.clip_count_back = 0;
		self.chosen_quality = 0;			# In SAM, quality can be specified either as mapping quality (mapq), or as alignment score (AS optional parameter). For instance, LAST outputs AS, and not mapq. If mapq == 255, then it will be checked if there is an AS entry in the given optional parameters. If AS is not present, then self.chose_quality will be equal to self.mapq, that is, 255.
		
		# (III) These parameters are modified during evaluation of mapping.
		self.evaluated = 0;
		self.is_correct_ref_and_orient = 0;
		self.is_duplicate = 0;
		self.is_best_of_duplicates = 0;
		self.actual_ref_reverse = 0;
		self.is_filtered_out = False;
		self.num_occurances_in_sam_file = 0;
		
		self.is_header_deformed = 0;
		self.actual_ref_header = '';
		self.trimmed_header = '';
		self.actual_ref_pos = 0;
		self.mapped_pos_with_shift = 0;
		self.min_distance = -1;

		self.num_correct_m_ops = 0;
		
		self.sam_basename = '';

		# For sanity checking - this is True if the SAM line contains all fields. It is false i.e. if file writing was interrupted before line was finished.
		self.line_fields_ok = False;

	def VerboseFormatLine(self):
		#line = 'qname = %s\tflag = %d\trname = %s\tpos = %d\tmapq = %d\tis_header_deformed = %d\tis_correct_ref_and_orient = %d\tis_duplicate = %d\tis_best_of_duplicates = %d\tactual_ref_pos = %d\tmapped_pos_with_shift = %d\tmin_distance = %d\tsecondary = %s\tactual_ref_reverse = %d\treverse = %d\tpaired = %s\t' % (
			#self.qname, self.flag, self.rname, self.pos, self. mapq, self.is_header_deformed, int(self.is_correct_ref_and_orient), int(self.is_duplicate), int(self.is_best_of_duplicates), self.actual_ref_pos, self.mapped_pos_with_shift, self.min_distance, str(self.IsSecondary()), self.actual_ref_reverse, self.IsReverse(), self.IsPaired());
		line = 'qname = %s\tactual_ref_pos = %d\tmapped_pos_with_shift = %d\tmin_distance = %d\tmapq = %d\tactual_ref_reverse = %d\treverse = %d\tis_correct_ref_and_orient = %d\tis_header_deformed = %d\tflag = %d\trname = %s\tpos = %d\tis_duplicate = %d\tis_filtered_out = %d\tis_best_of_duplicates = %d\tsecondary = %s\tpaired = %s\t' % (
			self.qname, self.actual_ref_pos, self.mapped_pos_with_shift, self.min_distance, self.mapq, self.actual_ref_reverse, self.IsReverse(), int(self.is_correct_ref_and_orient), self.is_header_deformed, self.flag, self.rname, self.pos, int(self.is_duplicate), int(self.is_filtered_out), int(self.is_best_of_duplicates), str(self.IsSecondary()), self.IsPaired());
		
		return line;

	def ParseLine(self, line, sam_basename=''):
		split_line = line.split('\t');
		
		if len(split_line) < 11:
			sys.stderr.write('ERROR: Line does not contain all mandatory SAM fields!\n');
			sys.stderr.write('Line: "%s"\n' % line);
			self.Clear();
			self.line_fields_ok = False;
			return;
		
		self.original_line = line;

		self.qname = split_line[0];
		self.flag = int(split_line[1]);
		self.rname = split_line[2];
		self.pos = int(split_line[3]);
		self.mapq = int(split_line[4]);
		self.cigar = split_line[5];
		self.mrnm = split_line[6];
		self.mpos = int(split_line[7]);
		self.isize = int(split_line[8]);
		self.seq = split_line[9];
		self.qual = split_line[10];
		
		### Process optional parameters.
		self.optional = {};
		i = 11;
		while i < len(split_line):
			split_optional = split_line[i].split(':');
			if (len(split_optional) < 2):
				i += 1;
				continue;
			self.optional[split_optional[0].strip()] = split_optional[-1].strip();		# Example of an optional parameter: AS:i:717 .
			i += 1;
		

		self.alignment_score = -1 if (('AS' in self.optional) == False) else int(self.optional['AS'])
		self.edit_distance = -1 if (('NM' in self.optional) == False) else int(self.optional['NM'])
		self.fp_filter = -1.0 if (('X5' in self.optional) == False) else float(self.optional['X5'])
		if (self.fp_filter == -1.0):
			self.fp_filter = -1.0 if (('ZF' in self.optional) == False) else float(self.optional['ZF'])

		self.optional_graphmap = {};
		if ('X3' in self.optional):

			params = self.optional['X3'];
			# print params;
			# print '';
			m = re.search(r'_supp\[(.*?)\]', params);			self.optional_graphmap['supp'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_AS\[(.*?)\]', params);				self.optional_graphmap['AS'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_AS_std\[(.*?)\]', params);			self.optional_graphmap['AS_std'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_AS_cov_bases\[(.*?)\]', params);	self.optional_graphmap['AS_cov_bases'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_AS_read_len\[(.*?)\]', params);	self.optional_graphmap['AS_read_len'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_AS_query_len\[(.*?)\]', params);	self.optional_graphmap['AS_query_len'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_num_kmers=(.*?)_', params);	self.optional_graphmap['num_kmers'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_cov_bases=(.*?)_', params);	self.optional_graphmap['cov_bases'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'lcs_length=(.*?)_', params);	self.optional_graphmap['lcs_length'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_match_rate=(.*?)_', params);	self.optional_graphmap['match_rate'] = float(m.group(1)) if (m and m.groups) else 1.0;
			m = re.search(r'_mismatch_rate=(.*?)\t', params);	self.optional_graphmap['mismatch_rate'] = float(m.group(1)) if (m and m.groups) else 0.0;
			m = re.search(r'_num_eq_ops=(.*?)_', params);	self.optional_graphmap['num_eq_ops'] = float(m.group(1)) if (m and m.groups) else 0.0;
			m = re.search(r'_num_x_ops=(.*?)_', params);	self.optional_graphmap['num_x_ops'] = float(m.group(1)) if (m and m.groups) else 0.0;
			m = re.search(r'_num_i_ops=(.*?)_', params);	self.optional_graphmap['num_i_ops'] = float(m.group(1)) if (m and m.groups) else 0.0;
			m = re.search(r'_num_d_ops=(.*?)_', params);	self.optional_graphmap['num_d_ops'] = float(m.group(1)) if (m and m.groups) else 0.0;

		self.chosen_quality = self.mapq + 0;
		# if (self.chosen_quality == 255):
		if (self.alignment_score != -1):
			self.chosen_quality = self.alignment_score;

		self.evalue = -1.0;
		if ('Z1' in self.optional):
			try:
				self.evalue = float(self.optional['Z1'].split('x')[0]);
			except:
				self.evalue = 5000000000.0;
		elif ('ZE' in self.optional):
			try:
				self.evalue = float(self.optional['ZE'].split('x')[0]);
			except:
				self.evalue = 5000000000.0;
		
		### Initialize evaluation values.
		self.evaluated = 0;
		self.is_correct_ref_and_orient = 0;
		self.is_duplicate = 0;
		self.is_best_of_duplicates = 0;
		self.actual_ref_reverse = 0;
		self.is_filtered_out = False;
		
		self.is_header_deformed = 0;
		self.actual_ref_header = '';
		self.trimmed_header = '';
		self.actual_ref_pos = 0;
		self.mapped_pos_with_shift = 0;		# This is assigned only after evaluation of the mapping position. It is equal to self.clipped_pos if self.qname is found in the reference SAM file, otherwise it is equal to 0.
		self.min_distance = -1;
		
		self.clipped_pos = self.pos;		# The position of the read subtracted by the amount of clipped bases at the begining of the read.
		self.clip_op_front = '';
		self.clip_count_front= 0;
		self.clip_op_back = '';
		self.clip_count_back = 0;
		self.num_occurances_in_sam_file = 1;

		self.num_correct_m_ops = 0;
		
		self.sam_basename = sam_basename;
		
		if (len(self.cigar) > 0 and self.cigar != '*'):
			m_front = re.match("^([\d]+)([SH])", self.cigar);
			m_back = re.match(".*[MIDNSHP=X]([\d]+)([SH])$", self.cigar);
			
			if (m_front):
				cigarcount = int(m_front.group(1));
				cigarop = m_front.group(2);
				self.clipped_pos -= cigarcount;
				self.clip_op_front = cigarop;
				self.clip_count_front = cigarcount;
				#print '(front) cigarcount = %d, cigarop = %s, sam_basename = %s, qname = %s' % (cigarcount, cigarop, self.sam_basename, self.qname);
			if (m_back):
				cigarcount = int(m_back.group(1));
				cigarop = m_back.group(2);
				self.clip_op_back = cigarop;
				self.clip_count_back = cigarcount;
				#print '(back) cigarcount = %d, cigarop = %s, sam_basename = %s, qname = %s' % (cigarcount, cigarop, self.sam_basename, self.qname);

		#if (self.pos != self.clipped_pos):
			#print 'Clipped position: %d (original: %d, clip_op = %s, clip_count = %d)' % (self.clipped_pos, self.pos, self.clip_op, self.clip_count);
		self.line_fields_ok = True;
	
	def Verbose(self):
		print 'qname = %s' % self.qname;
		print 'flag = %s' % self.flag;
		print 'rname = %s' % self.rname;
		print 'pos = %s' % self.pos;
		print 'mapq = %s' % self.mapq;
		print 'cigar = %s' % self.cigar;
		print 'mrnm = %s' % self.mrnm;
		print 'mpos = %s' % self.mpos;
		print 'isize = %s' % self.isize;
		print 'seq = %s' % self.seq;
		print 'qual = %s' % self.qual;
		print '(evaluated = %d)' % self.evaluated;
		print '(min_distance = %d)' % self.min_distance;
		print '(is_correct_ref_and_orient = %d)' % self.is_correct_ref_and_orient;
		print '(is_duplicate = %d)' % self.is_duplicate;
		print '(is_filtered_out = %d)' % self.is_filtered_out;
		print '(is_best_of_duplicates = %d)' % self.is_best_of_duplicates;
		print '(clipped_pos = %d)' % self.clipped_pos;
	
	def FormatAccuracy(self):
		line = '';
		
		#query.min_distance = ret_min_distance;
		#query.actual_ref_pos = ret_reference_pos;
		#query.actual_ref_reverse = ret_ref_reverse;
		#query.actual_ref_header = ret_ref_header;
		#query.mapped_pos_with_shift = ret_mapped_pos;
		
		line += 'distance = %d\t' % self.min_distance;
		line += 'header_hit = %s\t' % (str(self.rname == self.actual_ref_header));
		line += 'reverse_hit = %s\t' % (str(self.actual_ref_reverse == self.IsReverse()));
		line += 'qname = %s\t' % self.qname;
		line += 'mapped = %s\t' % str(self.IsMapped());
		#line = 'header_def = %s\t'
		line += 'ref_header = %s\t' % self.actual_ref_header;
		line += 'map_header = %s\t' % self.rname;
		line += 'ref_pos = %d\t' % self.actual_ref_pos;
		line += 'map_pos_clipped = %d\t' % self.clipped_pos;
		line += 'ref_reverse = %s\t' % str(self.actual_ref_reverse);
		line += 'map_reverse = %s\t' % str(self.IsReverse());
		#line += 'duplicate = %s\t' % (str(self.is_duplicate != 0));
		line += 'cigar_start = %s\t' % (self.cigar[0:5]);
		line += 'cigar_end = %s\t' % (self.cigar[(len(self.cigar)-5):]);
		# line += 'chosen_quality = %d\t' % (self.chosen_quality);
		line += 'mapq = %d\t' % (self.mapq);
		line += 'AS = %d\t' % (0 if (('AS' in self.optional) == False) else int(self.optional['AS']));
		line += 'EValue = %s\t' % self.evalue;
		line += 'num_occ = %d\t' % (self.num_occurances_in_sam_file);
		line += 'fp_filter = %s\t' % self.fp_filter;
		line += 'NM = %d\t' % self.edit_distance;
		line += 'len = %d\t' % (len(self.seq) - (self.clip_count_front if (self.clip_op_front == 'S') else 0) - (self.clip_count_back if (self.clip_op_back == 'S') else 0));
		line += 'correct_M = %d' % (self.num_correct_m_ops);
		
		return line;
	
		
	
	# Checks the SAM flag to see if the read is paired end.
	def IsPaired(self):
		return (True if (self.flag & 0x01) else False);
	
	# Checks the SAM flag to see if the alignment is reported as mapped.
	def IsMapped(self):
		return (False if (self.flag & 0x04) else True);
		#return (True if (self.flag & 0x04) else False);
	
	# Checks the SAM flag to see if the SEQ was reversed.
	def IsReverse(self):
		return (True if (self.flag & 0x10) else False);
	
	# Checks the SAM flag to see if this is a secondary alignment.
	def IsSecondary(self):
		return (True if (self.flag & 0x100) else False);

	# Splits the CIGAR string into individual operations, in the
	# same format as the original CIGAR string is in.
	# The result is returned as an array of a 2-element array, e.g.
	# [[12, 'M'], [3, 'D'], [4, 'M']].
	def SplitCigar(self):
		i = 0;
		cigarcount_string = '';
		cigar_operations = [];
		if (self.IsMapped() == False):
			return cigar_operations;
		while i < len(self.cigar):
			if (self.cigar[i] in CIGAR_OPERATIONS_EXTENDED):
				cigar_operations.append([int(cigarcount_string), self.cigar[i]]);
				cigarcount_string = '';
			else:
				cigarcount_string += self.cigar[i];
			i += 1;
		if (cigarcount_string != ''):
			sys.stderr.write('ERROR: Faulty CIGAR string!\n');
			sys.stderr.write('cigarcount_string: %s\n' % cigarcount_string);
			sys.stderr.write('i = %d\n' % i);
			sys.stderr.write('%s\n' % self.original_line);
			
			cigar_operations = [];
		return cigar_operations;
	
	# Splits the CIGAR string into individual operations. Unlike SplitCigar,
	# this function also converts the extended cigar format to the basic cigar
	# format. This includes counting the number of successive M, = and X operations.
	# The result is returned as an array of a 2-element array, e.g.
	# [[12, 'M'], [3, 'D'], [4, 'M']].
	def SplitCigarInBasicFormat(self):
		i = 0;
		cigarcount_string = '';
		cigar_operations = [];
		if (self.IsMapped() == False):
			return cigar_operations;
		while i < len(self.cigar):
			if (self.cigar[i] in CIGAR_OPERATIONS_EXTENDED):
				# Check if it is a match/mismatch operation and convert to M:
				if (self.cigar[i] in 'M=X'):
					# If it is a match/m[headers, sam_lines]ismatch and the list is empty, simply add it.
					if (len(cigar_operations) == 0):
						cigar_operations.append([int(cigarcount_string), 'M']);
					# If the list is not empty, check if the previous operation was an 'M'.
					elif (cigar_operations[-1][1] == 'M'):
						cigar_operations[-1][0] += int(cigarcount_string);
					# If the previous operation was not an M, then simply add it again.
					else:
						cigar_operations.append([int(cigarcount_string), 'M']);
				else:
					cigar_operations.append([int(cigarcount_string), self.cigar[i]]);
				cigarcount_string = '';
			else:
				cigarcount_string += self.cigar[i];
			i += 1;
		if (cigarcount_string != ''):
			sys.stderr.write('ERROR: Faulty CIGAR string!\n');
			sys.stderr.write('cigarcount_string: %s\n' % cigarcount_string);
			sys.stderr.write('i = %d\n' % i);
			sys.stderr.write('%s\n' % self.original_line);
			
			cigar_operations = [];
		return cigar_operations;
	
	# Sums the counts of M/I/S/=/X operations in the CIGAR string to determine the
	# length of SEQ. This is used for testing, debugging and sanity checking. If the
	# SAM line has been properly formatted, the result of this function should be
	# the same as len(self.seq).
	# Although this reflects on the seq length, if seq is hard clipped, this does
	# not report the length of the entire sequence (read) if it were unclipped.
	# For this purpose, check CalcReadLengthFromCigar function.
	# Edit: 14.04.2015. Added new parameter: switch_ins_and_dels.
	# It seems that PBSim reports the variations that it induced in the reference
	# instead of reporting the variations from the read's perspective. In other words,
	# when a PBSim's CIGAR has an insertion, that is an insertion to the reference and
	# not the read (vice versa goes for the deletions). Normally, CIGAR is reported
	# from the read's point of view.
	def CalcAlignmentLengthFromCigar(self, switch_ins_and_dels=False):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			# if ((switch_ins_and_dels == False and (cigar_op in 'MIS=X')) or
			# 	(switch_ins_and_dels == True and (cigar_op in 'MDS=X'))):
			if ((switch_ins_and_dels == False and (cigar_op in 'MDS=X')) or
				(switch_ins_and_dels == True and (cigar_op in 'MIS=X'))):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	# Sums the counts of M/I/S/H/=/X operations in the CIGAR string to determine the
	# length of the entire read. Unlike CalcAlignmentLengthFromCigar, the value returned
	# by this function may be greater or equal to the length of seq, because seq can be
	# hard clipped.
	def CalcReadLengthFromCigar(self, switch_ins_and_dels=False):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			# if (cigar_op in 'MISH=X'):
			if ((switch_ins_and_dels == False and (cigar_op in 'MISH=X')) or
				(switch_ins_and_dels == True and (cigar_op in 'MDSH=X'))):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	# Sums the counts of M/D/=/X operations in the CIGAR string to determine the
	# length of the reference covered by the read.
	def CalcReferenceLengthFromCigar(self, switch_ins_and_dels=False):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			# if (cigar_op in 'MD=X'):
			if ((switch_ins_and_dels == False and (cigar_op in 'MD=X')) or
				(switch_ins_and_dels == True and (cigar_op in 'MI=X'))):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	def CountCIGAREvents(self):
		split_cigar = self.SplitCigar();
		alignment_length = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigar_op in 'MD=X'):
				alignment_length += cigar_count;
			i += 1;
		return alignment_length;

	# Counts all CIGAR ops separately and returns them as a dics.
	def CountAlignmentOps(self):
		split_cigar = self.SplitCigar();
		op_counts = {'=': 0, 'X': 0, 'I': 0, 'D': 0};
		for val in split_cigar:
			[cig_count, cig_op] = val;
			try:
				op_counts[cig_op] += cig_count;
			except:
				op_counts[cig_op] = cig_count;
		return op_counts;

	# Given a base position within the read, find its position on the reference by inspecting CIGAR.
	def FindBasePositionOnReference(self, base_position_in_read):
		cigar_pos_list = self.CalcCigarStartingPositions(True);
		i = 0;
		while (i < len(cigar_pos_list)):
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_list[i];
			if (cigar_op in 'M=X'):
				if (pos_on_read == base_position_in_read):
					return pos_on_reference;
				elif (pos_on_read < base_position_in_read and (pos_on_read + cigar_count) > base_position_in_read):
					return (pos_on_reference + (base_position_in_read - pos_on_read));
			i += 1;
		return -1;

	# Given a base coordinate on the reference, find its position on the read.
	def FindBasePositionOnRead(self, base_position_in_reference):
		cigar_pos_list = self.CalcCigarStartingPositions(True);
		i = 0;
		while (i < len(cigar_pos_list)):
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_list[i];
			if (cigar_op in 'M=XD'):
				if (pos_on_reference == base_position_in_reference):
					return pos_on_read;
				elif (pos_on_reference < base_position_in_reference and (pos_on_reference + cigar_count) > base_position_in_reference):
					return (pos_on_read + (base_position_in_reference - pos_on_reference));
			i += 1;
		return -1;

	### Given a start and end position on the read, this function extracts all CIGAR events for bases inbetween. End position is inclusive.
	def GetCigarBetweenBases(self, start_pos, end_pos):
		cigar_pos_list = self.CalcCigarStartingPositions(True);
		
		start_event = -1;
		end_event = -1;

		for i in xrange(0, len(cigar_pos_list)):
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_list[i];
			if (cigar_op in 'M=XI'):
				if (pos_on_read == start_pos):
					start_event = i;
					break;
				elif (pos_on_read < start_pos and (pos_on_read + cigar_count) > start_pos):
					start_event = i;
					break;

		if (start_event == -1):
			return [];

		for i in xrange(start_event, len(cigar_pos_list)):
			[cigar_count, cigar_op, pos_on_reference, pos_on_read] = cigar_pos_list[i];
			if (cigar_op in 'M=XI'):
				if (pos_on_read == end_pos):
					end_event = i;
					break;
				elif (pos_on_read < end_pos and (pos_on_read + cigar_count) > end_pos):
					end_event = i;
					break;

		if (end_event == -1):
			return [];

		return cigar_pos_list[start_event:(end_event+1)];

	# Splits the CIGAR string into individual operations, and determines their starting positions
	# on the reference. I.e. this function identifies the positions of matches, mismatches and deletions.
	# This is used, for example, in comparison of alignment operations to the reference SAM file.
	# If separate_matches_in_individual_bases == True, then 'M', '=' and 'X' operations will
	# be split into individual operations of length 1 (i.e. if the input CIGAR operation was
	# 3M, this would generate 1M1M1M. This is a nice hack for quick checking coordinates of
	# individual matching/mismatching bases.)
	# Function returns an array of arrays:
	# [[cigar_count1, cigar_op1, pos_on_reference1, pos_on_read1], [[cigar_count2, cigar_op2, pos_on_reference2, pos_on_read2], ...]
	# Example usage:
	# import utility_sam;
	# [headers, sam_lines] = utility_sam.LoadSAM(sam_path);
	# for sam_line in sam_lines:
	#	[cigar_count, cigar_op, pos_on_ref, pos_on_read] = sam_line.CalcCigarStartingPositions();
	
	def CalcCigarStartingPositions(self, separate_matches_in_individual_bases=False, switch_ins_and_dels=False):
		#cigar_list = self.SplitCigar();
		cigar_list = self.SplitCigarInBasicFormat();
		cigar_pos_list = [];
		# 08.03.2015. Added the -1 to transform the reference coordinates from 1-based to 0-based.
		pos_on_reference = self.clipped_pos + 0 - 1;
		pos_on_read = 0;
		
		i = 0;
		while (i < len(cigar_list)):
			cigar_count = cigar_list[i][0];
			cigar_op = cigar_list[i][1];
			
			if (separate_matches_in_individual_bases == False):
				# Edit 10.04.2015. I saw that the pos_on_read value was missing here, so I added it too:
				cigar_pos_list.append([cigar_count + 0, cigar_op + '', pos_on_reference + 0, pos_on_read + 0]);
				#if (cigar_op in 'MDSH=X'):
					#pos_on_reference += cigar_count;
				
				# S and H are also used to walk down the reference, because pos_on_reference was initialized
				# with self.clipped_pos, which is calculated from the SAM pos field by subtracting the number
				# of clipped bases. Otherwise, S and H should not be used to increase pos_on_reference.
				if (cigar_op in 'MSH=X'):
					pos_on_reference += cigar_count;
					pos_on_read += cigar_count;
				elif ((switch_ins_and_dels == False and cigar_op == 'D') or  (switch_ins_and_dels == True and cigar_op == 'I')):
					pos_on_reference += cigar_count;
				elif ((switch_ins_and_dels == False and cigar_op == 'I') or  (switch_ins_and_dels == True and cigar_op == 'D')):
					pos_on_read += cigar_count;
					
			else:
				if ((cigar_op in 'M=X') == False):
					cigar_pos_list.append([cigar_count + 0, cigar_op + '', pos_on_reference + 0, pos_on_read + 0]);
					if (cigar_op in 'SH'):
						pos_on_reference += cigar_count;
						pos_on_read += cigar_count;
					elif ((switch_ins_and_dels == False and cigar_op == 'D') or  (switch_ins_and_dels == True and cigar_op == 'I')):
						pos_on_reference += cigar_count;
					elif ((switch_ins_and_dels == False and cigar_op == 'I') or  (switch_ins_and_dels == True and cigar_op == 'D')):
						pos_on_read += cigar_count;
				else:
					j = 0;
					while (j < cigar_count):
						cigar_pos_list.append([1, cigar_op, pos_on_reference + 0, pos_on_read + 0]);
						pos_on_reference += 1;
						pos_on_read += 1;
						j += 1;
			
			i += 1;
		
		return cigar_pos_list;

	def CalcNumMappedBases(self):
		num_mapped_bases = len(self.seq);
		if (self.clip_op_front == 'S'):
			num_mapped_bases -= self.clip_count_front;
		if (self.clip_op_back == 'S'):
			num_mapped_bases -= self.clip_count_back;
		return num_mapped_bases;

	def IsAlignmentSane(self):
		if (self.IsMapped() == False):
			return True;

		if (self.line_fields_ok == False):
			return False;

		seq_len = len(self.seq);

		split_cigar = self.SplitCigar();
		alignment_length_on_seq = 0;
		i = 0;
		while i < len(split_cigar):
			cigar_count = split_cigar[i][0];
			cigar_op = split_cigar[i][1];
			# From the SAM format specification:
			#     Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.
			if (cigar_op in 'MIS=X'):
				alignment_length_on_seq += cigar_count;

			if ((cigar_op == 'D' or cigar_op == 'I') and cigar_count > 100):
				return False;
			if ((i + 1) < len(split_cigar) and (cigar_op in 'DI') and (split_cigar[i + 1][1] in 'DISH')):
				return False;

			i += 1;

		if (seq_len != alignment_length_on_seq):
			return False;

		return True;

	def GetNumCigarOps(self):
		if (self.IsMapped() == False):
			return 0;

		if (self.line_fields_ok == False):
			return 0;

		split_cigar = self.SplitCigar();
		return len(split_cigar);

	def CountErroneousWindows(self, window_len, min_error_rate):
		num_windows = 0;

		cigars = self.CalcCigarStartingPositions(separate_matches_in_individual_bases=True);

		# for cigar in cigars:
		# 	[cigar_count, cigar_op, pos_on_ref, pos_on_read] = cigar;
		# 	print cigar;

		seq_len = len(self.seq);
		window_start_base = 0;
		window_end_base = 0;
		window_cig_id = 0;

		# while (cig_op[window_start_base])
		real_bases = [];
		real_bases_deletions = [];
		real_bases_insertions = [];
		i = 0;
		while (i < seq_len):
			real_bases.append('');
			real_bases_deletions.append(0);
			real_bases_insertions.append(0);
			i += 1;
		# real_bases = [['', 0]]*seq_len;
		# print real_bases;

		num_hardclipped_bases = 0;
		for cigar in cigars:
			[cigar_count, cigar_op, pos_on_ref, pos_on_read] = cigar;
			if (cigar_op in 'SIM=X'):
				i = 0;
				while (i < cigar_count):
					try:
						real_bases[pos_on_read + i - num_hardclipped_bases] = cigar_op;
					except Exception, e:
						sys.stderr.write(str(e) + '\n');
						print cigar;
						print seq_len, pos_on_read, num_hardclipped_bases, i, len(cigars);
						# print self.original_line;
						print self.qname;
						exit(1);
					if (cigar_op == 'I'):
						real_bases_insertions[pos_on_read + i - num_hardclipped_bases] = 1;
					i += 1;
			elif (cigar_op in 'D'):
				real_bases_deletions[pos_on_read - num_hardclipped_bases] = cigar_count;
			elif (cigar_op in 'H'):
				### Hardclipped bases would cause a problem with the position of the bases in the cigar operation list.
				### The pos_on_read is with the respect to the full read sequence, including the hardclipped bases which are not
				### present in the self.seq string.
				num_hardclipped_bases += cigar_count;
			else:	### Whatever else may occur.
				pass;

		# i = 0;
		# while (i < len(real_bases)):
		# 	print i, real_bases[i], real_bases_insertions[i], real_bases_deletions[i];
		# 	i += 1;
		# print '';
		# print seq_len;

		i = 0;
		while (i < seq_len):
			if ((real_bases[i] in 'SH') == False):
				break;
			i += 1;
		nonclip_start = i;
		i = seq_len - 1;
		while (i >= 0):
			if ((real_bases[i] in 'SH') == False):
				break;
			i -= 1;
		nonclip_end = i;

		all_window_ratios = [];
		all_windows_over_threshold = [];
		window_start = nonclip_start;
		window_end = window_start + window_len - 1;
		while (window_end < nonclip_end):
			window_ins = sum(real_bases_insertions[window_start:(window_end+1)]);
			window_dels = sum(real_bases_deletions[window_start:(window_end+1)]);
			window_errors = window_ins + window_dels;
			window_ratio = float(window_errors) / float(window_len);
			if (window_ratio > min_error_rate):
				all_windows_over_threshold.append(window_ratio);
				# print window_ratio, min_error_rate;
								# print 'Tu sam 1!\n';
			all_window_ratios.append(window_ratio);
			window_start += 1;
			window_end += 1;


		num_over_threshold = len(all_windows_over_threshold);
		num_windows = len(all_window_ratios);
		err_window_rate = (float(num_over_threshold) / float(num_windows)) if (num_windows > 0) else (-1.0);

		# print seq_len;
		# print len(real_bases);
		# print self.CalcReadLengthFromCigar();
		# print real_bases;
		# exit(1);

		# print window_start;
		# print window_end;
		# print window_errors;
		# print window_len;
		# print window_ratio;

		# print all_window_ratios;
		# print  '';
		# print all_windows_over_threshold;
		# print len(all_windows_over_threshold);
		# print 'threshold = %.2f' % (min_error_rate);
		# print 'num_over_threshold = %d' % (num_over_threshold);
		# print 'num_windows = %d' % (num_windows);
		# print 'ratio = %.2f' % (float(num_over_threshold) / float(num_windows));

		# # while ((window_start_base + window_len) < seq_len):
		# # 	if ()
		# while (True):
		# 	if (cig_op[window_cig_id] in 'SH'):
		# 		window_start_base += 1;
		# 		continue;

		# 	num_windows += 1;



		# 	# if (cig_op[winodw_cig_id] == 'H'):
		# 	# 	continue;

		# 	if ((window_start_base + window_len) >= seq_len):
		# 		break;
		# 	break;

		# if (self.clip_op_back == 'S'):
		# 	seq_len -= self.clip_count_back;

		# window_start = (self.clip_count_front) if (self.clip_op_front == 'S') else (0);
		# while (window_start < len(seq_len)):

		# 	window_start += 1;

		return [err_window_rate, num_over_threshold, num_windows];
		




def LoadSAM(sam_path, verbose=False):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_path);
		return [[], []];
	
	headers = [];
	sam_lines = [];
	
	i = 0;
	for line in fp_reference:
		i += 1;
		if (verbose == True):
			sys.stdout.write('\rParsing SAM line %d...' % (i));
		line = line.strip();
		if len(line) == 0:
			continue;
		if line[0] == '@':
			headers.append(line);
			continue;
		sam_line = SAMLine(line);
		sam_lines.append(sam_line);
	fp_reference.close();
	
	if (verbose == True):
		sys.stdout.write('done!\n');
	
	return [headers, sam_lines];

def LoadOnlySAMHeaders(sam_path, verbose=False):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_path);
		return [];
	
	headers = [];
	
	i = 0;
	for line in fp_reference:
		i += 1;
		if (verbose == True):
			sys.stdout.write('\rParsing SAM line %d...' % (i));
		line = line.strip();
		if len(line) == 0:
			continue;
		if line[0] == '@':
			headers.append(line);
			continue;
		else:
			break;
	fp_reference.close();
	
	if (verbose == True):
		sys.stdout.write('done!\n');
	
	return headers;

def HashSAMLines(sam_lines):
	ret = {};
	num_unique_lines = 0;
	num_lines = 0;
	for sam_line in sam_lines:
		qname = sam_line.qname;
		if (qname in ret):
				current_hash = ret[qname];

				should_be_counted = True;
				# We have to check if this sequence is actually a paired read of another sequence, or there is something funny in the data.
				for existing_sam in current_hash:
					if (sam_line.IsSecondary() == True) or (sam_line.IsSecondary() == existing_sam.IsSecondary() and sam_line.IsReverse() == existing_sam.IsReverse()):
						# This case should not be counted. It means that, either the alignment is marked as secondary which means there should be a primary alignment as well, or that there is more than one primary alignment, and the orientation is the same, which means that an aligner is trying to artificially boost-up their statistics.
						should_be_counted = False;
						break;
				if should_be_counted == True:
					num_unique_lines += 1;	# Count only unique sequences.
		else:
			# This is a new sequence, not hashed before. Create a new list and count the sequence.
			ret[qname] = [sam_line];
			if sam_line.IsSecondary() == False:
				num_unique_lines += 1;	# Count only unique sequences, but not secondary ones.

		num_lines += 1;	# Count only unique sequences.

	for key in ret.keys():
		ret[key].sort(reverse=True, key=lambda sam_line: sam_line.chosen_quality);
		# ret[key].sort(reverse=True, key=lambda sam_line: ((sam_line.IsSecondary() == False), sam_line.chosen_quality));

	return [ret, num_lines, num_unique_lines];



# Hashes the entries from a SAM file by their QNAME for faster lookup during comparison.
def HashSAM(sam_path):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_path);
		return [{}, 0];
	
	ret = {};
	
	num_unique_references = 0;
	num_references = 0;
	
	for line in fp_reference:
		line = line.strip();
		
		if len(line) == 0:
			continue;
		
		if line[0] == '@':
			continue;
		
		sam_line = SAMLine(line);
		
		modified_qname = sam_line.qname;
		#modified_qname = '/'.join(sam_line.qname.split('/')[:-1]);
		
		try:
			current_hash = ret[modified_qname];
			
			should_be_counted = True;
			
			# We have to check if this sequence is actually a paired read of another sequence, or there is something funny in the data.
			for existing_sam in current_hash:
				if (sam_line.IsSecondary() == True) or (sam_line.IsSecondary() == existing_sam.IsSecondary() and sam_line.IsReverse() == existing_sam.IsReverse()):
					# This case should not be counted. It means that, either the alignment is marked as secondary which means there should be a primary alignment as well, or that there is more than one primary alignment, and the orientation is the same, which means that an aligner is trying to artificially boost-up their statistics.
					should_be_counted = False;
					break;
				
			
			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print line;
				#for sam in ret[sam_line.qname]:
					#print 'is_secondary: %s, is_reverse = %s, num_unique_references: %d' % (str(sam.IsSecondary()), str(sam.IsReverse()), num_unique_references);
			
			if should_be_counted == True:
				num_unique_references += 1;	# Count only unique sequences.
				#print 'Tu sam 2!';
			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print '---';

			# At least one sequence with the same name has already been added.
			current_hash.append(sam_line);
		except:
			# This is a new sequence, unhashed before. Create a new list and count the sequence.
			ret[modified_qname] = [sam_line];
			if sam_line.IsSecondary() == False:
				num_unique_references += 1;	# Count only unique sequences, but not secondary ones.
				

			#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
				#print line;
				#if sam_line.IsSecondary() == False:
					#print 'Tu sam 1!';
					
				#print '---';

		#if sam_line.qname == 'gi|48994873|gb|U00096.2|-463960':
			#for sam in ret[sam_line.qname]:
				#print 'is_secondary: %s, is_reverse = %s, num_unique_references: %d' % (str(sam.IsSecondary()), str(sam.IsReverse()), num_unique_references);
			#print '---';
			
		num_references += 1;	# Count only unique sequences.
		
	fp_reference.close();

	for key in ret.keys():
		ret[key].sort(reverse=True, key=lambda sam_line: sam_line.chosen_quality);
		# ret[key].sort(reverse=True, key=lambda sam_line: ((sam_line.IsSecondary() == False), sam_line.chosen_quality));

	return [ret, num_references, num_unique_references];

# Hashes the entries from a SAM file by their QNAME for faster lookup during comparison.
def HashSAMWithFilter(sam_path, qname_hash_to_filter={}):
	try:
		fp_reference = open(sam_path, 'r');
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_path);
		return [{}, 0, 0];
	
	ret = {};
	
	num_unique_references = 0;
	num_references = 0;
	
	for line in fp_reference:
		line = line.strip();
		
		if len(line) == 0:
			continue;
		
		if line[0] == '@':
			continue;
		
		sam_line = SAMLine(line);
		
		modified_qname = sam_line.qname;
		#modified_qname = '/'.join(sam_line.qname.split('/')[:-1]);
		
		if (qname_hash_to_filter != {}):
			try:
				if (qname_hash_to_filter[modified_qname] == 1):
					sam_line.is_filtered_out = True;
			except:
					sam_line.is_filtered_out = False;
		
		try:
			current_hash = ret[modified_qname];
			
			should_be_counted = True;
			
			# We have to check if this sequence is actually a paired read of another sequence, or there is something funny in the data.
			for existing_sam in current_hash:
				if (sam_line.IsSecondary() == True) or (sam_line.IsSecondary() == existing_sam.IsSecondary() and sam_line.IsReverse() == existing_sam.IsReverse()):
					# This case should not be counted. It means that, either the alignment is marked as secondary which means there should be a primary alignment as well, or that there is more than one primary alignment, and the orientation is the same, which means that an aligner is trying to artificially boost-up their statistics.
					should_be_counted = False;
					break;
			
			if should_be_counted == True and sam_line.is_filtered_out == False:
				num_unique_references += 1;	# Count only unique sequences.

			# At least one sequence with the same name has already been added.
			current_hash.append(sam_line);
		except:
			# This is a new sequence, unhashed before. Create a new list and count the sequence.
			ret[modified_qname] = [sam_line];
			if sam_line.IsSecondary() == False and sam_line.is_filtered_out == False:
				num_unique_references += 1;	# Count only unique sequences, but not secondary ones.
		num_references += 1;	# Count only unique sequences.
		
	fp_reference.close();
	
	for key in ret.keys():
		ret[key].sort(reverse=True, key=lambda sam_line: sam_line.chosen_quality);
	
	return [ret, num_references, num_unique_references];

# This is deprecated, should be removed.
def CheckSamModified(sam_filename, out_path_prefix):
	sam_timestamp = str(os.path.getmtime(sam_filename));
	
	path_roc = out_path_prefix + '.roc';
	path_sum = out_path_prefix + '.sum';
	
	lines = [];
	
	try:
		fp_sum = open(path_sum, 'r');
		lines = fp_sum.readlines();
		fp_sum.close();
		
		fp_sum = open(path_roc, 'r');
		fp_sum.close();
	except IOError:
		return [True, sam_timestamp];

#	fp_roc.write('SAM timestamp: %s\n' % sam_timestamp);
	for line in lines:
		split_line = line.split(':');
		if split_line[0].strip() == 'SAM timestamp':
			if len(split_line) < 2:
				continue;
			last_sam_timestamp = split_line[1].strip();
			
			if last_sam_timestamp == sam_timestamp:
				return [False, sam_timestamp];
			
			break;

	return [True, sam_timestamp];

# Illumina and Roche reads that were generated using ART have different paired end notations in their headers. That
# makes sense if that simulates the original headers as would be produced by real technology.
# The problem is that, for some reason, ART trims the paired end notation (i.e. /1 and /2 for Illumina) from the qname
# field in the ground-truth SAM file. For this reason, function TrimHeader is used to remove the paired end notation from
# the Illumina and Roche headers.
def TrimHeader(header):
	ret = header;
	
	illumina_re = r'^(.*_[\d]+)-[12]';
	roche_re = r'^(.*-[\d]+)/[12]';
	match_header = None;
	if match_header == None:
		match_header = re.match(illumina_re, header);
	if match_header == None:
		match_header = re.match(roche_re, header);
		
	if match_header:
		ret = str(match_header.group(1));	
	
	return ret;

def GetExecutionTime(sam_file):
	time_file = sam_file[0:-3] + 'time';
	
	try:
		fp_time = open(time_file, 'r');
		execution_time = fp_time.readline();
		fp_time.close();
		
		return execution_time;
	except IOError:
		return 'Execution time measurements not found!'
	
def GetExecutionStats(sam_file):
	time_file = sam_file[0:-3] + 'memtime';
	
	try:
		fp_time = open(time_file, 'r');
		lines = fp_time.readlines();
		fp_time.close();
		
	except IOError:
		return 'Execution time measurements not found!'
	
	#for line in lines:
		#line = line.strip();
		#line
	
	return ('\t' + '\n\t'.join([line.strip() for line in lines]));

def ParseMemTime(sam_file):
	memtime_path = os.path.splitext(sam_file)[0] + '.memtime';
	fp = open(memtime_path, 'r');
	lines = [line.strip() for line in fp.readlines() if (len(line.strip()) > 0)];
	fp.close();

	cmdline = '';
	realtime = 0;
	cputime = 0;
	usertime = 0;
	systemtime = 0;
	maxrss = 0;
	rsscache = 0;
	time_unit = '';
	mem_unit = '';

	for line in lines:
		if (line.startswith('Command line:')):
			cmdline = line.split(':')[1].strip();
		elif (line.startswith('Real time:')):
			split_line = line.split(':')[1].strip().split(' ');
			realtime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('CPU time:')):
			split_line = line.split(':')[1].strip().split(' ');
			cputime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('User time:')):
			split_line = line.split(':')[1].strip().split(' ');
			usertime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('System time:')):
			split_line = line.split(':')[1].strip().split(' ');
			systemtime = float(split_line[0].strip());
			time_unit = split_line[1].strip();
		elif (line.startswith('Maximum RSS:')):
			split_line = line.split(':')[1].strip().split(' ');
			maxrss = float(split_line[0].strip());
			mem_unit = split_line[1].strip();
		# elif (line.startswith('')):
		# 	split_line = line.split(':')[1].strip().split(' ');
		# 	rsscache = float(split_line[0].strip());
		# 	mem_unit = split_line[1].strip();

	return [cmdline, realtime, cputime, usertime, systemtime, maxrss, time_unit, mem_unit];

def WriteSamLines(sam_lines, output_path):
	try:
		fp = open(output_path, 'w');
		i = 0;
		while i < len(sam_lines):
			sam_line = sam_lines[i];
#			temp_line = sam_line.VerboseFormatLine();
			temp_line = sam_line.FormatAccuracy();
			fp.write(temp_line + '\n');
			i += 1;
		fp.close();
	except IOError:
		sys.stderr.write('ERROR: Could not open file "%s" for writing!\n' % output_path);
		return;

def GetBasicStats(sam_lines, allowed_distance=1):
	true_positive = 0;
	false_positive = 0;
	not_mapped = 0;
	
	for sam_line in sam_lines:
		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1 and sam_line.min_distance <= allowed_distance):
				true_positive += 1;
			else:
				# print 'sam_line.min_distance > allowed_distance: %d > %d' % (sam_line.min_distance, allowed_distance);
				false_positive += 1;
		else:
			not_mapped = 0;
	
	return [true_positive, false_positive, not_mapped];

def GetDistanceHistogramStats(sam_lines, distance_limits):
	sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);
	
	sorted_distance_limits = sorted(distance_limits);
	ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	current_distance_index = 0;
	
	total_correct = 0;
	
	i = 0;
	while i < len(sorted_lines_by_distance):
		sam_line = sorted_lines_by_distance[i];
		
		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1):
				if (sam_line.min_distance < sorted_distance_limits[current_distance_index]):
					#ret_distance_counts[current_distance_index] += 1;
					total_correct += 1;
				else:
					ret_distance_counts[current_distance_index] = total_correct;
					current_distance_index += 1;
					if (current_distance_index >= len(sorted_distance_limits)):
						break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					continue;	# Do not increase i, because we want to check the current SAM line also.
				
		i += 1;

	if (current_distance_index > 0):
		ret_distance_counts[current_distance_index - 1] = total_correct;
		
	# If all the sam lines have been counted and the maximum observed distance is less than the maximum value of distance_limits parameter,
	# fill the rest of the array with the same value (the last value that was counted), to make the graph flatline.
	i = current_distance_index;
	while (i < len(sorted_distance_limits)):
		ret_distance_counts[i] = total_correct; # ret_distance_counts[current_distance_index];
		i += 1;

	# Nothing important, just format the X axis values for the counts.
	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];
	#return [sorted_distance_limits_shifted, ret_distance_counts];

	return [sorted_distance_limits, ret_distance_counts];

def GetDistanceHistogramStatsScaleDuplicates(sam_lines, distance_limits, scale_by_num_occurances=True):
	sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);
	
	sorted_distance_limits = sorted(distance_limits);
	ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	current_distance_index = 0;
	
	total_correct = 0.0;
	
	i = 0;
	while i < len(sorted_lines_by_distance):
		sam_line = sorted_lines_by_distance[i];
		
		if (sam_line.IsMapped()):
			#if (sam_line.is_correct_ref_and_orient == 1):
			if (sam_line.is_correct_ref_and_orient == 1 and sam_line.min_distance <= sorted_distance_limits[current_distance_index]):
				#ret_distance_counts[current_distance_index] += 1;
				if (scale_by_num_occurances == True):
					total_correct += (1.0 / float(sam_line.num_occurances_in_sam_file));
				else:
					#print 'total_correct = ', total_correct;
					total_correct += 1.0;
				
			else:
				if (sam_line.min_distance > sorted_distance_limits[current_distance_index]):
					ret_distance_counts[current_distance_index] = total_correct;
					current_distance_index += 1;
					if (current_distance_index >= len(sorted_distance_limits)):
						break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					continue;	# Do not increase i, because we want to check the current SAM line also.
				
		i += 1;

	if (current_distance_index > 0):
		ret_distance_counts[current_distance_index - 1] = total_correct;
		
	# If all the sam lines have been counted and the maximum observed distance is less than the maximum value of distance_limits parameter,
	# fill the rest of the array with the same value (the last value that was counted), to make the graph flatline.
	i = current_distance_index;
	while (i < len(sorted_distance_limits)):
		ret_distance_counts[i] = total_correct; # ret_distance_counts[current_distance_index];
		i += 1;

	# Nothing important, just format the X axis values for the counts.
	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];
	#return [sorted_distance_limits_shifted, ret_distance_counts];

	return [sorted_distance_limits, ret_distance_counts];

#def GetDistanceHistogramStats(sam_lines, distance_limits):
	#sorted_lines_by_distance = sorted(sam_lines, key=lambda sam_line: sam_line.min_distance);
	
	#sorted_distance_limits = sorted(distance_limits);
	#ret_distance_counts = [0 for distance_limit in sorted_distance_limits];

	#current_distance_index = 0;
	
	##for sam_line in sorted_lines_by_distance:
	#i = 0;
	#while i < len(sorted_lines_by_distance):
		#sam_line = sorted_lines_by_distance[i];
		
		#if (sam_line.IsMapped()):
			#if (sam_line.is_correct_ref_and_orient == 1):
				#if (sam_line.min_distance < sorted_distance_limits[current_distance_index]):
					#ret_distance_counts[current_distance_index] += 1;
				#else:
					#current_distance_index += 1;
					#if (current_distance_index >= len(sorted_distance_limits)):
						#break;
					#ret_distance_counts[current_distance_index] = 0 + ret_distance_counts[current_distance_index - 1];
					#continue;	# Do not increase i, because we want to check the current SAM line also.
				
		#i += 1;

	#i = current_distance_index;
	#while (i < len(sorted_distance_limits)):
		#ret_distance_counts[i] = ret_distance_counts[current_distance_index];
		#i += 1;

	#min_sorted_distance_limit = sorted_distance_limits[0];
	#sorted_distance_limits_shifted = [(distance_limit - min_sorted_distance_limit) for distance_limit in sorted_distance_limits];
	
	#return [sorted_distance_limits_shifted, ret_distance_counts];

def GetMapqHistogramStats(sam_lines, mapq_limits, allowed_distance):
	sorted_lines_by_mapq = sorted(sam_lines, reverse=True, key=lambda sam_line: sam_line.mapq);
	sorted_mapq_limits = sorted(mapq_limits, reverse=True);
	ret_mapq_counts = [0 for mapq_limit in sorted_mapq_limits];
	
	current_mapq = 0;
	total_correct = 0;
	#for sam_line in sorted_lines_by_mapq:
	i = 0;
	while i < len(sorted_lines_by_mapq):
		sam_line = sorted_lines_by_mapq[i];
		
		if (sam_line.IsMapped() and sam_line.is_duplicate == 0):
			if (sam_line.is_correct_ref_and_orient == 1):
				if (sam_line.mapq >= sorted_mapq_limits[current_mapq]):
					if (sam_line.min_distance < allowed_distance):
						#ret_mapq_counts[current_mapq] += 1;
						total_correct += 1;
				else:
					#print '[%d, %d] %d, total_correct = %d' % (current_mapq, sorted_mapq_limits[current_mapq], sam_line.mapq, total_correct);
					ret_mapq_counts[current_mapq] = total_correct;
					current_mapq += 1;
					if (current_mapq >= len(sorted_mapq_limits)):
						break;
					
					continue;	# Do not increase i, because we want to check the current SAM line also.
					#total_correct += 1;
					#ret_mapq_counts[current_mapq] = 0 + ret_mapq_counts[current_mapq - 1];
		i += 1;
		
	ret_mapq_counts[current_mapq] = total_correct;
	
	i = current_mapq;
	while (i < len(sorted_mapq_limits)):
		ret_mapq_counts[i] = ret_mapq_counts[current_mapq];
		i += 1;
	
	return [allowed_distance, sorted_mapq_limits, ret_mapq_counts];

def FindMultipleQnameEntries(sam_files):
	#duplicate_list = [];
	duplicate_hash = {};
	
	for sam_file in sam_files:
		try:
			fp_sam = open(sam_file, 'r');
		except IOError:
			sys.stderr.write('ERROR: Could not open file "%s" for reading!\n' % sam_file);
			return [{}, 0];
		
		occurance_hash = {};
		for line in fp_sam:
			line = line.strip();
			if len(line) == 0:
				continue;
			if line[0] == '@':
				continue;
			sam_line = SAMLine(line);
			
			try:
				occurance_hash[sam_line.qname] += 1;
				duplicate_hash[sam_line.qname] = 1;
			except:
				occurance_hash[sam_line.qname] = 1;
	
		fp_sam.close();
	
	fp = open('test.multiple', 'w');
	fp.write('\n'.join(duplicate_hash.keys()));
	fp.close();
	
	return duplicate_hash;



def CountMappedReads(sam_file):
	fp = None;
	
	try:
		fp = open(sam_file, 'r');
	except IOError:
		sys.stderr.write('[%s] ERROR: Could not open file "%s" for reading!\n' % ('CountMappedReads', sam_file));
		exit(1);
	
	num_alignments = 0;
	num_mapped_alignments = 0;
	sam_reads = {};
	sam_mapped_reads = {};
	num_mapped_bases = 0;	# Count of the number of non-clipped bases in an alignment. Only highest-scoring alignments are considered for each read,
	total_num_bases = 0;
	highest_scoring_alignment_for_read = {};	# Hash to keep track of the highest scoring alignment for a read.

	for line in fp:
		if (len(line) == 0 or line[0] == '@'):
			continue;
		
		sam_line = SAMLine(line.rstrip());
		
		# Count the number of occurances of each read in the SAM file (i.e. multiple mappings).
		try:
			sam_reads[sam_line.qname] += 1;
			#print sam_line.qname;
		except:
			sam_reads[sam_line.qname] = 1;
			
		if (sam_line.IsMapped() == True):
			# Count the occurances of mapped reads.
			try:
				sam_mapped_reads[sam_line.qname] += 1;
			except:
				sam_mapped_reads[sam_line.qname] = 1;
			
			try:
				if (highest_scoring_alignment_for_read[sam_line.qname].chosen_quality < sam_line.chosen_quality):
					highest_scoring_alignment_for_read[sam_line.qname] = sam_line;
			except:
				highest_scoring_alignment_for_read[sam_line.qname] = sam_line;
			
			num_mapped_alignments += 1;
		
		num_alignments += 1;
	
	fp.close();
	
	num_unique_reads = 0;
	for value in sam_reads.values():
		if (value == 1):
			num_unique_reads += 1;
	num_unique_mapped_reads = 0;
	for value in sam_mapped_reads.values():
		if (value == 1):
			num_unique_mapped_reads += 1;
	
	for best_read in highest_scoring_alignment_for_read.values():
		total_num_bases += len(sam_line.seq);
		num_mapped_bases += best_read.CalcNumMappedBases();

	num_mapped_reads = len(highest_scoring_alignment_for_read.keys());

	return [num_alignments, num_mapped_alignments, num_unique_reads, num_mapped_reads, num_mapped_bases];



### Counts the number of bases mapped to the same position in both SAMLines.
### To clarify the algorithm: two arrays are generated, each being of size of the read (before clipping).
### These will hold coordinates for each base on the reference.
### Then, for each of the two, fill the 'M=X' array positions with their coordinates on the reference.
### For the insertions, only a "-1" value is used, because this base is not actualy present on the reference.
### Both query_sam and ref_sam need to have a "-1" at the same position to be counted as a true inserted base.
### Deletions are not counted in this case.
def CompareBasePositions(query_sam, ref_sam, switch_ins_and_dels=False):
	qsam_ref_coords = [None] * query_sam.CalcReadLengthFromCigar();
	rsam_ref_coords = [None] * ref_sam.CalcReadLengthFromCigar();

	num_mapped_bases = len(qsam_ref_coords) - query_sam.clip_count_front - query_sam.clip_count_back;
	num_ref_bases = len(rsam_ref_coords) - ref_sam.clip_count_front - ref_sam.clip_count_back;

	if (len(qsam_ref_coords) < len(rsam_ref_coords) or
		query_sam.IsMapped() == False or ref_sam.IsMapped() == False or query_sam.rname != ref_sam.rname or query_sam.IsReverse() != ref_sam.IsReverse()):
		# sys.stderr.write('Warning: Mappers output does not conform to SAM format specification! CIGAR field does not specify sequence of equal length as in the input FASTA file. Possibly hard clipping operations are missing.\n');	
		return [0, num_mapped_bases, num_ref_bases];

	# Find the positions of all the query bases.
	query_cigpos = query_sam.CalcCigarStartingPositions(False);
	for cigpos in query_cigpos:
		[cig_count, cig_op, pos_on_ref, pos_on_query] = cigpos;
		if (cig_op in 'M=X'):
			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = range(pos_on_ref, (pos_on_ref + cig_count));
		elif (cig_op == 'I'):
			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-1]*cig_count;
		elif (cig_op == 'S'):
			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-3]*cig_count;
		elif (cig_op == 'H'):
			qsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-4]*cig_count;

	# Find the positions of all the reference bases.
	ref_cigpos = ref_sam.CalcCigarStartingPositions(False, switch_ins_and_dels);
	for cigpos in ref_cigpos:
		[cig_count, cig_op, pos_on_ref, pos_on_query] = cigpos;
		if (cig_op in 'M=X'):
			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = range(pos_on_ref, (pos_on_ref + cig_count));
		elif ((switch_ins_and_dels == False and cig_op == 'I') or (switch_ins_and_dels == True and cig_op == 'D')):
			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-1]*cig_count;
		elif (cig_op == 'S'):
			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-3]*cig_count;
		elif (cig_op == 'H'):
			rsam_ref_coords[pos_on_query:(pos_on_query + cig_count)] = [-4]*cig_count;

	# Count the number of correctly mapped bases ('M' and 'I' CIGAR operations):
	num_correct_bases = 0;
	i = 0;
	while (i < len(qsam_ref_coords)):
		# Skip the clipped bases:
		if (qsam_ref_coords[i] != -3 and qsam_ref_coords[i]!= -4):
			# Count the equal bases
			if (qsam_ref_coords[i] == rsam_ref_coords[i]):
				num_correct_bases += 1;
			# else:
			# 	# Just debug output.
			# 	print 'qsam[i] = %d\trsam[i] = %d\tseq[i] = %s' % (qsam_ref_coords[i], rsam_ref_coords[i], ref_sam.seq[i]);
		i += 1;

	return [num_correct_bases, num_mapped_bases, num_ref_bases];

def CountCorrectlyMappedBases(hashed_sam_lines, hashed_reference_sam, out_summary_prefix='', sam_basename='', switch_ins_and_dels=False):
	# if (use_strict == False):
	# 	return [0.0, 0.0, 0, 0];

	fp_out = None;
	out_file = out_summary_prefix + '.csv';
	if (out_summary_prefix != ''):
		try:
			fp_out = open(out_file, 'w');
		except IOError:
			sys.stderr.write('[%s] ERROR: Could not open file "%s" for writing!\n' % (__name__, out_file));
			exit(1);
	sys.stderr.write('Starting to count the number of correctly mapped bases in the tested SAM file!\n');


	total_ref_bases = 0;
	for qname in hashed_reference_sam.keys():
		ref_sam = hashed_reference_sam[qname][0];
		# total_ref_bases += (ref_sam.CalcReadLengthFromCigar() - ref_sam.clip_count_front - ref_sam.clip_count_back);
		num_ref_bases = len(ref_sam.seq);
		# Remove the counts of soft-clipped bases from the read length. We ignore hard clipped bases, because they are
		# not present in the SEQ field anyway.
		num_ref_bases -= (ref_sam.clip_count_front if ref_sam.clip_op_front == 'S' else 0);
		num_ref_bases -= (ref_sam.clip_count_back if ref_sam.clip_op_back == 'S' else 0);
		total_ref_bases += num_ref_bases;

	sum_correct_bases = 0;
	sum_mapped_bases = 0;
	sum_ref_bases = 0;

	i = 0;
	for qname in hashed_sam_lines.keys():
		i += 1;
		if ((i % 1000) == 0):
			sys.stderr.write('\rLine %d' % (i));
			sys.stderr.flush();

		sam_line = hashed_sam_lines[qname][0];

		# TODO: THIS NEEDS TO BE REMOVED OR IMPLEMENTED SOMEHOW DIFFERENTLY!!
		# The point of this was that, BLASR doesn't conform to the SAM standard, and makes it difficult to
		# uniformly evaluate the results!
		if 'blasr' in sam_basename.lower():
			qname = '/'.join(qname.split('/')[:-1]);
			if sam_line.clip_count_front != 0 or sam_line.clip_count_back != 0:
				sys.stderr.write('BLASR CIGAR contains clipping! Please revise clipped_pos! Read: "%s".\n' % sam_line.qname);

		if (sam_line.IsMapped() == False):
			continue;
		if ((qname in hashed_reference_sam) == False):
			sys.stderr.write('\tERROR: Reference SAM does not contain qname "%s"!\n' % (qname));
			continue;

		sam_reference = hashed_reference_sam[qname][0];

		[num_correct_bases, num_mapped_bases, num_ref_bases] = CompareBasePositions(sam_line, sam_reference, switch_ins_and_dels);
		sum_correct_bases += num_correct_bases;
		sum_mapped_bases += num_mapped_bases;
		sum_ref_bases += num_ref_bases;

		# if (float(num_correct_bases) / float(num_mapped_bases) < 0.75):
		# 	print 'Original line:';
		# 	print sam_line.original_line;
		# 	print '';
		# 	print 'Reference SAM line:';
		# 	print sam_reference.original_line;
		# 	print '';
		# 	print '';

	precision = (100.0 * float(sum_correct_bases) / float(sum_mapped_bases)) if (sum_mapped_bases > 0) else 0.0;
	recall = (100.0 * float(sum_correct_bases) / float(total_ref_bases)) if (total_ref_bases) else 0.0;

	if (out_summary_prefix != ''):
		fp_out.write('percent_correct_m\tnum_correct_m\tnum_m_ops_in_reference\n');
		fp_out.write('%.2f\t%.2f\t%.2f\t%.2f\n' % (precision, recall, sum_correct_bases, sum_mapped_bases));
		fp_out.close();
	sys.stderr.write('\n');

	return [precision, recall, sum_correct_bases, sum_mapped_bases, total_ref_bases];



if __name__ == "__main__":
	pass;
