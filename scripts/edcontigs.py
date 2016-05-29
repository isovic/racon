#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));
import sys;
sys.path.append(SCRIPT_PATH + '/../codebase/samscripts/src/');

import subprocess;
import operator;

import utility_sam;
import fastqparser;

import copy;

from lis import *;

EDLIB_PATH = '%s/../tools/edlib/src/aligner' % (SCRIPT_PATH);
KMERCOMP_PATH = '%s/../tools/edlib/src/kmercomp' % (SCRIPT_PATH);
NUCMER_PATH = 'nucmer';

DRY_RUN = False;

import traceback;
from time import gmtime, strftime
### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());

    sys.stderr.write('[%s] %s\n' % (timestamp, message))
    sys.stderr.flush();
    if (fp_log != None):
        fp_log.write('[%s] %s\n' % (timestamp, message))
        fp_log.flush();
def execute_command(command, fp_log, dry_run=True):
    if (dry_run == True):
        log('Executing (dryrun): "%s".' % (command), fp_log);
        log('\n', fp_log);
        return 0;
    if (dry_run == False):
        log('Executing: "%s".' % (command), fp_log);
        rc = subprocess.call(command, shell=True);
        if (rc != 0):
            log('ERROR: subprocess call returned error code: %d.' % (rc), fp_log);
            log('Traceback:', fp_log);
            traceback.print_stack(fp_log);
            exit(1);
        return rc;

def execute_command_with_ret(dry_run, command):
	sys.stderr.write('Executing command: "%s"\n' % command);
	if (dry_run == False):
		p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	[output, err] = p.communicate()
	rc = p.returncode
	sys.stderr.write('\n');
	return [rc, output, err];

def parse_edlib_scores(rstdout):
	lines = rstdout.split('\n');
	[STATE_INIT, STATE_SCORES, STATE_END] = [0, 1, 2];
	current_state = STATE_INIT;
	scores = {};
	for line in lines:
		line = line.strip();
		if (len(line) == 0): continue;
		if (current_state != STATE_SCORES and line.startswith('<query number>:')):
			current_state = STATE_SCORES;
		elif (current_state == STATE_SCORES and line.startswith('Cpu time of searching:')):
			current_state = STATE_END;
		elif (current_state == STATE_SCORES):
			split_line = line.split(' ');
			qid = int(split_line[0][1:-1]);			# Remove the '#' and ':' chars.
			score = int(split_line[1]);			# The actual score for the query.
			scores[qid] = score;
	return scores;

def get_circular_score(ref_path, contig_path, temp_folder):
	if (not os.path.exists(temp_folder)):
		os.makedirs(temp_folder);

	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	circularized_fwd_path = '%s/circ-fwd.fa' % (temp_folder);
	circularized_rev_path = '%s/circ-rev.fa' % (temp_folder);

	fp_fwd = open(circularized_fwd_path, 'w');
	fp_rev = open(circularized_rev_path, 'w');

	for i in xrange(0, len(seqs_ref)):
		rev_seq = fastqparser.revcomp_seq(seqs_ref[i]);
		rev_qual = quals_ref[i][::-1];
		# if (len(quals_ref) > 0):
		# 	fp_fwd.write('@%s\n%s%s\n+\n%s%s\n' % (headers_ref[i], seqs_ref[i], seqs_ref[i], quals_ref[i], quals_ref[i]));
		# 	fp_rev.write('@%s\n%s%s\n+\n%s%s\n' % (headers_ref[i], rev_seq, rev_seq, rev_qual, rev_qual));
		# else:
		fp_fwd.write('>%s\n%s%s\n' % (headers_ref[i], seqs_ref[i], seqs_ref[i]));
		fp_rev.write('>%s\n%s%s\n' % (headers_ref[i], rev_seq, rev_seq));

	fp_fwd.close();
	fp_rev.close();

	# sys.stdout.write('Aligning the fwd orientation...\n');
	# command = '%s %s %s -m HW' % (EDLIB_PATH, contig_path, circularized_fwd_path);
	# [rc_fwd, rstdout_fwd, rstderr_fwd] = execute_command_with_ret(DRY_RUN, command);
	# scores_fwd = parse_edlib_scores(rstdout_fwd);
	# for i in xrange(0, len(scores_fwd)):
	# 	sys.stdout.write('[%d] %d %s\n' % (i, scores_fwd[i], 'fwd'));
	# sys.stdout.write('\n');

	sys.stdout.write('Aligning the rev orientation...\n');
	command = '%s %s %s -m HW' % (EDLIB_PATH, contig_path, circularized_rev_path);
	[rc_rev, rstdout_rev, rstderr_rev] = execute_command_with_ret(DRY_RUN, command);
	scores_rev = parse_edlib_scores(rstdout_rev);
	for i in xrange(0, len(scores_rev)):
		sys.stdout.write('[%d] %d %s\n' % (i, scores_rev[i], 'rev'));
	sys.stdout.write('\n');

	# sys.stdout.write('Maximum scores of both orientations:\n');
	# for i in xrange(0, len(scores_fwd)):
	# 	score_max = scores_fwd[i] if (scores_fwd[i] >= scores_rev[i]) else scores_rev[i];
	# 	orient = 'fwd' if (scores_fwd[i] >= scores_rev[i]) else 'rev';
	# 	sys.stdout.write('[%d] %d %s\n' % (i, score_max, orient));

def eval_contigs(ref_path, contig_path, temp_folder, generate_kmer_spectrum=False):
	if (not os.path.exists(temp_folder)):
		os.makedirs(temp_folder);

	[headers_contigs, seqs_contigs, quals_contigs] = fastqparser.read_fastq(contig_path);
	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	ref_hash = hash_headers(headers_ref);
	contig_hash = hash_headers(headers_contigs);

	single_contig_path = '%s/singlecontig.fasta' % (temp_folder);
	for i in xrange(0, len(seqs_contigs)):
		contig_name = headers_contigs[i].split()[0];
		contig_seq = seqs_contigs[i];

		fp_contig = open(single_contig_path, 'w');
		fp_contig.write('>%s\n%s\n' % (contig_name, seqs_contigs[i]));
		fp_contig.close();

		nucmer_out_prefix = '%s/nucmer' % (temp_folder);
		sys.stderr.write('\n');
		sys.stderr.write('Running MUMmer on contig: "%s"\n' % (contig_name));
		command = '%s --maxmatch --extend -p %s %s %s; delta-filter -r -q %s.delta > %s.filt.delta; show-coords -r -c %s.filt.delta > %s.filt.coords' % \
					(NUCMER_PATH, nucmer_out_prefix, ref_path, single_contig_path, nucmer_out_prefix, nucmer_out_prefix, nucmer_out_prefix, nucmer_out_prefix);
		# execute_command(command, None, False);
		[rc, rstdout, rstderr] = execute_command_with_ret(DRY_RUN, command);

		sys.stderr.write('\n');
		sys.stderr.write('Parsing the coords file.\n');
		# fp = open('/home/isovic/work/eclipse-workspace/git/consise/temp2-mummer/test-data/out/nucmer.coords2', 'r');
		coords_path = '%s.filt.coords' % (nucmer_out_prefix);
		fp = open(coords_path, 'r');
		lines = fp.readlines();
		fp.close();
		coords = parse_coords_lines(lines, contig_name, seqs_ref, ref_hash, seqs_contigs, contig_hash);
		print 'coords: "%s"' % (coords);
		print 'lines:\n', lines;
		sys.stdout.flush();
		[rstart, rend, qstart, qend, is_fwd, rname, qname] = coords;
		extract_seqs_for_edlib(temp_folder, '.%d' % (i), ref_path, contig_path, rstart, rend, qstart, qend, is_fwd, rname, qname, generate_kmer_spectrum=generate_kmer_spectrum);
		sys.stderr.write('\n');



def parse_coords_lines(lines, contig_name, seqs_ref, ref_hash, seqs_contigs, contig_hash):
	state_coord_lines = False;
	num_coord_lines = 0;
	rstart = 0;
	rend = 0;
	qstart = 0;
	qend = 0;
	all_frags = [];
	num_fwd = 0;

	for line in lines:
		line = line.strip();
		if (len(line) == 0): continue;
		if (line == '='*len(line)):
			state_coord_lines = True;
			print 'state_coord_lines = ', state_coord_lines;
			continue;
		if (state_coord_lines == False):
			continue;
		if (line.endswith(contig_name) == False):
			continue;

		print line;

		sline = line.split();
		print sline;
		[curr_rstart, curr_rend, curr_qstart, curr_qend, curr_rname, curr_qname] = [int(sline[0]), int(sline[1]), int(sline[3]), int(sline[4]), sline[-2], sline[-1]];
		fwd = True if (curr_qstart <= curr_qend) else False;
		print 'fwd = ', fwd;
		if (len(all_frags) > 0 and curr_rname != all_frags[-1][-2]):
			sys.stderr.write('ERROR: Fragments of one contig are not aligned to the same reference! Possible structural variant? Exiting.\n');
			exit(1);
		all_frags.append([curr_rstart, curr_rend, curr_qstart, curr_qend, fwd, curr_rname, curr_qname]);
		num_fwd += (1 if (fwd) else 0);
		num_coord_lines += 1;

	# for frag in all_frags:
	# 	print frag;
	# print num_fwd;
	correct_orient = True if (num_fwd > (len(all_frags) / 2)) else False;
	print 'correct_orient = ', correct_orient;
	all_frags = [val for val in all_frags if val[4] == correct_orient];
	print 'Printing frags:';
	for frag in all_frags:
		print frag;
	sys.stdout.flush();

	if (len(all_frags) == 0):
		return [rstart, rend, qstart, qend, True, '', ''];

	### This is experimental:
	# num_frags = len(all_frags);
	# ext_frags = copy.copy(all_frags);
	# for i in xrange(0, num_frags):
	# 	frag = ext_frags[i];
	# 	contig_len = len(seqs_contigs[contig_hash[frag[6]]]);
	# 	ref_len = len(seqs_ref[ref_hash[frag[5]]]);
	# 	# all_frags.append([frag[0] + ref_len, frag[1] + ref_len, frag[2] + contig_len, frag[3] + contig_len, frag[4], frag[5], frag[6]]);
	# 	frag[2] = contig_len - frag[2];
	# 	frag[3] = contig_len - frag[3];
	# 	ext_frags.append([frag[0] + ref_len, frag[1] + ref_len, frag[2] + 0, frag[3] + 0, frag[4], frag[5], frag[6]]);
	# 	# all_frags.append([frag[0] + ref_len, frag[1] + ref_len, contig_len - frag[2] - 1, contig_len - frag[3] - 1, frag[4], frag[5], frag[6]]);
	# print 'contig_len = %d' % (contig_len);
	# print 'ref_len = %d' % (ref_len);
	# print '';
	# for i in xrange(0, len(ext_frags)):
	# 	if (i == num_frags): print '-';
	# 	print ext_frags[i];
	# print '';

	# lis_in = [];
	# for i in xrange(0, len(ext_frags)):
	# 	lis_in.append(ext_frags[i][2]);
	# 	print lis_in[-1];
	# print '';
	# lis_out = longest_increasing_subsequence(lis_in);
	# for i in xrange(0, len(lis_out)):
	# 	print lis_out[i];



	# joined_frags = [];
	# if (correct_orient == True):
	# 	for frag in all_frags:
	# 		if (len(joined_frags) == 0):
	# 			joined_frags.append([frag[0] + 0, frag[1] + 0, frag[2] + 0, frag[3] + 0, frag[4], frag[5], frag[6]]);
	# 		elif (frag[0] >= joined_frags[-1][0] and frag[2] >= joined_frags[-1][2]):
	# 			joined_frags[-1][1] = frag[1] + 0;
	# 			joined_frags[-1][3] = frag[3] + 0;
	# 		else:
	# 			joined_frags.append([frag[0] + 0, frag[1] + 0, frag[2] + 0, frag[3] + 0, frag[4], frag[5], frag[6]]);


	# print '';
	# for frag in joined_frags:
	# 	print frag;
	# print;



	rstart = all_frags[0][0];
	rend = all_frags[-1][1];
	qstart = all_frags[0][2];
	qend = all_frags[-1][3];
	rname = all_frags[0][-2];
	qname = all_frags[0][-1];

	return [rstart, rend, qstart, qend, correct_orient, rname, qname];

def hash_headers(headers):
	ret = {};
	for i in xrange(0, len(headers)):
		ret[headers[i]] = i;
		ret[headers[i].split()[0]] = i;
	return ret;

def extract_seqs_for_edlib(temp_folder, temp_suffix, ref_path, contig_path, rstart, rend, qstart, qend, is_fwd, rname, qname, generate_kmer_spectrum=False):
	if (not os.path.exists(temp_folder)):
		os.makedirs(temp_folder);

	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);
	[headers_contig, seqs_contig, quals_contig] = fastqparser.read_fastq(contig_path);

	ref_hash = hash_headers(headers_ref);
	contig_hash = hash_headers(headers_contig);

	print ref_hash;

	if ((rname in ref_hash) == False):
		sys.stderr.write('ERROR: Reference name "%s" not found in file "%s"! Exiting.\n' % (rname, ref_path));
		exit(1);
	if ((qname in contig_hash) == False):
		sys.stderr.write('ERROR: Contig name "%s" not found in file "%s"! Exiting.\n' % (qname, contig_path));
		exit(1);

	if (rend < rstart):
		sys.stderr.write('ERROR: Reference end should come before reference start (it is expected that the ref is forward oriented), but ref_start = %d, ref_end = %d. Exiting.\n' % (ref_start, ref_end));
		exit(1);

	rid = ref_hash[rname];
	ref_header = headers_ref[rid];
	ref_seq = seqs_ref[rid][(rstart-1):(rend)];	# Coordinates are 1-based.

	qid = contig_hash[qname];
	contig_header = headers_contig[qid];
	contig_seq = '';
	if (is_fwd):
		if (qend >= qstart):
			contig_seq = seqs_contig[qid][(qstart-1):(qend)];
		else:
			contig_seq = seqs_contig[qid][(qstart-1):] + seqs_contig[qid][0:(qend)];
	else:
		if (qend > qstart):
			contig_seq = seqs_contig[qid][(qend-1):] + seqs_contig[qid][0:(qstart)]
		else:
			contig_seq = seqs_contig[qid][(qend-1):(qstart)];

		contig_seq = fastqparser.revcomp_seq(contig_seq);

	nw_ref_path = '%s/nw-ref%s.fasta' % (temp_folder, temp_suffix);
	nw_contig_path = '%s/nw-contig%s.fasta' % (temp_folder, temp_suffix);
	nw_kmer_comp_path ='%s/nw-kmers%s.spect' % (temp_folder, temp_suffix);

	fp_nw_ref = open(nw_ref_path, 'w');
	fp_nw_contig = open(nw_contig_path, 'w');
	fp_nw_ref.write('>%s\n%s\n' % (ref_header, ref_seq));
	fp_nw_contig.write('>%s\n%s\n' % (contig_header, contig_seq));
	fp_nw_ref.close();
	fp_nw_contig.close();

	sys.stderr.write('Running Edlib to determine the edit distance...\n');
	command = '%s %s %s -m NW' % (EDLIB_PATH, nw_contig_path, nw_ref_path);

	[rc, rstdout, rstderr] = execute_command_with_ret(DRY_RUN, command);
	# execute_command(command, None, False);
	scores = parse_edlib_scores(rstdout);
	unaligned_len = len(seqs_ref[rid]) - len(ref_seq);
	if (len(scores) == 0): sys.stderr.write('ERROR: len(scores) == 0!\nreturn code: %d\nrstdout:\n%s\n' % (rc, rstdout));
	sys.stderr.write('Final edit distance: %d, aligned edit distance: %d, unaligned ref len: %d, aligned ref len: %d, aligned contig len: %d\n' % ((scores[0] + unaligned_len), scores[0], unaligned_len, len(ref_seq), len(contig_seq)));
#	for i in xrange(0, len(scores)):
#		sys.stdout.write('[%d] edit dist: %d\tunaligned len: %d\n' % (i, scores[i], unaligned_len));

	sys.stdout.write('\n');


	if (generate_kmer_spectrum == True):
		sys.stderr.write('Generating the kmer spectrum.\n');
		command = '%s -o %s %s %s' % (KMERCOMP_PATH, nw_kmer_comp_path, nw_contig_path, nw_ref_path);
		[rc, rstdout, rstderr] = execute_command_with_ret(DRY_RUN, command);
		sys.stderr.write('Stdout:\n%s\nStderr:\n%s\n' % (rstdout, rstderr));
		sys.stderr.write('Done generating the kmer spectrum!\n');

#	[rc_rev, rstdout_rev, rstderr_rev] = execute_command_with_ret(DRY_RUN, command);
#	scores_rev = parse_edlib_scores(rstdout_rev);
#	for i in xrange(0, len(scores_rev)):
#		sys.stdout.write('[%d] %d %s\n' % (i, scores_rev[i], 'rev'));
#	sys.stdout.write('\n');



def main():
	# ./edcontigs.py temp2-mummer/ecoli_K12_MG1655_U00096.3.fasta temp-ecoli_pacbio/consensus-ecoli_pacbio-poa-v3-pileup.fasta

	if (len(sys.argv) < 3):
		sys.stderr.write('Tool for calculating the edit distance of two sequences using Edlib.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('\t%s <ref.fa> <contigs.fa> [options]\n' % (sys.argv[0]));
		exit(1);

	ref_path = sys.argv[1];
	contigs_path = sys.argv[2];

	generate_kmer_spectrum = False;
	temp_path = None;

	for i in xrange(3, len(sys.argv)):
		arg = sys.argv[i];
		if (arg == '--temp-path'):
			temp_path = sys.argv[i+1];
			i += 1;
		elif (arg == '--spectrum'):
			generate_kmer_spectrum = True;

	if (temp_path == None):
		if (len(os.path.dirname(contigs_path)) != 0):
			temp_path = os.path.dirname(contigs_path) + '/temp-ed';
		else:
			temp_path = './temp-ed';

	eval_contigs(ref_path, contigs_path, temp_path, generate_kmer_spectrum=generate_kmer_spectrum);

if __name__ == "__main__":
	main();
