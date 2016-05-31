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
def log(message, fp_log, silent=False):
    if (silent == True): return;

    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
    if (message != ''):
    	prefix = '[%s] ' % (timestamp);
    else:
    	prefix = '';

    sys.stderr.write('%s%s\n' % (prefix, message));
    sys.stderr.flush();
    if (fp_log != None and fp_log != sys.stderr):
        fp_log.write('%s%s\n' % (prefix, message));
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

def execute_command_with_ret(dry_run, command, silent=False):
	if (silent == False):
		sys.stderr.write('Executing command: "%s"\n' % command);
	if (dry_run == False):
		p = subprocess.Popen(command, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE);
	[output, err] = p.communicate()
	rc = p.returncode
	if (silent == False):
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

def eval_contigs(ref_path, contig_path, temp_folder, generate_kmer_spectrum=False, silent=False):
	if (not os.path.exists(temp_folder)):
		os.makedirs(temp_folder);

	[headers_contigs, seqs_contigs, quals_contigs] = fastqparser.read_fastq(contig_path);
	[headers_ref, seqs_ref, quals_ref] = fastqparser.read_fastq(ref_path);

	ref_hash = hash_headers(headers_ref);
	contig_hash = hash_headers(headers_contigs);

	avg_accuracy_overall = 0.0;
	avg_id_overall = 0.0;

	single_contig_path = '%s/singlecontig.fasta' % (temp_folder);
	for i in xrange(0, len(seqs_contigs)):
		contig_name = headers_contigs[i].split()[0];
		contig_seq = seqs_contigs[i];

		fp_contig = open(single_contig_path, 'w');
		fp_contig.write('>%s\n%s\n' % (contig_name, seqs_contigs[i]));
		fp_contig.close();

		### Run nucmer to align the contig to the reference, also, filter the delta file and generate alignment coordinates.
		nucmer_out_prefix = '%s/nucmer' % (temp_folder);
		log('', sys.stderr, silent=silent);
		log('Running MUMmer on contig: "%s"' % (contig_name), sys.stderr, silent=silent);
		command = '%s --maxmatch --extend -p %s %s %s; delta-filter -r -q %s.delta > %s.filt.delta; show-coords -r -c %s.filt.delta > %s.filt.coords' % \
					(NUCMER_PATH, nucmer_out_prefix, ref_path, single_contig_path, nucmer_out_prefix, nucmer_out_prefix, nucmer_out_prefix, nucmer_out_prefix);
		[rc, rstdout, rstderr] = execute_command_with_ret(DRY_RUN, command, silent=True);

		### Load the coordinates.
		log('Parsing the coords file.', sys.stderr, silent=silent);
		# fp = open('/home/isovic/work/eclipse-workspace/git/consise/temp2-mummer/test-data/out/nucmer.coords2', 'r');
		coords_path = '%s.filt.coords' % (nucmer_out_prefix);
		fp = open(coords_path, 'r');
		lines = fp.readlines();
		fp.close();

		frags = parse_coords_lines(lines, contig_name, seqs_ref, ref_hash, seqs_contigs, contig_hash);

		avg_accuracy_contig = 0.0;
		avg_id_contig = 0.0;

		log('Running Edlib to determine the edit distance for each fragment...', sys.stderr, silent=silent);
		for frag in frags:
			# print frag;

			[rstart, rend, qstart, qend, fwd, rname, qname, identity] = frag;
			ref_seq = seqs_ref[ref_hash[rname]];
			[nw_ref, nw_contig] = extract_seqs_for_edlib(ref_seq, contig_seq, rstart, rend, qstart, qend);

			temp_suffix = '.%d' % (i);
			nw_ref_path = '%s/nw-ref%s.fasta' % (temp_folder, temp_suffix);
			nw_contig_path = '%s/nw-contig%s.fasta' % (temp_folder, temp_suffix);

			fp_nw_ref = open(nw_ref_path, 'w');
			fp_nw_contig = open(nw_contig_path, 'w');
			fp_nw_ref.write('>%s\n%s\n' % (rname, nw_ref));
			fp_nw_contig.write('>%s\n%s\n' % (qname, nw_contig));
			fp_nw_ref.close();
			fp_nw_contig.close();

			command = '%s %s %s -m NW' % (EDLIB_PATH, nw_contig_path, nw_ref_path);
			[rc, rstdout, rstderr] = execute_command_with_ret(DRY_RUN, command, silent=True);
			scores = parse_edlib_scores(rstdout);
			if (len(scores) == 0): log('ERROR: len(scores) == 0!\nreturn code: %d\nrstdout:\n%s' % (rc, rstdout), sys.stderr);
			# sys.stderr.write('Final edit distance: %d\n' % (scores[0]));

			normalized_edit_dist = float(scores[0]) / float(abs(qend - qstart + 1));
			accuracy = (1.0 - normalized_edit_dist);

			frag.append(scores[0]);
			frag.append(100.0 * normalized_edit_dist);
			frag.append(100.0 * accuracy);
			# print frag;

			avg_accuracy_contig += accuracy;
			avg_id_contig += frag[7];

		avg_accuracy_contig /= float(len(frags));
		avg_id_contig /= float(len(frags));
		log('Average ID for contig "%s": %f%%' % (contig_name, avg_id_contig), sys.stderr, silent=silent);
		log('Average accuracy for contig "%s": %f%%' % (contig_name, 100.0*avg_accuracy_contig), sys.stderr, silent=silent);
		log('', sys.stderr, silent=silent);

		avg_accuracy_overall += avg_accuracy_contig;
		avg_id_overall += avg_id_contig;

	avg_accuracy_overall /= float(len(seqs_contigs));
	avg_id_overall /= float(len(seqs_contigs));

	log('Draft assembly: "%s"\n' % (contig_path), sys.stderr, silent=silent);
	log('Overall average ID for the draft assembly: %f%%' % (avg_id_overall), sys.stderr, silent=silent);
	log('Overall average accuracy for the draft assembly: %f%%' % (100.0*avg_accuracy_overall), sys.stderr, silent=silent);
	log('', sys.stderr);

	# sys.stdout.write('=============================================\n');
	sys.stdout.write('================= Summary ===================\n');
	sys.stdout.write('Draft assembly: "%s"\n' % (contig_path));
	sys.stdout.write('Overall average ID for the draft assembly: %f%%\n' % (avg_id_overall));
	sys.stdout.write('Overall average accuracy for the draft assembly: %f%%\n' % (100.0*avg_accuracy_overall));
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
			continue;
		if (state_coord_lines == False):
			continue;
		if (line.endswith(contig_name) == False):
			continue;

		sline = line.split();
		[curr_rstart, curr_rend, curr_qstart, curr_qend, curr_rname, curr_qname, curr_identity] = [int(sline[0]), int(sline[1]), int(sline[3]), int(sline[4]), sline[-2], sline[-1], float(sline[9])];
		fwd = True if (curr_qstart <= curr_qend) else False;
		# print 'fwd = ', fwd;
		# if (len(all_frags) > 0 and curr_rname != all_frags[-1][-2]):
		# 	sys.stderr.write('ERROR: Fragments of one contig are not aligned to the same reference! Possible structural variant? Exiting.\n');
		# 	exit(1);
		all_frags.append([curr_rstart, curr_rend, curr_qstart, curr_qend, fwd, curr_rname, curr_qname, curr_identity]);
		num_fwd += (1 if (fwd) else 0);
		num_coord_lines += 1;


	return all_frags;

def hash_headers(headers):
	ret = {};
	for i in xrange(0, len(headers)):
		ret[headers[i]] = i;
		ret[headers[i].split()[0]] = i;
	return ret;

def extract_seqs_for_edlib(ref_seq, contig_seq, rstart, rend, qstart, qend):
	if (rend > rstart):
		nw_ref = ref_seq[(rstart-1):(rend+1-1)];	# +1 because the end base is inclusive, and -1 because it's 1-based.
	else:
		nw_ref = fastqparser.revcomp_seq(ref_seq[(rend-1):(rstart+1-1)]);

	if (qend > qstart):
		nw_contig = contig_seq[(qstart-1):(qend+1-1)];	# +1 because the end base is inclusive, and -1 because it's 1-based.
	else:
		nw_contig = fastqparser.revcomp_seq(contig_seq[(qend-1):(qstart+1-1)]);

	return [nw_ref, nw_contig];

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
	silent_mode = False;

	for i in xrange(3, len(sys.argv)):
		arg = sys.argv[i];
		if (arg == '--temp-path'):
			temp_path = sys.argv[i+1];
			i += 1;
		elif (arg == '--spectrum'):
			generate_kmer_spectrum = True;
		elif (arg == '--silent'):
			silent_mode = True;

	if (temp_path == None):
		if (len(os.path.dirname(contigs_path)) != 0):
			temp_path = os.path.dirname(contigs_path) + '/temp-ed';
		else:
			temp_path = './temp-ed';

	eval_contigs(ref_path, contigs_path, temp_path, generate_kmer_spectrum=generate_kmer_spectrum, silent=silent_mode);

if __name__ == "__main__":
	main();
