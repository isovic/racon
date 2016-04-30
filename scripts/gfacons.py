#! /usr/bin/python

import os
SCRIPT_PATH = os.path.dirname(os.path.realpath(__file__));

import sys;
sys.path.append(SCRIPT_PATH + '/../codebase/samscripts/src/');

import subprocess;

import fastqparser;

# VERBOSE_DEBUG = True;
VERBOSE_DEBUG = False;

import traceback;
from time import gmtime, strftime
### Logs messages to STDERR and an output log file if provided (opened elsewhere).
def log(message, fp_log):
    timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());

    sys.stderr.write('[%s] %s\n' % (timestamp, message))
    sys.stderr.flush();
    if (fp_log != None and fp_log != sys.stderr):
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

class GoldenPath:
	def __init__(self, line=None):
		if (line == None):
			self.utg_name = '';
			self.utg_start = 0;
			self.read_name = '';
			self.rstart = 0;
			self.rend = 0;
			self.orient = '';
			self.inc_len = 0;
		else:
			split_line = line.split('\t');
			keyword = split_line[0];
			self.utg_name = split_line[1];
			self.utg_start = int(split_line[2]);
			buff = split_line[3].split(':');
			self.orient = split_line[4];
			self.inc_len = int(split_line[5]);
			self.read_name = buff[0];
			buff_pos = buff[1].split('-');
			self.rstart = int(buff_pos[0]);
			self.rend = int(buff_pos[1]);

class Contig:
	def __init__(self, name, seq):
		self.name = name;
		self.seq = seq;
		self.gp = [];	# Golden path from the output of Miniasm. This is a list of GoldenPath objects.

def parse_gfa(gfa_path):
	try:
		fp = open(gfa_path, 'r');
	except Exception, e:
		sys.stderr.write('ERROR: Could not open file "%s" for reading! Exiting.\n' % (gfa_path));
		sys.stderr.write(str(e));
		exit(1);

	contig_dict = {};
	contigs = [];

	for line in fp:
		line = line.strip();
		if (len(line) == 0): continue;
		sline = line.split('\t');
		if (sline[0] == 'S'):
			[keyword, seqname, seq, seqlen] = sline;
			contigs.append(Contig(seqname, seq));
			contig_dict[seqname] = contigs[-1];
		elif (sline[0] == 'a'):
			utg_name = sline[1];
			if (not utg_name in contig_dict):
				sys.stderr.write('ERROR: utg_name not found in previously loaded contig_dict! Malformed GFA file? utg_name = "%s". Exiting.\n' % (utg_name));
				exit(1);
			contig_dict[utg_name].gp.append(GoldenPath(line));

	return contigs;

def verbose_gfa_contigs(contigs):
    for contig in contigs:
    	print contig.name;
    	for i in xrange(0, len(contig.gp)):
    		print '\t' + contig.gp[i].read_name;

# Nekoliko stvari za isprobati:
# 1. Napraviti alignment ili overlap readova na referencu (raw contig), da nadjem koji se readovi preklapaju. Kada popravljam samo jedan read, onda zapisati samo one s kojima se preklapa u file.
# 2. Umjesto popravljanja svih readova u pathu popraviti samo one koji su najdalje preklopljeni da se smanji broj operacija.

# 3. Neovisno o ovome gore, trebao bih: popraviti region selection da smanjim IF-ove tako da stavim self hitove izvan counting dijela, i trebao bih IndexSpacedHashFast srediti da bude MFS to LFS a ne obrnuto.

def correct_layout_reads(contig, reads_path, num_threads, out_consensus_path, ref_path=None):
	[headers_reads, seqs_reads, quals_reads] = fastqparser.read_fastq(reads_path);
	seq_hash = {};
	for i in xrange(0, len(headers_reads)):
		seq_hash[headers_reads[i]] = [headers_reads[i], seqs_reads[i], quals_reads[i]];
		seq_hash[headers_reads[i].split()[0]] = [headers_reads[i], seqs_reads[i], quals_reads[i]];		# Handle the case when read name was split at a whitespace.
		seq_hash[headers_reads[i].split(':')[0]] = [headers_reads[i], seqs_reads[i], quals_reads[i]];		# Handle the case of nanopore reads in FASTQ files, where ':' is used as a delimiter instead of space.

	single_read_path = '%s/../temp/singleread.fastq' % (SCRIPT_PATH);
	single_read_sam = '%s/../temp/singleread.sam' % (SCRIPT_PATH);
	single_read_consensus = '%s/../temp/singleread.cons.fa' % (SCRIPT_PATH);
	all_corrected_reads_path = '%s/../temp/reads.cons.fa' % (SCRIPT_PATH);
	graphmap_bin = '%s/../tools/graphmap/bin/Linux-x64/graphmap' % (SCRIPT_PATH);
	all_corrected_reads_sam = '%s/../temp/reads.cons.sam' % (SCRIPT_PATH);
	aln_reads_to_raw_contig = '%s/../temp/reads_to_raw_contig.sam' % (SCRIPT_PATH);

	raw_contig_temp_path = '%s/../temp/contig.raw.fa' % (SCRIPT_PATH);
	final_consensus_contig_path = '%s/../temp/contig.consensus.fa' % (SCRIPT_PATH);

	# ### Alignment svih readova na contig, da se odredi koji se alignaju s kojima.
	# command = '%s align -r %s -d %s -o %s --evalue 0 --mapq 40 -v 1' % (graphmap_bin, raw_contig_temp_path, reads_path, aln_reads_to_raw_contig)
	# execute_command(command, sys.stderr, dry_run=False);
	# Ovdje naci overlap

	### Touch the file to clear it.
	fp = open(all_corrected_reads_path, 'w');
	fp.close();

	for i in xrange(0, len(contig.gp)):
		timestamp = strftime("%Y/%m/%d %H:%M:%S", gmtime());
		sys.stderr.write('[%s] Correcting read %d/%d.\n' % (timestamp, (i + 1), len(contig.gp)));
		# if (i < 3): continue;

		gp = contig.gp[i];

		[h, s, q] = seq_hash[gp.read_name];

		fp = open(single_read_path, 'w');
		fp.write('@%s\n%s\n+\n%s\n' % (h, s, q));
		fp.close();

		command = '%s align -r %s -d %s -o %s --no-self-hits --mapq -1 --rebuild-index -v 1' % (graphmap_bin, single_read_path, reads_path, single_read_sam)
		execute_command(command, sys.stderr, dry_run=False);

		# command = '%s/../bin/consise --align 1 -M 5 -X -4 -G -8 -E -6 -w 500 --bq 10.0 --ovl-margin 0.0 --msa poa -b 200 -t 4 %s %s %s' % (SCRIPT_PATH, single_read_path, single_read_sam, single_read_consensus);
		command = '%s/../bin/consise --align 1 -M 1 -X -1 -G -1 -E -1 -w 500 --bq 10.0 --ovl-margin 0.0 --msa poa -b 200 -t %d %s %s %s' % (SCRIPT_PATH, num_threads, single_read_path, single_read_sam, single_read_consensus);
		# command = '%s/../bin/consise --align 1 -M 5 -X -4 -G -8 -E -6 -w 500 --bq 10.0 --ovl-margin 0.0 --msa poa -b 200 -t %d %s %s %s' % (SCRIPT_PATH, num_threads, single_read_path, single_read_sam, single_read_consensus);
		execute_command(command, sys.stderr, dry_run=False);

		if (ref_path != None and ref_path != '-'):
			command = 'dnadiff -p %s/../temp/dnadiff/consread %s %s;' % (SCRIPT_PATH, ref_path, single_read_consensus);
			command += 'grep "TotalBases" %s/../temp/dnadiff/consread.report;' % (SCRIPT_PATH);
			command += 'grep "AlignedBases" %s/../temp/dnadiff/consread.report;' % (SCRIPT_PATH);
			command += 'grep "AvgIdentity" %s/../temp/dnadiff/consread.report;' % (SCRIPT_PATH);
			execute_command(command, sys.stderr, dry_run=False);

		command = 'cat %s >> %s' % (single_read_consensus, all_corrected_reads_path);
		execute_command(command, sys.stderr, dry_run=False);

	fp = open(raw_contig_temp_path, 'w');
	fp.write('>%s\n%s\n' % (contig.name, contig.seq));
	fp.close();

	command = '%s align -r %s -d %s -o %s --no-self-hits --mapq -1 --rebuild-index -v 1' % (graphmap_bin, raw_contig_temp_path, all_corrected_reads_path, all_corrected_reads_sam)
	execute_command(command, sys.stderr, dry_run=False);

	command = '%s/../bin/consise --align 1 -M 5 -X -4 -G -8 -E -6 -w 500 --ovl-margin 0.0 --msa poa -b 200 -t %d %s %s %s' % (SCRIPT_PATH, num_threads, raw_contig_temp_path, all_corrected_reads_sam, final_consensus_contig_path);
	execute_command(command, sys.stderr, dry_run=False);

	command = 'cat %s >> %s' % (final_consensus_contig_path, out_consensus_path);
	execute_command(command, sys.stderr, dry_run=False);

	if (ref_path != None and ref_path != '-'):
		command = 'dnadiff -p %s/../temp/dnadiff/cons %s %s;' % (SCRIPT_PATH, ref_path, final_consensus_contig_path);
		command += 'grep "TotalBases" %s/../temp/dnadiff/cons.report;' % (SCRIPT_PATH);
		command += 'grep "AlignedBases" %s/../temp/dnadiff/cons.report;' % (SCRIPT_PATH);
		command += 'grep "AvgIdentity" %s/../temp/dnadiff/cons.report;' % (SCRIPT_PATH);
		execute_command(command, sys.stderr, dry_run=False);



def main():
	if (len(sys.argv) < 5 or len(sys.argv) > 6):
		sys.stderr.write('Script for producing a consensus from GFA assembly layouts.\n');
		sys.stderr.write('Usage:\n');
		sys.stderr.write('  %s <layout.gfa> <reads.fastq> <consensus.fasta> num_threads [<reference.fa>]\n' % (sys.argv[0]));
		exit(1);

	gfa_file = sys.argv[1];
	reads_file = sys.argv[2];
	out_cons_file = sys.argv[3];
	num_threads = int(sys.argv[4]);
	reference_path = None;
	if (len(sys.argv) >= 6):
		reference_path = sys.argv[5];

	contigs = parse_gfa(gfa_file);
	verbose_gfa_contigs(contigs);

	fp = open(out_cons_file, 'w');
	fp.close();
	for contig in contigs:
		correct_layout_reads(contig, reads_file, num_threads, out_cons_file, reference_path);
if __name__ == "__main__":
	main();
