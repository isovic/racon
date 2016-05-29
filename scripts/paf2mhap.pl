#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;

my %opts = ();
getopts("p", \%opts);
my $is_100 = defined($opts{p});

die("Usage: paf2mhap.pl [-p] <ref.fa> <in.fa> <in.paf>\n") if (@ARGV == 0);

warn("Parsing ref FASTA to create the name<=>id table...\n");
my %hashref;
my $fnref = shift(@ARGV);
open(FHref, $fnref =~ /\.gz$/? "gzip -dc {} |" : $fnref) || die;
my $cntref = 0;
while (<FHref>) {
	if (/^>(\S+)/) {
		$hashref{$1} = ++$cntref unless defined($hashref{$1});
	}
}
close(FHref);

warn("Parsing reads FASTA to create the name<=>id table...\n");
my %hash;
my $fn = shift(@ARGV);
open(FH, $fn =~ /\.gz$/? "gzip -dc {} |" : $fn) || die;
my $cnt = 0;
while (<FH>) {
	if (/^>(\S+)/) {
		$hash{$1} = ++$cnt unless defined($hash{$1});
	}
}
close(FH);

warn("Converting PAF to MHAP format...\n");
while (<>) {
	chomp;
	my @t = split;
	next if ($t[0] eq $t[5]); # NB: ignore self matches
	my $cnt = /cm:i:(\d+)/? $1 : 0;
	my $r = $t[9] / $t[10];
	$r = sprintf("%.4f", $is_100? 100. * $r : $r);
	die if !defined($hash{$t[0]}) || !defined($hashref{$t[5]});
	print(join(" ", $hash{$t[0]}, $hashref{$t[5]}, $r, $cnt, 0, @t[2,3,1], $t[4] eq '+'? 0 : 1, @t[7,8,6]), "\n");
}
