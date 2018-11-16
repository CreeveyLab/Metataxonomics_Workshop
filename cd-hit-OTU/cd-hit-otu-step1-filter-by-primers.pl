#!/usr/bin/perl

my $fasta = shift;
my $primer_file = shift;

if (not defined($primer_file)) {

   print <<EOD;

Usage $0 fasta-file-of-raw-reads primer_file

This script read a file of primer sequences,
scan a fasta file of raw reads, and remove the reads that don't match the primers

format of primer_file
AGTGCGTAGTG[ACTG]CAGC[AC]GCCGCGGTAA

EOD
   die;
}

my @primers = ();
open(TMP, $primer_file) || die "can not open $primer_file\n";
while($ll=<TMP>){
  $ll =~ s/\s//g;
  push(@primers, $ll);
}
close(TMP);

my $seq_total = 0;
my $seq_good = 0;
my $seq_bad = 0;

open(TMP, $fasta) || die "can not open $fasta\n";
while($ll=<TMP>){
  if ($ll =~ /^>/){
    if ($seq) {
      if (match_primers($seq)) {
        print $des, $seq, "\n";
        $seq_good++; $seq_total++;
      }
      else {
        $seq_bad++; $seq_total++;
      }
    }
    $des = $ll;
    $seq = "";
  }
  else {
    $ll =~ s/\s//g;
    $seq = $seq . $ll;
  } 
}
close(TMP);


print <<EOD;
Total seqs:	$seq_total
Good seqs:	$seq_good
Filtered seqs:	$seq_bad
EOD

sub match_primers {
  my $match_flag = 0;
  my $seq = shift;
  my ($i, $j, $k);

  foreach $i (@primers) {
    if ($seq =~ /^$i/i) {
      $match_flag = 1;
      last;
    }
  }
  return $match_flag;
}
