#!/usr/bin/perl

use Getopt::Std;

getopts("i:o:p:f:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $fasta      = $opts{i};
my $output     = $opts{o};
my $prefix     = $opts{p};
my $cutoff_len_p = $opts{f}; $cutoff_len_p = 0.8   unless ($cutoff_len_p);

my $primer_file = "";
if (-e $prefix) { $primer_file = $prefix; $prefix=0;}
else { $prefix = 6 unless ($prefix) };


my @primers = ();
if ($primer_file) {
  open(TMP, $primer_file) || die "can not open $primer_file\n";
  while($ll=<TMP>){
    $ll =~ s/\s//g;
    push(@primers, $ll);
  }
  close(TMP);
}

my $n=0;
my @len = ();
my $seq = "";
my %head = ();
my %len_dist = ();

open(TMP, $fasta) || die "can not open $fasta\n";
while($ll=<TMP>){
  if ($ll =~ /^>/){
    if ($seq) {
      $n++;
      push(@len, length($seq));
      $len_dist{length($seq)}++;
      if ($prefix >0) {
        my $head = uc(substr($seq,0,$prefix));
        $head{$head}++;
      }
    }
    $seq = "";
  }
  else {
    $ll =~ s/\s//g;
    $seq = $seq . $ll;
  }  
}
    if ($seq) {
      $n++;
      push(@len, length($seq));
      $len_dist{length($seq)}++;
      if ($prefix >0) {
        my $head = uc(substr($seq,0,$prefix));
        $head{$head}++;
      }
    }
close(TMP);
my @head_key;
my $consensus;
if ($prefix >0) {
  @head_key = keys %head;
  @head_key = sort {$head{$b} <=> $head{$a}} @head_key;
  #foreach $i (@head_key) { print "$i\t$head{$i}\n"; }
  $consensus = $head_key[0];
}

my @len_key = keys %len_dist;
@len_key = sort {$len_dist{$b} <=> $len_dist{$a}} @len_key;
my $lend = $len_key[0];

my $len_t = 0;
for ($i=0; $i<$n; $i++) {
  $len_t += $len[$i];
}
my $len_ave = $len_t/$n;
my $len_med = $len[int($n/2)];

my $dev = 0;
for ($i=0; $i<$n; $i++){
  $dev += ($len[$i]-$len_ave) * ($len[$i]-$len_ave);
}
$dev = sqrt($dev/$n);

my $len_cutoff_upper = $len_med;
my $len_cutoff_lower = int($len_med * $cutoff_len_p);
my $des = "";
my $len1;
$seq = "";

my $trimmed_letter = 0;
my $trimmed_seq = 0;
my $seq_wrong_head = 0;
my $seq_short = 0;
my $seq_N = 0;
my $seq_good = 0;

open(TMP, $fasta) || die "can not open $fasta\n";
open(OUT, ">$output") || die "can not write to $output\n";
while($ll=<TMP>){
  if ($ll =~ /^>/){
    if ($seq) {
      $len1 = length($seq);
      if ($len1 > $len_cutoff_upper) {
        $seq = substr($seq,0,$len_cutoff_upper);
        $trimmed_seq++;
        $trimmed_letter += $len1 - $len_cutoff_upper;
      }
      if    (($primer_file and (not match_primers($seq))) or 
             (($prifix >0) and (substr($seq,0,$prefix) ne $consensus)) ) {$seq_wrong_head++;}
      elsif ($seq =~ /N|n/) {$seq_N++;}
      elsif ($len1 < $len_cutoff_lower) {$seq_short++;}
      else { print OUT $des,$seq,"\n"; $seq_good++; }      
    }
    $des = $ll;
    $seq = "";
  }
  else {
    $ll =~ s/\s//g;
    $seq = $seq . $ll;
  } 
}
    if ($seq) {
      $len1 = length($seq);
      if ($len1 > $len_cutoff_upper) {
        $seq = substr($seq,0,$len_cutoff_upper);
        $trimmed_seq++;
        $trimmed_letter += $len1 - $len_cutoff_upper;
      }
      if    (($primer_file and (not match_primers($seq))) or 
             (($prifix >0) and (substr($seq,0,$prefix) ne $consensus)) ) {$seq_wrong_head++;}
      elsif ($seq =~ /N|n/) {$seq_N++;}
      elsif ($len1 < $len_cutoff_lower) {$seq_short++;}
      else { print OUT $des,$seq,"\n"; $seq_good++; }      
    }
close(TMP);
close(OUT);



$percent_trimmed_letter = int(100000*$trimmed_letter/$len_t)/1000;
$percent_trimmed_seq    = int(100000*$trimmed_seq/$n)/       1000;

print <<EOD;

Total seq:	$n
Total letters:	$len_t
Average_len:	$len_ave
Med_len:	$len_med
Standard dev:	$dev

trimmed_letters:	$trimmed_letter
trimmed_letters(%):	$percent_trimmed_letter
trimmed_seqs:		$trimmed_seq
trimmed_seqs(%):	$percent_trimmed_seq

prefix consensus: $consensus
filtered seqs with wrong prefix or primers: $seq_wrong_head
filtered seqs with ambiguous base calls: $seq_N
long seqs trimmed to: $len_cutoff_upper
short seqs < $len_cutoff_lower: $seq_short
good seqs left: $seq_good
EOD

if ($prefix>0) {
  print "\n\nPrefix\tcount\n";
  foreach $i (@head_key) { print "$i\t$head{$i}\n"; }
}
print "\n\nLength\tcount\n";
@len_key = sort {$b <=> $a} @len_key;
foreach $i (@len_key) { print "$i\t$len_dist{$i}\n";}



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





sub usage {
<<EOD

Usage $0 -i fasta-file-of-raw-reads -o output -f length-cutoff-lower-fraction  -p prefix-length/primers_file

        lenght-cutoff-lower-fraction default 0.8
        prefix-length default 6

This script scan a fasta file of raw reads
        (1) depends on third parameter 
            (1a) if a primers_file is provided, 
                 read primers from this file, remove the reads don't match the primers
            (1b) if a prefix-length (a digit number) is provided
	         get the consensus of prefix of the all reads 
                 remove the reads without this prefix
        (2) reads with letter 'N' are removed
        (3) calculate the median length of reads, this is the upper length cutoff.
            Long reads are trimmed to this length.
            This cutoff * lenght-cutoff-lower-fraction is the lower length cutoff.
            Short reads below this length are removed.

format of primer_file, each line is a primer sequence, where [AT] mean either A or T
AGTGCGTAGTG[ACTG]CAGC[AC]GCCGCGGTAA

EOD
}
######### END usage
