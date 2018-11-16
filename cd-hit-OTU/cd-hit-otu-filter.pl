#!/usr/bin/perl

use Getopt::Std;

getopts("i:o:p:t:P:",\%opts);
die usage() unless ($opts{i} and $opts{o} and defined($opts{P}));

my $pair_flag    = $opts{P};

my $fastq        = $opts{i};
my $output       = $opts{o};
my $prefix       = $opts{p};
my $cutoff_len_p = $opts{t}; 

my $fastq2     = "";
my $output2    = "";
my $prefix2    = "";
my $cutoff_len_p2 = "";

my @primers = ();
my @primers2= ();

###### parse parameters
###### read prefix or primers
if ($pair_flag) {
  ($fastq,        $fastq2)        = split(/\,/, $fastq);
  ($output,       $output2)       = split(/\,/, $output);
  ($prefix,       $prefix2)       = split(/\,/, $prefix);
  ($cutoff_len_p, $cutoff_len_p2) = split(/\,/, $cutoff_len_p);

  $cutoff_len_p2 = 1.0   unless ($cutoff_len_p2);

  if (-e $prefix2) { 
    @primers2  = read_primer_file($prefix2);
  }
  else { 
    $prefix2 = 6 unless (defined($prefix2)); 
    @primers2 = read_consensus($fastq2, $prefix2);
  }
} ## END if ($pair_flag)

$cutoff_len_p = 1.0   unless ($cutoff_len_p);
if (-e $prefix) { 
  @primers  = read_primer_file($prefix);
}
else { 
  $prefix = 6 unless (defined($prefix)); 
  @primers = read_consensus($fastq, $prefix);
}

my ($score_v, $offset)=check_fastq_version($fastq);
my $e_cutoff  = 0.1; #worst base cutoff for error probability
my $q_cutoff  = int(10 * log(1/$e_cutoff) / log(10)); #worst base score
my $q_cutoffn = 20; ## binomial thought that if p=0.1 and n=20, then it is likely (p=0.8) to have > 2 errors
my $q_cutoffnp= 30; ## since we will do pair correction, we tolerant more errors for paired reads


my ($ll, $lla, $llb, $id, $ida, $idb, $seq, $seqa, $seqb, $qua, $quaa, $quab);
my ($len, $lena, $lenb);
my ($i, $j, $k);
my $seq_no = 0;
my $seq_N  = 0;
my $seq_low_q = 0;
my $seq_non_prefix = 0;
my $seq_good = 0;
my @score_array  = ();
my @score_array2 = ();
my $score_count  = 0;

######## scan files
if ($pair_flag) {
  open(TTTa, $fastq)  || die "can not open $fastq\n";
  open(TTTb, $fastq2) || die "can not open $fastq2\n";
  open(OUTa, "> $output")  || die "can not write $output\n";
  open(OUTb, "> $output2") || die "can not write $output2\n";
  while(($lla=<TTTa>) and ($llb=<TTTb>)){
    chop($lla); $lla =~ s/\s.+$//; $ida = substr($lla,1);
    chop($llb); $llb =~ s/\s.+$//; $idb = substr($llb,1);
    die "$ida from $fastq doesn't match $idb from $fastq2" unless ($ida eq $idb);

    $seqa = <TTTa>; $seqa =~ s/\s+$//g;
    $seqb = <TTTb>; $seqb =~ s/\s+$//g;
    $lla  = <TTTa>; #read ID
    $llb  = <TTTb>; #read ID
    $quaa = <TTTa>; $quaa =~ s/\s+$//g;
    $quab = <TTTb>; $quab =~ s/\s+$//g;
    $lena = length($seqa);
    $lenb = length($seqb);

    $seq_no++;

    if (($seqa =~ /N|n/) or ($seqb =~ /N|n/)) {
      $seq_N++;
      next;
    }

    if (@primers) {
      if (not match_primers($seqa, @primers )) { $seq_non_prefix++; next; }
      if (not match_primers($seqb, @primers2)) { $seq_non_prefix++; next; }
    }

    my $cutoffa = ($cutoff_len_p < 1.001) ? int($lena * $cutoff_len_p ) : $cutoff_len_p;
    my $cutoffb = ($cutoff_len_p2< 1.001) ? int($lenb * $cutoff_len_p2) : $cutoff_len_p2;

    for ($i=0, $j=0; $i<$cutoffa; $i++) {
      my $c1=ord(substr($quaa,$i,1))-$offset;
      $j++ if ($c1 < $q_cutoff);
    }
    if ($j > $q_cutoffnp) {$seq_low_q++;  next;}

    for ($i=0, $j=0; $i<$cutoffb; $i++) {
      my $c1=ord(substr($quab,$i,1))-$offset;
      $j++ if ($c1 < $q_cutoff);
    }
    if ($j > $q_cutoffnp) {$seq_low_q++;  next;}

    if ($cutoffa < $lena) { $seqa = substr($seqa, 0, $cutoffa); $quaa = substr($quaa, 0, $cutoffa); $lena = $cutoffa; }
    if ($cutoffb < $lenb) { $seqb = substr($seqb, 0, $cutoffb); $quab = substr($quab, 0, $cutoffb); $lenb = $cutoffb; }

    print OUTa ">$ida\n$seqa\n";
    print OUTb ">$idb\n$seqb\n";

    if (rand() < 1/16 ) {
      for ($i=0; $i<$lena; $i++) { my $c1 = ord(substr($quaa,$i,1))-$offset; $score_array[$i]  += $c1; }
      for ($i=0; $i<$lenb; $i++) { my $c1 = ord(substr($quab,$i,1))-$offset; $score_array2[$i] += $c1; }
      $score_count++;
    }
    $seq_good++;
  }
  close(TTTa);
  close(TTTb);
  close(OUTa);
  close(OUTb);
  process_score("$output.err", $score_count, @score_array);
  process_score("$output2.err", $score_count, @score_array2);

} #### END ($pair_flag)
else {
  open(TTTa, $fastq)  || die "can not open $fastq\n";
  open(OUTa, "> $output")  || die "can not write $output\n";
  while($lla=<TTTa>){
    chop($lla); $lla =~ s/\s.+$//; $ida = substr($lla,1);

    $seqa = <TTTa>; $seqa =~ s/\s+$//g;
    $lla  = <TTTa>; #read ID
    $quaa = <TTTa>; $quaa =~ s/\s+$//g;
    $lena = length($seqa);

    $seq_no++;

    if ($seqa =~ /N|n/) { $seq_N++; next; }

    if (@primers) {
      if (not match_primers($seqa, @primers )) { $seq_non_prefix++; next; }
    }

    my $cutoffa = ($cutoff_len_p < 1.001) ? int($lena * $cutoff_len_p ) : $cutoff_len_p;

    for ($i=0, $j=0; $i<$cutoffa; $i++) {
      my $c1=ord(substr($quaa,$i,1))-$offset;
      $j++ if ($c1 < $q_cutoff);
    }
    if ($j > $q_cutoffn) {$seq_low_q++;  next;}

    if ($cutoffa < $lena) { $seqa = substr($seqa, 0, $cutoffa); $quaa = substr($quaa, 0, $cutoffa); $lena = $cutoffa; }

    print OUTa ">$ida\n$seqa\n";

    if (rand() < 1/16 ) {
      for ($i=0; $i<$lena; $i++) { my $c1 = ord(substr($quaa,$i,1))-$offset; $score_array[$i]  += $c1; }
      $score_count++;
    }
    $seq_good++;
  }
  close(TTTa);
  close(OUTa);

  process_score("$output.err", $score_count, @score_array);
} #### END else


my $primer_line = join("/", @primers);
print <<EOD;
quality score version: $score_v
quality score offset: $offset
Total seq:	$seq_no
filtered seqs with ambiguous base calls: $seq_N
filtered seqs with wrong prefix or primers: $seq_non_prefix
      primers: $primer_line
filtered seqs with low-quality bases: $seq_low_q
good seqs left: $seq_good
EOD


sub match_primers {
  my ($seq, @primers) = @_;
  my $match_flag = 0;
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
This script filter and trim fastq files
        (1) reads with letter 'N' are removed
        (2) trim long sequences, see option -t
        (3) tilter based on prefix or primer sequences, see option -p
        (4) scan quality scores of filtered reads can calculate error probability
        (5) output are in fasta format
it also generates average (averaged from good sequences) error probability for each base, see .err files

Usage $0 optioins
  for single reads
         -P 0
         -i fastq-file-of-raw-reads
         -o output-fasta-file
         -t trim_cutoff, default 1.0 (means no trimming)
            if cutoff is a integer number > 1 (like 200), the program will trim sequences to this length
            if cutoff is a fraction (like 0.8), the program will keep fraction of this reads
         -p prefix-length/primers_file
            (a) if a primers_file is provided, 
                 read primers from this file, remove the reads don't match the primers
            (b) if a prefix-length (a digit number) is provided, default 6
                 get the consensus of prefix of the all reads 
                 remove the reads without this prefix

  for paired reads
         -P 1
         -i fastq-file-of-raw-reads-of-1-end,fastq-file-of-raw-reads-of-another-end
         -o output-fasta-file-of-1-end,output-fasta-file-of-another-end
         -t trim_cutoff-of-1-end,trim_cutoff-of-another-end default 1.0/1.0 (means no trimming)
            see description for single reads
         -p prefix-length-of-1-end,prefix-length-of-another-end/primers-file-of-1-end,primers-file-of-another-end
            see description for single reads


EOD
}
######### END usage










sub read_primer_file {
  my $file = shift;
  my @p = ();
  my $ll;

  open(PPP, $file) || die "can not open $file";
  while($ll=<PPP>){
    $ll =~ s/\s//g;
    push(@p, $ll);
  }
  close(PPP);
  return @p;
}
######## END read_primer_file


sub read_consensus {
  my $file = shift;
  my $prefix = shift;
  my @p = ();
  return @p if ($prefix <=0);

  my $n = 0;
  my $ll;
  my %prefix_count = ();
  my ($i, $j, $k);

  open(QQQ, $file) || die "can not open $file\n";
  while($ll=<QQQ>){
    $ll = <QQQ>;
    my $p1 = uc(substr($ll,0,$prefix));
    $prefix_count{$p1}++;
    $n++;
    $ll = <QQQ>;
    $ll = <QQQ>;
  }
  close(QQQ);

  my @prefix_key = keys %prefix_count;;
  @prefix_key = sort {$prefix_count{$b} <=> $prefix_count{$a}} @prefix_key;

  $j = 0;
  $k = 0;
  foreach $i (@prefix_key){
    push(@p, $i);
    $j += $prefix_count{$i};
    $k ++;
    last if ($j > $n * 0.8); #keep consensus that make up 80% of reads
    last if ($k >= 4);        #most 4 consensus
  }

  return @p;
}
########## END read_consensus



sub process_score {
  my ($output, $score_count, @scores) = @_;
  my $len = $#scores+1;
  my ($i, $j, $k);
  open(SSS, "> $output") || die "can not write $output\n";

  for ($i=0; $i<$len; $i++){
    my $ave = int($scores[$i]/$score_count);
    my $e = 1 / (10 ** ($ave/10));
    print SSS "$i\t$ave\t$e\n";
  }
  close(SSS);
}
########## END process score


sub check_fastq_version {
  my $f=shift;
  my ($i, $j, $k, $ll);

  $i=0;
  $j=0;
  open(TMP, $f) || die "can not open $f\n";
  while($ll=<TMP>){
    if ($ll =~ /^\@/) {
      $ll = <TMP>; #### read sequence
      $ll = <TMP>; #### read ID
      $ll = <TMP>; #### read quality score
      $j += ord(substr($ll,0,1)); # just check first base
      $i++;
      last if ($i>1000); #1000 reads is enough
    }
  }
  close(TMP);

  $k = $j/$i;
  my $offset = ($k>80) ? 64 : 33;
  my $v      = ($k>80) ? "Illumina 1.3/1.5" : "Illumina 1.8";
  return ($v, $offset);
}
#### END check_fastq_version
