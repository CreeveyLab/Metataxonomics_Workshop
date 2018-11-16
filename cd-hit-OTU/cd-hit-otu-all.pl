#!/usr/bin/perl

use Getopt::Std;
my $script_name = $0;
my $script_dir = $0;
   $script_dir =~ s/[^\/]+$//;
   chop($script_dir);
   $script_dir = "./" unless ($script_dir);

getopts("i:o:p:c:e:s:f:m:",\%opts);
die usage() unless ($opts{i} and $opts{o});

my $fasta_raw  = $opts{i};
my $dir        = $opts{o};
my $prefix     = $opts{p}; $prefix     = 6      unless ($prefix);
my $cutoff_p   = $opts{f}; $cutoff_p   = 0.8    unless ($cutoff_p);
my $err        = $opts{e}; $err        = 0.005  unless ($err);
my $sd1        = $opts{s}; $sd1        = 1      unless ($sd1);
my $otu_cutoff = $opts{c}; $otu_cutoff = 0.97   unless ($otu_cutoff);
my $chimera_f  = $opts{m}; $chimera_f  = "true" unless ($chimera_f);

my ($str, $cmd, $ll);
$cmd = `mkdir -p $dir`;
my $f2 = "$dir/fltd";


#### Trim
$str = "$script_dir/cd-hit-otu-step1-filter-trim.pl -i $fasta_raw -o $f2 -f $cutoff_p -p $prefix > $dir/trim.log";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

#### dup and chimaric checking
$str = "$script_dir/cdhit-dup/cd-hit-dup -i $f2 -o $f2.dup -d 0 -m false -f $chimera_f > $dir/dup.log";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

my $clstr = "$f2.dup.clstr";
my $clstr2 = "$f2.dup2.clstr";

my $largest_chimaric = 0;
my $largest = 0;
my $ave_len = 0;
my ($i, $j, $k);

  ($i, $largest, $largest_rep, $largest_faa, $ave_len) = clstr_info0($clstr, "$f2.dup");
  $homo_sites = check_homo_sites($largest_faa, 3);
  $ave_len = length($largest_faa);

  $str = "$script_dir/cd-hit-otu-simulation-pyro.pl $err $ave_len $homo_sites $largest $sd1";

print "\n\nrunning command\n$str\n";
$cmd = `$str`;
print $cmd;
open(OUT1, "> $dir/simulation.log") || die "can not open to write";
print OUT1 $cmd;
close(OUT1);
## simulation

my $cutoff = 0;
my @lls = split(/\n/,$cmd);
foreach $ll (@lls) {
  if ($ll =~ /Suggested cutoff value to remove small cluster is (\d+)/) {
    $cutoff = $1; last;
  }
}
$cutoff++;
die "can not parse cluster cutoff" unless ($cutoff);


## combine and sort (regular + chimaric) clusters
my $t1 = "$f2.dupall";
$str = "cat $clstr $clstr2 | $script_dir/cd-hit/clstr_sort_by.pl > $t1.clstr";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

## select reps from combined clusters
$str="$script_dir/clstr_sort_trim_rep.pl $t1.clstr $f2 2 > $t1-rep.fa";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

## cluster at $c1 => 99.5 or something that allow 1 mismatch or 1-2 gaps
my $c1 = ($ave_len-1) / $ave_len - 0.001;
$str = "$script_dir/cd-hit/cd-hit-est -i $t1-rep.fa -o $t1.nr2nd -c $c1 -n 10 -l 11 -p 1 -d 0 -g 1 -b 3 > $t1.nr2nd.log";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

$str = "$script_dir/cd-hit/clstr_rev.pl $t1.clstr $t1.nr2nd.clstr | $script_dir/cd-hit/clstr_sort_by.pl > $t1.nr2nd-all.clstr";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;

## select reps from $t1.nr2nd-all.clstr
$str = "$script_dir/clstr_select_rep.pl size $cutoff 999999999 < $t1.nr2nd-all.clstr > $t1-pri-rep.ids";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;
print $cmd;

$str = "$script_dir/fetch_fasta_by_ids.pl $t1-pri-rep.ids $t1-rep.fa > $t1-pri-rep.fa";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;
print $cmd;

if (-s $clstr2) {
  ## save chimaric ids
  $str = "$script_dir/clstr_select_rep.pl size 1 999999999 < $clstr2 > $dir/chimaric.ids";
  print "\n\nrunning command\n$str\n";
  $cmd = `$str`;

  ## exclude chimaric reads from $t1-pri-rep.fa
  $str = "$script_dir/fetch_fasta_exclude_ids.pl $dir/chimaric.ids $t1-pri-rep.fa > $t1-pri-rep-good.fa";
  print "\n\nrunning command\n$str\n";
  $cmd = `$str`;
  print $cmd;

  $str = "$script_dir/cd-hit/cd-hit-est -i $t1-pri-rep-good.fa -o $dir/OTU -c $otu_cutoff -n 8 -l 11 -p 1 -d 0 -g 1 -b 5 -G 0 -aS 0.8 > $dir/OTU.log";
}
else {
  $str = "$script_dir/cd-hit/cd-hit-est -i $t1-pri-rep.fa      -o $dir/OTU -c $otu_cutoff -n 8 -l 11 -p 1 -d 0 -g 1 -b 5 -G 0 -aS 0.8 > $dir/OTU.log";
}

print "\n\nrunning command\n$str\n";
$cmd = `$str`;
print $cmd;

$str = "$script_dir/cd-hit/clstr_rev.pl $t1.nr2nd-all.clstr $dir/OTU.clstr > $dir/OTU.nr2nd.clstr";
print "\n\nrunning command\n$str\n";
$cmd = `$str`;
print $cmd;

my ($tu,$ts,$cu,$cs)=times(); my $tt=$tu+$ts+$cu+$cs;
print "\ntotal cpu time $tt\n";
$cmd = `echo >> $dir/OTU.log`;
$cmd = `echo total cpu time $tt >> $dir/OTU.log`;


sub clstr_info1 {
  my $clstr_file = shift;

  my $largest = 0;
  my $N = 0;
  my $total_len = 0;
  my $ll;
  my $no = 0;
  open(TMP, $clstr_file) || die "can not open $clstr_file";
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($no) {
        if ($no > $largest) {$largest=$no;}
      }
      $no = 0;
    }
    else {
      $N++;
      $no++;
      chop($ll);
      if ($ll =~ /(\d+)(aa|nt),/) {
        $total_len += $1;
      }
    }
  }
      if ($no) {
        if ($no > $largest) {$largest=$no;}
      }
  close(TMP);
  my $ave_len = int($total_len/$N);
  return($N, $largest, $ave_len);
}
########## END clstr_info1


sub usage {
<<EOF
Usage:

$script_name -i input_fasta_file -o output_dir -f length-cutoff-lower-fraction  -p prefix-length/primers_file -e sequence_error_rate -c OTU_cutoff

Parameters:
-i input fasta file 
-o output dir
-f length-cutoff-lower-fraction default 0.8
-p prefix-length/primers_file default 6
-e sequence error rate default 0.005
-c OTU cutoff default 0.97
-m whether to perform chimera checking (true/false), default true

length-cutoff-lower-fraction & prefix-length/primers_file are used in filtering and trimming raw reads, where
        if a primers_file is provided,
          read primers from this file, remove the reads don't match the primers
        if a prefix-length (a digit number) is provided
          get the consensus of prefix of the all reads
          remove the reads without this prefix
        calculate the median length of reads, this is the upper length cutoff.
          Long reads are trimmed to this length.
          This cutoff * lenght-cutoff-lower-fraction is the lower length cutoff.
          Short reads below this length are removed.

format of primer_file, each line is a primer sequence, where [AT] mean either A or T
AGTGCGTAGTG[ACTG]CAGC[AC]GCCGCGGTAA


EOF
}
###### END usage


sub clstr_info0 {
  my $clstr_file = shift;
  my $fasta_file = shift;

  my $largest = 0;
  my $N = 0;
  my $total_len = 0;
  my $ll;
  my $no = 0;
  my $rep = "";
  my $largest_rep = "";
  open(TMP, $clstr_file) || die "can not open $clstr_file";
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($no) {
        if ($no > $largest) {$largest=$no; $largest_rep=$rep;}
      }
      $no = 0;
    }
    else {
      $N++;
      $no++;
      chop($ll);
      if ($ll =~ /(\d+)(aa|nt),/) {
        $total_len += $1;
      }

      if ($ll =~ /\*$/) {
        if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
          $rep = $3;
        }
      }
    }
  }
  close(TMP);
  my $ave_len = int($total_len/$N);
  my $largest_faa = "";

  open(TMP, $fasta_file) || die "can not open $fasta_file";
  while($ll=<TMP>){
    if ($ll =~ /^>/) {
      if ($ll =~ /^>(\S+)/) {
        my $id = $1;
        if ($id eq $largest_rep) {
          while($ll=<TMP>){
            last if ($ll =~ /^>/);
            $ll =~ s/\s//g;
            $largest_faa .= $ll;
          }
          last;
        }
      }
    }
  }
  close(TMP);
  return($N, $largest, $largest_rep, $largest_faa, $ave_len);
}
########## END clstr_info0

sub check_homo_sites {
  my $seq = shift;
  my $homo_cutoff=shift;;
     $homo_cutoff = 3 unless ($homo_cutoff);

  $seq =~ s/\s//g;
  my ($i, $j, $k);
  my $homo_no=0;
  my $len = length($seq);

  for ($i=0; $i<$len-1; $i++){
    my $c=substr($seq,$i,1);

    for ($k=0, $j=$i+1; $j<$len; $j++) {
      my $c1=substr($seq,$j,1);
      last if ($c1 ne $c);
      $k++;
    }

    if ($k >=$homo_cutoff-1 ) {
      $homo_no++;
      $i+=$k;
    }
  }
  return $homo_no;
}
################ END check_homo_sites
























