#!/usr/bin/perl

my $len= $ARGV[0];
my $err= $ARGV[1];
my $obv_pri = $ARGV[2];
my $sd1 = $ARGV[3];
   $sd1 = 5 unless (defined($sd1));

#######################
#### usage ############
#######################
if (!$len or !$err or !$obv_pri) {
  print <<EOD;


Simulation to select the cluster size cutoff...
Please give read length, error rate per base and size of the largest cluster

Usage: $0 read-length per-base-error-rate size-of-the-largest-cluster
	read-length is the average length of tags
	per-base-error-rate is the error rate, the 454 pyrosequencing error rate is about 0.005
		if the low-quality reads were filtered out, then the error rate is about 0.0025
	size-of-the-largest-cluster is the size of the largest cluster at 1st cluster run
		by cd-hit-454 at 100% with parameters like "-c 1.0 -n 10 -b 1"

Program developed by Weizhong Li's lab: http://weizhongli-lab.org/

EOD
  exit;
  ##### exit
}


my ($i, $j, $k);
my %error_model = (
  'x' => 0.30,  # fraction of deletion error
                # suppose correct base is G
  'A' => 2/30,  # fraction of substitution error to base A, about 0.0667
  'C' => 2/30,  # fraction of substitution error to base C, about 0.0667
  'T' => 2/30,  # fraction of substitution error to base T, about 0.0667
                # suppose correct base is G
  'a' => 2/30,  # fraction of insertion error to base A, about 0.0667
  'c' => 2/30,  # fraction of insertion error to base C, about 0.0667
  't' => 2/30,  # fraction of insertion error to base T, about 0.0667
  'g' => 0.30,  # fraction of insertion error to base G, 
                # Homopolymer extension, represents big portion of insertion error
  # totoal sum is 1.0
  # error rate is based on
  # S. M. Huse, J. A. Huber, H. G. Morrison et al., Genome Biol 8 (7), R143 (2007).
  # B. Niu, L. Fu, S. Sun et al., BMC Bioinformatics 11, 187 (2010).
);

my @error_model_1000=();
my $error_model_N = 0;
foreach $i (keys %error_model) {
  $j = int($error_model{$i} * 1000);
  for ($k=0; $k<$j; $k++) {push(@error_model_1000,$i);}
}
$error_model_N = $#error_model_1000+1;


my $N = 50;  ####

my $cov1 = $obv_pri;
my @min_array=();
my @max_array=();

srand(124);
print STDERR "simulating ";

while(1){
  my ($ave_cls, $sd_cls, $ave_pri, $ave_sec, $ave_secb, $sd_pri, $sd_sec, $sd_secb) = sim1($len, $err, $cov1, $N);

  #print "$cov1\t$ave_pri\t$sd_pri\t$ave_sec\t$sd_sec\n";

  if ($ave_pri+$sd_pri*$sd1 < $obv_pri) {
    @min_array=($cov1, $ave_pri, $ave_secb, $sd_pri, $sd_secb);
  }
  if ($ave_pri-$sd_pri*$sd1 > $obv_pri) {
    @max_array=($cov1, $ave_pri, $ave_secb, $sd_pri, $sd_secb);
    last;
  }
  $cov1++;
  $cov1 = int($cov1*1.05);
  print STDERR ".";
}

#print "min\t", join("\t", @min_array), "\n";
#print "min\t", join("\t", @max_array), "\n";

my $max_sec = int($max_array[2] + $sd1*$max_array[4])+1;
my $cutoff = $max_sec;

my $missed_low_cov = $cutoff;
my $p = (1-$err)**$len;

while(1){
  my ($p_l, $p_le, $p_e, $p_ge, $p_g) = nkp_full($missed_low_cov,$cutoff,$p);


  if    ($p_ge > 1-1e-6){ $str_e6  = ">= $missed_low_cov to be missed is <1e-6"; last;}
  elsif ($p_ge > 0.9999){ $str_9999= ">= $missed_low_cov to be missed is <0.0001"; }
  elsif ($p_ge > 0.999) { $str_999 = ">= $missed_low_cov to be missed is <0.001"; }
  elsif ($p_ge > 0.99)  { $str_99  = ">= $missed_low_cov to be missed is <0.01"; }
  elsif ($p_ge > 0.9)   { $str_9   = ">= $missed_low_cov to be missed is <0.1"; }
  elsif ($p_ge > 0.5)   { $str_5   = ">= $missed_low_cov to be missed is <0.5"; }

  #print "$missed_low_cov\t$cutoff\t$p_ge\n";
  $missed_low_cov++;
  print STDERR ".";
}
  print STDERR "\n";

print <<EOD;

=================================================================================================
Read length is $len
Per-base error rate is $err
The probability for a read to be error-free is $p
The observed size of the largest cluster is $obv_pri
The estimated upper bound abundance for this tag is $max_array[0]
	(note! this abundance is different from the abundance of this species, which
	may has several copies of this tag with or without variations) 
The estimated size of the largest secondary cluster is $max_sec

Suggested cutoff value to remove small cluster is $cutoff
Representative sequences for large clusters will be used in next cluster run for OTU calculation
At this cutoff, the tags with abundance <$cutoff will be missed, 
In addition,

The probability for tags with abundance $str_5
The probability for tags with abundance $str_9
The probability for tags with abundance $str_99
The probability for tags with abundance $str_999
The probability for tags with abundance $str_9999
The probability for tags with abundance $str_e6

=================================================================================================

EOD

######
###### Given read length, per-base error rate, number of reads (coverage)
###### this function do N rendom tests and find
###### size of the primary cluster (the error-free cluster), standard deviation
###### size of the largest secondary cluster (the largest cluster with 1 error), standard deviation
######
sub sim1 {
  my ($len, $err, $cov1, $N) = @_;
  my ($ave_pri, $ave_sec, $sd_pri, $sd_sec);
  my ($i, $j, $k, $N1);
  
  #make pre-read array
  my $per_read_model_N = 0;
  my @per_read_model_10000 = ();
  
  for ($i=0; $i<=$len; $i++){
    my $p1 = nkp($len, $len-$i, (1-$err)); # when $i=0; no error
    $j = int($p1 * 10000);
    for ($k=0; $k<$j; $k++) {push(@per_read_model_10000,$i); $per_read_model_N++;}
    last if ($p1 < 1/$cov1/1000.0);
  } 

  my @cls=();
  my @pri=();
  my @sec=();
  my @secb=();
  my $total_cls=0;
  my $total_pri=0;
  my $total_sec=0;
  my $total_secb=0;

  for ($N1=0; $N1<$N; $N1++){
    my %err_count=();
    my %err_type = ();
    my $no_cls = 0;
    for ($i=0; $i<$cov1; $i++, $des++) {

      # for each read
      my $no_error_this_read = $per_read_model_10000[int(rand()*$per_read_model_N)];
      my $e_str = "0";
      for ($j=0; $j<$no_error_this_read; $j++){
        $k = int(rand()*$len); #error position
        my $c = $error_model_1000[int(rand()*$error_model_N)];
        $c = "$c.$k";
        $e_str .= $c;
      }
      if (defined($err_type{$e_str})) { $err_type{$e_str}++; }
      else                 { $no_cls++; $err_type{$e_str}=1; }
    }
    my $max_pri = $err_type{"0"};
    my $max_sec = 0;
    foreach $c (keys %err_type){
      $j = $err_type{$c};
      next if ($j == $max_pri);
      if ($j > $max_sec) {$max_sec=$j;}
    }
    my $max_secb = 0;
    foreach $c (keys %err_type){
      next unless ($c =~ /A|T|C/); #only wrong base call
      $j = $err_type{$c};
      next if ($j == $max_pri);
      if ($j > $max_secb) {$max_secb=$j;}
    }
    push(@cls, $no_cls);
    push(@pri, $max_pri);
    push(@sec, $max_sec);
    push(@secb, $max_secb);
    $total_pri += $max_pri;
    $total_sec += $max_sec;
    $total_secb += $max_secb;
    $total_cls += $no_cls;
  }

  my $ave_cls = int(10*$total_cls/$N)/10;
  my $ave_pri = int(10*$total_pri/$N)/10;
  my $ave_sec = int(10*$total_sec/$N)/10;
  my $ave_secb = int(10*$total_secb/$N)/10;

  my $sd_cls  = 0;
  my $sd_pri  = 0;
  my $sd_sec  = 0;
  my $sd_secb  = 0;
  for ($N1=0; $N1<$N; $N1++){
    $sd_cls += ($cls[$N1]-$ave_cls)*($cls[$N1]-$ave_cls);
    $sd_pri += ($pri[$N1]-$ave_pri)*($pri[$N1]-$ave_pri);
    $sd_sec += ($sec[$N1]-$ave_sec)*($sec[$N1]-$ave_sec);
    $sd_secb += ($secb[$N1]-$ave_secb)*($secb[$N1]-$ave_secb);
  }
  $sd_cls = int(10*sqrt($sd_cls/$N))/10;
  $sd_pri = int(10*sqrt($sd_pri/$N))/10;
  $sd_sec = int(10*sqrt($sd_sec/$N))/10;
  $sd_secb = int(10*sqrt($sd_secb/$N))/10;

  return ($ave_cls, $sd_cls, $ave_pri, $ave_sec, $ave_secb, $sd_pri, $sd_sec, $sd_secb);
}
############# sim1



# Binomial distribution
# introduction		http://en.wikipedia.org/wiki/Binomial_distribution
# online version	http://stattrek.com/Tables/Binomial.aspx
# 
sub nkp_full {
  my ($n, $k, $p) = @_;
  my $i;

  my ($p_l, $p_le, $p_e, $p_ge, $p_g);

  $p_l  = 0; #less
  $p_le = 0; #less or equal
  $p_ge = 0; #greater or equal
  $p_g  = 0; #greater

  for ($i=0; $i<=$n; $i++) {
    my $p1 = nkp($n, $i, $p);
    if ($i< $k) {$p_l  += $p1;}
    if ($i<=$k) {$p_le += $p1;}
    if ($i==$k) {$p_e   = $p1;}
    if ($i>=$k) {$p_ge += $p1;}
    if ($i> $k) {$p_g  += $p1;}
  }
  return ($p_l, $p_le, $p_e, $p_ge, $p_g);
}

sub nkp {
  my ($n, $k, $p) = @_;
  my $P0;

  if ($n == $k) {
    #(p**n)
    $P0 =  ($p**$n);
  }
  elsif ($k == 0) {
    $P0 = (1-$p)**$n;
  }
  else {
    #
    #  P0 = Cnk * (p**k) * ((1-p) ** (n-k))
    my @array1 = ();
    my @array2 = ();
    my @array3 = ();
    my $i;
    for ($i=1; $i<=$n;      $i++) { push(@array1, $i    ); }
    for ($i=1; $i<=$k;      $i++) { push(@array2, $i    ); }
    for ($i=1; $i<=($n-$k); $i++) { push(@array2, $i    ); }
    for ($i=1; $i<=$k;      $i++) { push(@array3, $p    ); }
    for ($i=1; $i<=($n-$k); $i++) { push(@array3, (1-$p)); }
    @array1 = sort {$b<=>$a} @array1;
    @array2 = sort {$b<=>$a} @array2;
    @array3 = sort {$a<=>$b} @array3;

    $P0=1;
    for ($i=0; $i<$n; $i++) {
      $P0 = $P0 * $array1[$i]*$array3[$i]/$array2[$i];
    }
  }
  return $P0;
}


