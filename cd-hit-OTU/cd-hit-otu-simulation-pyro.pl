#!/usr/bin/perl

my $err     = shift;
my $len     = shift;
my $homo_no = shift;
my $obv_pri = shift;
my $sd1     = shift;
   $sd1     = 2 unless (defined($sd1));

#######################
#### usage ############
#######################
if (!$err or !$obv_pri) {
  print <<EOD;


Simulation to select the cluster size cutoff...

Usage: $0 err_file size-of-the-largest-cluster
	size-of-the-largest-cluster is the size of the largest cluster at 1st cluster run

Program developed by Weizhong Li's lab: http://weizhongli-lab.org/
EOD
  exit;
  ##### exit
}

my ($i, $j, $k);

my $p0 = (1-$err) ** $len;
my $p1 = nkp($len, 1, $err);
my $p2 = nkp($len, 2, $err);
my $p3 = nkp($len, 3, $err);
my $fr = 2/3; #### 1/3 are homo polymer based errors : 1/3 for insert, 1/3 for del
my $worst_base_span = $homo_no;

srand(124);

my ($cov1, $ave_pri, $sd_pri, $ave_sec, $sd_sec, $ave_3rd, $sd_3rd);
$cov1 = int($obv_pri / $p0 + 0.5);

  $ave_pri = $cov1 * $p0;
  $sd_pri  = sqrt($cov1*$p0*(1-$p0));

    ($ave_2nd, $sd_2nd) = simu_2nd($cov1*$p1*0.33, $len, 3);
    ($ave_3rd, $sd_3rd) = simu_3rd($cov1*$p2*0.33, $len, 3);

  if ($worst_base_span) { ####
    ($ave_2ndh, $sd_2ndh) = simu_2nd($cov1*$p1*$fr, $worst_base_span, 2);
    ($ave_3rdh, $sd_3rdh) = simu_3rd($cov1*$p2*$fr, $worst_base_span, 2);

    if ($ave_3rdh > $ave_3rd) {
       $ave_3rd = $ave_3rdh;
       $sd_3rd  = $sd_3rdh;
    }
  }


my $max_3rd = int($ave_3rd - $sd_3rd*$sd1);
my $cutoff = ($max_3rd > 3) ? $max_3rd : 3; 
my $missed_low_cov = int($cutoff / $p0 + 0.5);

print STDERR "\n";
print <<EOD;

=================================================================================================
The probability for a read to have 0 error  is $p0
The probability for a read to have 1 error  is $p1
The probability for a read to have 2 errors is $p2
The probability for a read to have 3 errors is $p3
Homoploymer sites: $homo_no
The observed size of the largest cluster is $obv_pri
The estimated upper bound abundance for this tag is $cov1
	(note! this abundance is different from the abundance of this species, which
	may has several copies of this tag with or without variations) 

Suggested cutoff value to remove small cluster is $cutoff
Representative sequences for large clusters will be used in next cluster run for OTU calculation
At this cutoff, the tags with abundance <$cutoff will be missed, 
In addition,

The probability for tags with abundance $missed_low_cov is about 0.5

=================================================================================================

EOD

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

sub scan_fastq {
  my $file = shift;
  my $no = 0;
  my $len = 0;
  my @score = ();
  my @count = ();

  my ($i, $j, $k, $ll);
  open(TMP, $file) || die "can not open $file";
  while($ll=<TMP>){
    if ($ll =~ /^\@/) {
      $ll = <TMP>; #### read sequence
      $ll =~ s/\s//g;
      my $len1 = length($ll);
      if ($len1 > $len) {$len = $len1;}

      $ll = <TMP>; die unless ($ll =~ /^\+/); #### read ID
      $ll = <TMP>; #### read quality score

      for ($i=0; $i<$len1; $i++) {
        my $c1 = ord(substr($ll,$i,1)) - 33;
        $score[$i] += $c1;
        $count[$i]++;
      }
      $no++;

      $i = int(rand()*10)*8; #skip some lines
      $j = 0;
      for ($j=0; $j<$i; $j++) {
        $ll=<TMP>;
        last unless ($ll);
      }
    }
  }
  close(TMP);

  my $p1 = 1;
  my $p2 = 1;
  for ($i=0; $i<$len; $i++) {
    $ave = $score[$i]/$count[$i];
    $e = 1 / (10 ** ($ave/10));
    $p = 1-$e;
    $p1 *= $p;

    $score[$i] = $e;
  }

  #p1 is the probability for a error-free read
  #p2 is the probability that only last base is wrong (last is the worst qualty base)
  $p2 = $p1 / (1-$score[-1]) * $score[-1];
  return ($len, $p1, $p2);
}
############# END scan_fastq

sub scan_err {
  my $file = shift;
  my ($i, $j, $k, $ll);
  my ($p0, $p1, $p2, $p3, $pw);

  my $PCR_f = 3;
  my @scores = ();
  $p0 = 1;
  my $score_sum=0;
  open(TMP, $file) || die "can not open $file";
  while($ll=<TMP>){
    chop($ll);
    my ($t1, $t2, $e) = split(/\t/,$ll);

    if (rand() < 1/$PCR_f) {
      $e += $PCR_e * $PCR_f;
    }

    push(@scores, $e);
    $score_sum += $e;
    $p0 *= (1-$e);
  }
  close(TMP);

  @scores = sort {$b <=> $a} @scores;
  my $len = $#scores+1;
  #now score in decrsing order

  my ($ei, $ej, $ek);
  $p1 = 0;
  for ($i=0; $i<$len; $i++) {
    $ei = $scores[$i];
    $p1 += $ei * $p0/(1-$ei);
  }

  $p2 = 0;
  for ($i=0; $i<$len; $i++) {
    $ei = $scores[$i];
    for ($j=$i+1; $j<$len; $j++) {
      $ej = $scores[$j];
      $p2 += $ei * $ej * $p0 / (1-$ei) / (1-$ej);
    }
  }

  $p3 = 0;
  for ($i=0; $i<$len; $i++) {
    $ei = $scores[$i];
    for ($j=$i+1; $j<$len; $j++) {
      $ej = $scores[$j];
      for ($k = $j+1; $k<$len; $k++) {
        $ek = $scores[$k];
        $p3 += $ei * $ej * $ek *  $p0 / (1-$ei) / (1-$ej) / (1-$ek);
      }
    }
  }

    my $lastn = 0;
    my $sum1=0;
    foreach $i (@scores) {
      $sum1 += $i;
      $lastn ++;
      last if ($sum1 > $score_sum*0.5);
      last if ($lastn > 20);
    }
    my $fraction_last10 = $sum1/$score_sum;
    return($p0, $p1, $p2, $p3, $fraction_last10, $lastn);
}
###### END scan_err

sub simu_2nd {
  my $N = shift;
  my $len = shift;
  my $m = shift; # type of errors, for wrong base call is 3
  my ($i00, $i, $j, $k, $r1, $r2, $r);

  my %types=();

  for ($i=0; $i<$N; $i++){
    $r1 = int(rand($len));
    $types{$r1}++;
  }

  my $sum=0;
  my $type_n=0;
  foreach $r (keys %types) {
    $type_n++;
    $sum += $types{$r};
  }
  my $ave = $sum/$type_n;

  $sum=0;
  foreach $r (keys %types) {
    $sum += ($types{$r}-$ave) *
            ($types{$r}-$ave);
  }
  $sd = sqrt($sum / $type_n);

  #ave and sd are only position based
  $ave = $ave / $m;
  $sd  = $sd  / $m;
  return($ave, $sd);
}
########## END simu_3rd


sub simu_3rd {
  my $N = shift;
  my $len = shift;
  my $m = shift; # type of errors, for wrong base call is 3
  my ($i00, $i, $j, $k, $r1, $r2, $r);

  my %types=();

  for ($i=0; $i<$N; $i++){
    $r1 = int(rand($len));
    $r2 = int(rand($len));
    if ($r1 > $r2) {
      $r=$r1; $r1=$r2; $r2=$r;
    }
    $types{"$r1.$r2"}++;
  }

  my $sum=0;
  my $type_n=0;
  foreach $r (keys %types) {
    $type_n++;
    $sum += $types{$r};
  }
  my $ave = $sum/$type_n;

  $sum=0;
  foreach $r (keys %types) {
    $sum += ($types{$r}-$ave) * 
            ($types{$r}-$ave);
  }
  $sd = sqrt($sum / $type_n);

  #ave and sd are only position based
  $ave = $ave / $m / $m;
  $sd  = $sd  / $m / $m;
  return($ave, $sd);
}
########## END simu_3rd












