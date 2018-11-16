#!/usr/bin/perl

my $delimiter = shift;
my $clstr = shift;
my $out1 = "$clstr.otu.txt";
my $out2 = "$clstr.sample.txt";
my $out3 = "$clstr.NAT.txt";

my $sample_no = 0;
my %sample_id = ();
my %count = ();
my %count_t = ();
my %count_s = ();
my %NAT_n0=();
my @NAT_array = qw/10 20 30 40 50 60 70 80 90 95 99 100/;
my %count_NAT = ();

my $OTU=0;
open(TMP, $clstr) || die "can not open $clstr";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $OTU++;
  }
  else {
    chop($ll);
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      my $rep = $2;
      my $sample_id = "count";
      if ($delimiter) {
        $sample_id = (split(/$delimiter/, $rep))[0];
        $sample_id{$sample_id}=1;
      }
      $count{$OTU}{$sample_id}++;
      $count_t{$OTU}++;
      $count_s{$sample_id}++;
    }
    else {
      die "format error $ll";
    }
  }
}
close(TMP);

my @sample_ids = sort keys %sample_id;



#calculate NAT50 and NAT80

foreach $sample_id (@sample_ids) {
  foreach $N (@NAT_array) {
    $N0 = $count_s{$sample_id} * $N / 100;

    my @otu_c = ();
    for ($i=1; $i<=$OTU; $i++){
      $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
      push(@otu_c, $k);
    }
    @otu_c = sort {$b <=> $a} @otu_c;

    $k = 0;
    for ($i=0; $i<$OTU; $i++) {
      $k+= $otu_c[$i];
      if ($k >= $N0){
        $count_NAT{$N}{$sample_id} = $i+1;
        last;
      }
    }
  }
}

#sort by NAT50
@sample_ids = sort {$count_NAT{'50'}{$a} <=> $count_NAT{'50'}{$b}} @sample_ids;


open(OUT1, "> $out1") || die "can not write $out1";
print OUT1 "OTU";
foreach $sample_id (@sample_ids){
  print OUT1 "\t$sample_id";
}
print OUT1 "\ttotal\n";
for ($i=1; $i<=$OTU; $i++){
  print OUT1 "OTU",$i;
  foreach $sample_id (@sample_ids){
    $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
    print OUT1 "\t$k";
  }
  print OUT1 "\t$count_t{$i}";
  print OUT1 "\n";
}
close(OUT1);


open(OUT2, "> $out2") || die "can not write $out2";
print OUT2 "Sample";
for ($i=1; $i<=$OTU; $i++){
  print OUT2 "\tOTU",$i;
}
print OUT2 "\ttotal\ttotal_otu\n";
foreach $sample_id (@sample_ids){
  print OUT2 $sample_id;
  for ($i=1; $i<=$OTU; $i++){
    $k = $count{$i}{$sample_id}? $count{$i}{$sample_id} : 0;
    print OUT2 "\t$k";
  }
  print OUT2 "\t$count_s{$sample_id}";
  print OUT2 "\t$count_NAT{100}{$sample_id}";
  print OUT2 "\n";
}
close(OUT2);


open(OUT3, "> $out3") || die "can not write $out3";
foreach $sample_id (@sample_ids){
  print OUT3 "\t$sample_id";
} 
print OUT3 "\n";
foreach $N (@NAT_array) {
  print OUT3 $N;
  foreach $sample_id (@sample_ids){
    $k = $count_NAT{$N}{$sample_id}? $count_NAT{$N}{$sample_id} : 0;
    print OUT3 "\t$k";
  } 
  print OUT3 "\n";
} 
close(OUT3);

