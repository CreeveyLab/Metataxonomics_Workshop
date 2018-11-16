#!/usr/bin/perl

my $delimiter = shift;
my $clstr = shift;
my $core_cutoff = shift;

my $out1  = "$clstr.smp";
my $out_core = "$clstr.core.txt";
my %count = ();
my $count_t = ();
my $sample_no = 0;

my @core_otus = ();
my %core_count = ();
my %core_total = ();
my %smp_ids = ();

my $OTU=0;
open(TMP, $clstr) || die "Can not open $clstr";
open(OUT1, "> $out1") || die "Can not write to $out1";
open(OUT2, "> $out_core") || die "Can not write to $out_core";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    if ($count_t) {
      print OUT1 ">OTU $OTU\t$sample_no samples\t$count_t reads\n";
      my @ids = keys %count;
         @ids = sort {$count{$b} <=> $count{$a}} @ids;
      $i=0;
      foreach $sample_id (@ids){
        $i++; 
        print OUT1 "$i\t$sample_id\t$count{$sample_id}\n";
      }

      if (($sample_no >= $core_cutoff) and ($core_cutoff)) {
        push(@core_otus, $OTU);
        $core_total{$OTU} = $count_t;
        foreach $sample_id (@ids){
          $core_count{$OTU}{$sample_id} = $count{$sample_id};
          $smp_ids{$sample_id} = 1;
        }
      }
    }
    $OTU++;
    %count=();
    $count_t=0;
    $sample_no=0;
  }
  else {
    chop($ll);
    if ($ll =~ /\d+(aa|nt), >(.+)\.\.\./) {
      my $rep = $2;
      if ($delimiter) {
        $sample_id = (split(/$delimiter/, $rep))[0];
      }
      if (not defined($count{$sample_id})) {
         $sample_no++;
         $count{$sample_id}=0;
      }
      $count{$sample_id}++;
      $count_t++;
    }
    else {
      die "format error $ll";
    }
  }
}

    if ($count_t) {
      print OUT1 ">OTU $OTU\t$sample_no samples\t$count_t reads\n";
      my @ids = keys %count;
         @ids = sort {$count{$b} <=> $count{$a}} @ids;
      $i=0;
      foreach $sample_id (@ids){
        $i++; 
        print OUT1 "$i\t$sample_id\t$count{$sample_id}\n";
      }

      if (($sample_no >= $core_cutoff) and ($core_cutoff)) {
        push(@core_otus, $OTU);
        $core_total{$OTU} = $count_t;
        foreach $sample_id (@ids){
          $core_count{$OTU}{$sample_id} = $count{$sample_id};
          $smp_ids{$sample_id} = 1;
        }
      }
    }
close(TMP);
close(OUT1);


if ($core_cutoff) {
  print OUT2 "OTU";
  my @smp_ids = sort keys %smp_ids;

  foreach $sample_id (@smp_ids) {
    print OUT2 "\t$sample_id";
  }
  print OUT2 "\ttotal\n";


  foreach $OTU (@core_otus) {
    print OUT2 "OTU",$OTU;
    foreach $sample_id (@smp_ids) {
      $k = $core_count{$OTU}{$sample_id}? $core_count{$OTU}{$sample_id} : 0;
      print OUT2 "\t$k";
    }
    print OUT2 "\t$core_total{$OTU}\n";
  }

}
close(OUT2);
