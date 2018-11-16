#!/usr/bin/perl

my $clstr = shift;
my $fasta = shift;
my $cutoff = shift;
my $pri_fasta_out = shift;
my $non_pri_fasta_out = shift;

my %pri_rep_count = ();
my %pri_rep_len = ();
my %non_pri = ();

my ($i, $j, $k);
open(TMP, $clstr) || die "can not open $clstr";
my $rep = "";
my $no = 0;
my $len = 0;
my @members = ();

while($ll=<TMP>){
  if ($ll =~ /^>/) {
    if (($no >= $cutoff) and ($len>0)) {
      $pri_rep_count{$rep} = $no;
      $pri_rep_len{$rep}   = $len;
    }
    elsif ($len>0){ #non pri
      push(@members, $rep);
      foreach $i (@members) {$non_pri{$i}=1;}
    }
    $rep = "";
    $no = 0;
    $len = 0;
    @members = ();
  }
  else {
    chop($ll);
    if ($ll =~ /\*$/) {
      $rep = "";
      if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
        $len = $1;
        $rep = $3;
      }
      else {
        die "format error $ll";
      }
    }
    else {
      if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
        push(@members, $3);
      }
    }
    $no++;
  }
}
    if (($no >= $cutoff) and ($len>0)) {
      $pri_rep_count{$rep} = $no;
      $pri_rep_len{$rep}   = $len;
    }
    elsif ($len>0){ #non pri
      push(@members, $rep);
      foreach $i (@members) {$non_pri{$i}=1;}
    }
close(TMP);


my @reps = keys %pri_rep_count;
   # sort by abundance
   @reps = sort {$pri_rep_count{$b} <=> $pri_rep_count{$a}} @reps;
my $rep_no = $#reps+1;

#make most abundance seq the longest
for ($i=0; $i<$rep_no-1; $i++) {
  $len1 = $pri_rep_len{$reps[$i]};

  for ($j=$i+1; $j<$rep_no; $j++){
    my $repj = $reps[$j];
    if ($pri_rep_len{$repj} > $len1) {
      $pri_rep_len{$repj} = $len1;
    }
  }
}
   
# retrive fasta for pri
my %pri_rep_des = ();
my %pri_rep_seq = ();
my $id;
my $min_len = $pri_rep_len{$reps[-1]};

my $seq = "";
my $des = "";

open(OOO, "> $non_pri_fasta_out") || die "can not write to $non_pri_fasta_out";
open(FASTA, $fasta) || die "can not open $fasta";
while($ll = <FASTA>) {
  if ($ll =~ /^>/) {
    if ($seq) {
      $id = substr($des,1);
      chop($id); $id =~ s/\s.+$//;

      if ($pri_rep_count{$id}) {
        $pri_rep_des{$id} = $des;
        $pri_rep_seq{$id} = $seq;
      }
      elsif ($non_pri{$id}) {
        print OOO $des;
        print OOO substr($seq, 0, $min_len), "\n"; #all seq in non-pri-fasta is shorter then pri fasta
      }
    }
    $des = $ll;
    $seq = "";
  }
  else {
    $ll =~ s/\s//g;
    $seq .= $ll;
  }
}
    if ($seq) {
      $id = substr($des,1);
      chop($id); $id =~ s/\s.+$//;

      if ($pri_rep_count{$id}) {
        $pri_rep_des{$id} = $des;
        $pri_rep_seq{$id} = $seq;
      }
      elsif ($non_pri{$id}) {
        print OOO $des;
        print OOO substr($seq, 0, $min_len); #all seq in non-pri-fasta is shorter then pri fasta
      }
    }

close(FASTA);
close(OOO);


open(OUT, "> $pri_fasta_out") || die "can not write to $pri_fasta_out";
foreach $id (@reps) {
  print OUT $pri_rep_des{$id};
  print OUT substr($pri_rep_seq{$id},0, $pri_rep_len{$id});
  print OUT "\n";
}
close(OUT);

