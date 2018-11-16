#!/usr/bin/perl

my $clstr = shift;
my $fasta = shift;
my $cutoff = shift;
my $pri_fasta_out = shift;
my $non_pri_fasta_out = shift;

my %pri_rep_count = ();
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

   
# retrive fasta for pri
my %pri_rep_des = ();
my %pri_rep_seq = ();
my $id;
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
        print OOO $seq, "\n";
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
        print OOO $seq, "\n";
      }
    }

close(FASTA);
close(OOO);


open(OUT, "> $pri_fasta_out") || die "can not write to $pri_fasta_out";
foreach $id (@reps) {
  print OUT $pri_rep_des{$id};
  print OUT $pri_rep_seq{$id};
  print OUT "\n";
}
close(OUT);

