#!/usr/bin/perl

my $clstr_ref=shift;
my $clstr_sbj=shift;
my $clstr_ref_no = 0;
my $clstr_sbj_no = 0;
my %id_2_clstr_ref=();
my %id_2_clstr_sbj=();
my %id_all = ();

my $clstr_no = -1;
open(TMP, $clstr_ref) || die "can not open $clstr_ref";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $clstr_no++;
  }
  else {
    chop($ll);
    if ($ll =~ /\*$/) { ; }
    else {
      if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
        $id_2_clstr_ref{$3} = $clstr_no;
        $id_all{$3}=1;
      }
    }
  }
}
close(TMP);
$clstr_ref_no = $clstr_no+1;


my $clstr_no = -1;
open(TMP, $clstr_sbj) || die "can not open $clstr_sbj";
while($ll=<TMP>){
  if ($ll =~ /^>/) {
    $clstr_no++;
  }
  else {
    chop($ll);
    if ($ll =~ /(\d+)(aa|nt), >(.+)\.\.\./) {
      $id_2_clstr_sbj{$3} = $clstr_no;
      $id_all{$3}=1;
    }
  }
}
close(TMP);
$clstr_sbj_no = $clstr_no+1;


my @sorted_sbj = 0..($clstr_sbj_no-1);
my %sbj_idx = ();

for ($i=0; $i<$clstr_sbj_no; $i++){
  my $idx = 0;
  my $n = 0;
  foreach $id (keys %id_2_clstr_sbj){
    next unless ($id_2_clstr_sbj{$id} == $i);
    $n++;
    $idx += $id_2_clstr_ref{$id};
  }
  if ($n>0) { $sbj_idx{$i} = $idx/$n; }  
  else {$sbj_idx{$i} = 0;}
}

@sorted_sbj = sort { $sbj_idx{$a} <=> $sbj_idx{$b} } @sorted_sbj;

my %mat=();
foreach $id (keys %id_all){
  if (defined($id_2_clstr_ref{$id})) {
    if (not defined($id_2_clstr_sbj{$id})) {
      $id_2_clstr_sbj{$id} = $clstr_sbj_no;
    }
  }
  elsif (defined($id_2_clstr_sbj{$id})) {
    if (not defined($id_2_clstr_ref{$id})) {
      $id_2_clstr_ref{$id} = $clstr_ref_no;
    }
  }
  $mat[$id_2_clstr_ref{$id}][$id_2_clstr_sbj{$id}]++;
}

push(@sorted_sbj, $clstr_sbj_no);

print "count";
for ($j=0; $j<=$clstr_sbj_no; $j++) {
  print "\tsbj", $j+1;
}
print "\n";


for ($i=0; $i<=$clstr_ref_no; $i++) {
  print "ref", $i+1;

  for ($j=0; $j<=$clstr_sbj_no; $j++) {
    $j1 = $sorted_sbj[$j];
    $k = defined($mat[$i][$j1]) ? $mat[$i][$j1] : 0;
    print "\t$k";
  }
  print "\n";
}
