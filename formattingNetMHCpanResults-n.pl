#!/usr/bin/perl -w
use strict;
die "\nUsage: perl $0   <affinity_NetMHCpan.raw.txt>   <outdir>   <prefix>\n\n" unless (@ARGV == 3);

my $outdir = $ARGV[1];
my $prefix = $ARGV[2];
$ARGV[2] =~ s/.tumor.MaxQuant//g;

my $rawHeader = `grep Pos $ARGV[0] | head -1`;
chomp $rawHeader;
$rawHeader =~ s/^\s+|\s+$//g;
my @head = split /\s+/, $rawHeader;
my %header;
for(my $i=0; $i<@head; $i++) {
    $header{$head[$i]} = $i;
}

open (OUT, "> $outdir/$prefix.NetMHCpan_out.txt") || die $!;
print OUT join("\t", @head)."\n";

my %hash;
my %pepSEQ;
open (AF, "$ARGV[0]") || die $!;
while (<AF>) {
    chomp;
    next if ($_ !~ /^\s+/);
#next if ($_ !~ /HLA/ and $_ !~ /H-2/);
	next if ($_ !~ /HLA/);
	$_ =~ s/^\s+|\<\=//g;
    $_ .= "\tNB" if ($_ !~ /SB|WB/);
    my @inf = split /\s+/, $_;
    print OUT join("\t", @inf)."\n";
    
    my $peptide = $inf[$header{"Peptide"}];
    $pepSEQ{$peptide} = 1;
    push @{$hash{$peptide}{MHC}}, $inf[$header{"MHC"}];
    push @{$hash{$peptide}{Aff}}, $inf[$header{"Aff(nM)"}];
    push @{$hash{$peptide}{Bind}}, $inf[$header{"BindLevel"}];
    push @{$hash{$peptide}{Rank_EL}}, $inf[$header{"%Rank_EL"}];
    push @{$hash{$peptide}{Rank_BA}}, $inf[$header{"%Rank_BA"}];
}
close AF;
close OUT;

my $pepNum = scalar(keys %pepSEQ);
my %affHLAcount;
my ($affNum, $affNumRedundancy) = (0, 0);
open (FEATURES, "> $outdir/$prefix.peptides.affinity_features.txt") || die $!;
open (STRONGEST, "> $outdir/$prefix.mostBindingHLA.affinity.txt") || die $!;
print FEATURES join("\t", "#Peptide", "Length", "BinderNum", "MHC", "Aff.IC50nM", "Rank_EL(%)", "Rank_BA(%)", "BinderType")."\n";
print STRONGEST join("\t", qw(Peptide  MHC  Aff.IC50nM.min  Rank_EL.min(%)  Rank_BA.min(%)  BinderType BatchID))."\n";
foreach my $pep (sort{$a cmp $b} keys %hash) {
    my $mhc = join(",", @{$hash{$pep}{MHC}});
    my $aff = join(",",@{$hash{$pep}{Aff}});
    my $rankEL = join(",",@{$hash{$pep}{Rank_EL}});
    my $rankBA = join(",",@{$hash{$pep}{Rank_BA}});
    my $bindType = join(",", @{$hash{$pep}{Bind}});
    my $binderNum = 0;
    for (my $i=0; $i<@{$hash{$pep}{Bind}}; $i++) {
        $binderNum ++ if ($hash{$pep}{Bind}[$i] =~ /WB|SB/);
    }
    $affNum ++ if ($bindType =~ /SB|WB/);
    print FEATURES join("\t", $pep, length($pep), $binderNum, $mhc, $aff, $rankEL, $rankBA, $bindType)."\n";
    
    my $minMHC = $hash{$pep}{MHC}[0];
    my $minIC50 = $hash{$pep}{Aff}[0];
    my $minBindRankEL = $hash{$pep}{Rank_EL}[0];
    my $minBindRankBA = $hash{$pep}{Rank_BA}[0];
    my $minBindLevel = $hash{$pep}{Bind}[0];
    for (my $i=0; $i<@{$hash{$pep}{Rank_EL}}; $i++) {
        if ($minBindRankEL > $hash{$pep}{Rank_EL}[$i]) {
            $minMHC = $hash{$pep}{MHC}[$i];
            $minIC50 = $hash{$pep}{Aff}[$i];
            $minBindRankEL = $hash{$pep}{Rank_EL}[$i];
            $minBindRankBA = $hash{$pep}{Rank_BA}[$i];
            $minBindLevel = $hash{$pep}{Bind}[$i];
        }
    }
    $affHLAcount{$minMHC} ++ if ($minBindLevel =~ /SB|WB/);
    print STRONGEST join("\t", $pep, $minMHC, $minIC50, $minBindRankEL, $minBindRankBA, $minBindLevel,$ARGV[2])."\n";
    
}
close FEATURES;
close STRONGEST;

open (STAT, "> $outdir/$prefix.affinity.peptide.count.stat.txt") || die $!;
my @header = qw(BatchID	Item  affPEPcount  affPEPratio(%));
print STAT join("\t", @header)."\n";

my $total = 0;
foreach my $hla (sort{$a cmp $b} keys %affHLAcount) {
    print STAT join("\t",$ARGV[2], $hla, $affHLAcount{$hla}, sprintf("%.6f", $affHLAcount{$hla}/$pepNum*100))."\n";
    $total += $affHLAcount{$hla};
}
print STAT join("\t", $ARGV[2],"Non_Binder",$pepNum-$total,sprintf("%.6f", ($pepNum-$total)/$pepNum*100))."\n";
print STAT join("\t", $ARGV[2],"Total", $total, sprintf("%.6f", $total/$pepNum*100))."\n";
close STAT;
