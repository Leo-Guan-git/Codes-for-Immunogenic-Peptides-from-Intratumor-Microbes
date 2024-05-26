# @Author: guanxiangyu
# @Date:   2024-05-26 16:27:45
# @Last Modified by:   guanxiangyu
# @Last Modified time: 2024-05-26 16:40:47
#!/usr/bin/perl -w
use strict;
use POSIX;
use File::Basename;
use Data::Dumper;

@ARGV == 3 or die "Usage: perl $0 <Fastq.list> <Configs> <Outdir>\n";
my ($fqlist,$config,$outdir) = @ARGV;
system("mkdir -p $outdir") unless (-e "$outdir");
#######################################################################
my (%configs,%scripts);
my (%fqinputs,%cleanfqs,%splitfqs,%cleanfastq,%phred,%hlatype,%unmapfqs);
my (%bam,%realnbam,%realnbam2,%bqsrbam,%rnabam);
my (%snvvcf,%indelvcf,%fusionresult,%exprresult);

my @MegaBOLThead = qw{SampleName Read1 Read2};
my @chrs;
########################### Main Parameters ###########################

&readConfigs($config);

open LIST,"<$configs{reference_block}" or die "$!";
while(my $line = <LIST>){
	chomp $line;
	push @chrs, $line;
}
close LIST;

&readFqlist($fqlist);
&makeDir(\%fqinputs) if($configs{alignment} =~ /true/i);
&makeDir(\%bqsrbam) if($configs{alignment} =~ /false/i);
&runsplit(\%fqinputs) if ($configs{alignment} =~ /true/i and $configs{split_num} > 1);

if ($configs{fastqclean} =~ /true/i){
	if($configs{split_num} > 1){&FastqClean(\%splitfqs);}else{&FastqClean(\%fqinputs);}
}else{
	&FastqLink(\%fqinputs);
}

##allign BQSR HaplotypeCaller
&MegaBOLT(%cleanfqs) if($configs{alignment} =~ /true/i);

##QC
&megaboltqc(%bqsrbam) if($configs{megaboltqc} =~ /true/i);

# metaphlan
&metaphlan(%unmapfqs) if($configs{metaphlan} =~ /true/i);
# Kraken2
&kraken2(%unmapfqs) if($configs{kraken} =~ /true/i);

#&HLAtyping()
if(-e $configs{hlalist}){
	&readHLAlist();
}
else{
	&hlasomatic(%bqsrbam) if($configs{hlasomatic} =~ /true/i);
	&HLAminer(%cleanfqs) if($configs{hlaminer} =~ /true/i);
}

# Generate day info
chomp(my $day = `date +%Y%m%d`);
## run scripts
&edgeList($day);
## Generate run.sh
open RUNSH, ">$outdir/shell_run/run.$configs{projectname}.$day.sh" or die "$!";
print RUNSH "$configs{monitor} taskmonitor --q1 $configs{queue} --q2 $configs{queue2}";
print RUNSH " --P1 $configs{priority} --P2 $configs{priority2}";
print RUNSH " -p $configs{projectname} -i $outdir/shell_run/edge.$configs{projectname}.$day.list -f 1\n";
close RUNSH;

########################### Functions ###########################
sub readConfigs{
	my ($configfile) = @_;
	open CFG, "<$configfile" or die "$!";
	while(my $line = <CFG>){
		chomp $line;
		next if ($line =~ /^#/);
		$configs{$1} = $2 if ($line =~ /^(.*)=(.*);/);
	}
	close CFG;
}

sub readFqlist{
	my ($inputs) = @_;
	open INPUT, "<$inputs" or die "$!";
	my $num = 0;
	while(my $line = <INPUT>){
		chomp $line;
		my ($sampid,$dtype,$phred,$lane,$reads) = split /\s+/,$line;
		$dtype = "Tumor" if($dtype =~ /Tumor/i);
		$dtype = "Normal" if($dtype =~ /Normal/i);
		if ($configs{alignment} =~ /false/i){
			$bqsrbam{$sampid}{$dtype} = $reads; ## bam file input
		}else{
			my ($fq1,$fq2) = split /\,/, $reads;
			push @{$fqinputs{$sampid}{$dtype}{$lane}{fq1}}, $fq1;
			push @{$fqinputs{$sampid}{$dtype}{$lane}{fq2}}, $fq2;
		}
		$phred{$sampid}{$dtype} = $phred;
	}
	close INPUT;
}

sub makeDir{
	my ($hash) = @_;
	my %hash = %{$hash};
	system("mkdir -p $outdir/shell_run") unless (-e "$outdir/shell_run");
	for my $sampid(keys %hash){
		system("mkdir -p $outdir/$sampid/0.shell") unless (-e "$outdir/$sampid/0.shell");
		system("mkdir -p $outdir/$sampid/1.clean") unless (-e "$outdir/$sampid/1.clean");
		system("mkdir -p $outdir/$sampid/2.megaBOLT") unless (-e "$outdir/$sampid/2.megaBOLT");
		system("mkdir -p $outdir/$sampid/3.microbiome") unless (-e "mkdir -p $outdir/$sampid/3.microbiome");
		system("mkdir -p $outdir/$sampid/4.hlasomatic") unless (-e "mkdir -p $outdir/$sampid/4.hlasomatic");
		## shell
		system("mkdir -p $outdir/$sampid/0.shell/a01.clean") unless (-e "$outdir/$sampid/0.shell/a01.clean");
		system("mkdir -p $outdir/$sampid/0.shell/a02.megaBOLT") unless (-e "$outdir/$sampid/0.shell/a02.megaBOLT");
		system("mkdir -p $outdir/$sampid/0.shell/a03.metaphlan") unless (-e "$outdir/$sampid/0.shell/a03.metaphlan");
		system("mkdir -p $outdir/$sampid/0.shell/a04.kraken2")unless(-e "$outdir/$sampid/0.shell/a04.kraken2");
		system("mkdir -p $outdir/$sampid/0.shell/a05.hla")unless(-e "$outdir/$sampid/0.shell/a05.hla");
		##qc
		system("mkdir -p $outdir/$sampid/2.megaBOLT/QC") unless (-e "$outdir/$sampid/2.megaBOLT/QC");
				
		if($configs{alignment} =~ /true/i){
			for my $dtype(keys %{$hash{$sampid}}){
				for my $lane(keys %{$hash{$sampid}{$dtype}}){
					system("mkdir -p $outdir/$sampid/1.clean/$dtype/$lane") unless (-e "$outdir/$sampid/1.clean/$dtype/$lane");
				}
			}
		}
	}
}

sub generateShell{
	my ($output_shell, $content, $finish_string) = @_;
	my $shell_name = basename($output_shell);
	$finish_string ||= "Still_waters_run_deep";
	open OUT, ">$output_shell" or die "Cannot open file $output_shell:$!";
	print OUT "#!/bin/bash\n";
	print OUT "echo hostname: `hostname`\n";
	print OUT "echo ==========start at : `date` ==========\n";
	print OUT "$content && \\\n" if($content);
	print OUT "echo ==========end at : `date` ========== && \\\n";
	print OUT "echo $finish_string 1>&2 && \\\n";
	print OUT "echo $finish_string > $output_shell.sign\n";
	print OUT "qstat -j $shell_name > $output_shell.log\n";
	close OUT;
}

sub runsplit{
	my ($hash) = @_;
	my %hash = %{$hash};
	for my $sampid(sort keys %hash){
		for my $dtype(sort keys %{$hash{$sampid}}){
			for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
				my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
				my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
				for my $i(0..$#fq1){
					my $tmp=substr($fq1[$i],0,-3);my $readsnum;
					if(-e "$tmp.fqStat.txt"){
						$readsnum=`cat $tmp.fqStat.txt | grep 'ReadNum' | cut -f2`;chomp($readsnum);
						$splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} = "true";
					}else{
						push @{$splitfqs{$sampid}{$dtype}{$lane}{fq1}}, $fq1[$i];
						push @{$splitfqs{$sampid}{$dtype}{$lane}{fq2}}, $fq2[$i];
						$splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} = "false";
						next;
					}
					my $splitrow=POSIX::ceil ($readsnum/$configs{split_num});
					$splitrow *=4;
					my ($unzipfile, $unzipfile1);
					if(($splitrow*($configs{split_num}-1)/4)<$readsnum){
						my $shell="$outdir/$sampid/0.shell/a01.clean/a01.fqsplit_$sampid\_$dtype\_$lane\_$i.1.sh";
						my $cmd;
						if($fq1[$i]=~/\.gz/){
							$unzipfile=basename $fq1[$i];
							$unzipfile=substr($unzipfile,0,-3);
							$cmd = "gunzip -c $fq1[$i] > $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile && \\\n";
							$cmd.= "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile -d -a 2";
							$cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_ && \\\n";
							$cmd.= "rm -rf $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile";
						}else{
							$unzipfile=$fq1[$i];
							$cmd = "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile -d -a 2";
							$cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_";
						}
						generateShell($shell,$cmd);
						$scripts{split}{$sampid}{$dtype}{$lane}{$i.1}=$shell;
						$shell="$outdir/$sampid/0.shell/a01.clean/a01.fqsplit_$sampid\_$dtype\_$lane\_$i.2.sh";
						undef $cmd;
						if($fq2[$i]){
							if($fq2[$i]=~/\.gz/){
								$unzipfile1=basename $fq2[$i];
								$unzipfile1=substr($unzipfile1,0,-3);
								$cmd = "gunzip -c $fq2[$i] > $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 && \\\n";
								$cmd.= "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 -d -a 2";
								$cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_ && \\\n";
								$cmd.= "rm -rf $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1";
							}else{
								$unzipfile1=$fq2[$i];
								$cmd = "split -$splitrow $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1 -d -a 2";
								$cmd.= " $outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_";
							}
						}
						generateShell($shell,$cmd);
						$scripts{split}{$sampid}{$dtype}{$lane}{$i.2}=$shell;
					}else{
						print "split err: the split num is too high for $fq1[$i]";
						print " and $fq2[$i]" if($fq2[$1]);
						print "\n";
					}
					my $s = $configs{split_num} -1;
					for my $n(0..$s){
					my ($splitfq1,$splitfq2);
					if($n < 10){
						$splitfq1 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_0$n";
						$splitfq2 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_0$n" if($unzipfile1);
					}else{
						$splitfq1 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile\_$n";
						$splitfq2 = "$outdir/$sampid/1.clean/$dtype/$lane/$unzipfile1\_$n" if($unzipfile1);
					}
						push @{$splitfqs{$sampid}{$dtype}{$lane}{fq1}}, $splitfq1;
						push @{$splitfqs{$sampid}{$dtype}{$lane}{fq2}}, $splitfq2;
					}
				}
			}
		}
	}
}

sub FastqClean{
	my ($hash) = @_;
	my %hash = %{$hash};
	for my $sampid(sort keys %hash){
		for my $dtype(sort keys %{$hash{$sampid}}){
			if ($dtype =~ /Tumor/i or $dtype =~ /Normal/i){
				$cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$dtype.csv";
				open ELIST, ">$outdir/$sampid/2.megaBOLT/$dtype.csv" or die "$!";
			}
			for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
				my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
				my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
				for my $i(0..$#fq1){
					my ($out_fq1, $out_fq2);
					$out_fq1 = "$outdir/$sampid/1.clean/$dtype/$lane/fq$i/clean_1.$i.fq.gz";
					$out_fq2 = "$outdir/$sampid/1.clean/$dtype/$lane/fq$i/clean_2.$i.fq.gz";
					system("mkdir -p $outdir/$sampid/1.clean/$dtype/$lane/fq$i") unless (-e "$outdir/$sampid/1.clean/$dtype/$lane/fq$i");
					my $shell = "$outdir/$sampid/0.shell/a01.clean/a01.cleanfq_$sampid\_$dtype\_$lane\_$i.sh";
					my $cmd;
					if($configs{soapnukeclean} =~ /false/){
						$cmd = "$configs{fastp} -i $fq1[$i] -o $out_fq1";
						$cmd.= " -I $fq2[$i] -O $out_fq2" if($fq2[$i]);
						$cmd.= " -j $outdir/$sampid/1.clean/$dtype/$lane/fq$i/$sampid.$dtype.$lane.fastp.json";
						$cmd.= " -h $outdir/$sampid/1.clean/$dtype/$lane/fq$i/$sampid.$dtype.$lane.fastp.html";
						$cmd.= " $configs{fastpparameter}";
					}else{
						$cmd = "$configs{soapnuketool} filter -1 $fq1[$i] -C clean_1.$i.fq.gz";
						$cmd.= " -2 $fq2[$i] -D clean_2.$i.fq.gz" if($fq2[$i]);
						$cmd.= " $configs{soapnukeparameter} -o $outdir/$sampid/1.clean/$dtype/$lane/fq$i";
					}
					if($configs{split_num} > 1 && $hash{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
						$cmd.= "&& \\\nrm -rf $fq1[$i]";
						$cmd.= " && \\\nrm -rf $fq2[$i]" if($fq2[$i]);

					}
					generateShell($shell,$cmd);
					$scripts{clean}{$sampid}{$dtype}{$lane}{$i}=$shell;
					push @{$cleanfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
					push @{$cleanfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
				}
			}
			my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
			my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
			print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
			print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
			print ELIST "\n";
			close ELIST;
		}
	}
}

sub FastqLink{
	my ($hash) = @_;
	my %hash = %{$hash};
	for my $sampid(sort keys %hash){
		for my $dtype(sort keys %{$hash{$sampid}}){
			if ($dtype =~ /Tumor/i or $dtype =~ /Normal/i){
				$cleanfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$dtype.csv";
				open ELIST, ">$outdir/$sampid/2.megaBOLT/$dtype.csv" or die "$!";
			}
			for my $lane(sort keys %{$hash{$sampid}{$dtype}}){
				my @fq1 = @{$hash{$sampid}{$dtype}{$lane}{fq1}};
				my @fq2 = @{$hash{$sampid}{$dtype}{$lane}{fq2}};
				for my $i(0..$#fq1){
					my ($out_fq1, $out_fq2);
					$out_fq1 = "$outdir/$sampid/1.clean/$dtype/$lane/clean_1.$i.fq.gz";
					$out_fq2 = "$outdir/$sampid/1.clean/$dtype/$lane/clean_2.$i.fq.gz" if($fq2[$i]);
					`ln -s $fq1[$i] $out_fq1`;
					`ln -s $fq2[$i] $out_fq2` if($fq2[$i]);
					push @{$cleanfastq{$sampid}{$dtype}{fq1}}, $out_fq1;
					push @{$cleanfastq{$sampid}{$dtype}{fq2}}, $out_fq2;
				}
			}
			my $mega_fq1 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq1}});
			my $mega_fq2 = join(",", @{$cleanfastq{$sampid}{$dtype}{fq2}});
			print ELIST join("\t","$sampid\-$dtype",$mega_fq1);
			print ELIST "\t$mega_fq2" unless($mega_fq2 =~ /,{2,}/);
			print ELIST "\n";
			close ELIST;
		}
	}
}

sub MegaBOLT{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $Suffix = "mm2";
			$Suffix = "bwa" if ($configs{BWA} =~ /true/i);
			$bqsrbam{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.bam";
			$unmapfqs{$sampid}{$dtype} = "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/unmap.fq.gz";
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/a02.megaBOLT_$dtype.sh";
			# my $cmd = "MegaBOLT --type basic --runtype WGS --outputdir $outdir/$sampid/2.megaBOLT/ --bwa 1";
			my $cmd = "MegaBOLT --type basic --runtype WGS --outputdir $outdir/$sampid/2.megaBOLT/";
			$cmd.= " --bwa 1" if ($configs{BWA} =~ /true/i); 
			# $cmd.= " --list $outdir/$sampid/2.megaBOLT/Tumor/Tumor.csv --list2 $outdir/$sampid/2.megaBOLT/Normal/Normal.csv";
			$cmd.= " --list $outdir/$sampid/2.megaBOLT/$dtype.csv";
			$cmd.= " --ref $configs{reference}" if($configs{reference});
			$cmd.= " --knownSites $configs{Cosmic}" if($configs{Cosmic});
			$cmd.= " --knownSites $configs{GATKdbsnp}" if($configs{GATKdbsnp});
			$cmd.= " --knownSites $configs{'1kg_phase1_snp'}" if($configs{'1kg_phase1_snp'});
			$cmd.= " --knownSites $configs{Mills_indel}" if($configs{Mills_indel});
			$cmd.= " --vcf $configs{GATKdbsnp}" if($configs{GATKdbsnp});
			if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/i){
				$cmd.= " && \\\n";
				my $mydir = "$outdir/$sampid/5.somatic/tumor_only";
				system("mkdir -p $mydir") unless(-e $mydir);
				$cmd.= "MegaBOLT --type mutect2 --mutect2-input $bqsrbam{$sampid}{$dtype} --tumor $sampid-Tumor";
				$cmd.= " --ref $configs{reference}" if($configs{reference});
				$cmd.= " --vcf $configs{GATKdbsnp}" if($configs{GATKdbsnp});
				$cmd.= " --knownSites $configs{Cosmic}" if($configs{Cosmic});
				$cmd.= " --knownSites $configs{GATKdbsnp}" if($configs{GATKdbsnp});
				$cmd.= " --knownSites $configs{'1kg_phase1_snp'}" if($configs{'1kg_phase1_snp'});
				$cmd.= " --knownSites $configs{Mills_indel}" if($configs{Mills_indel});		
				$cmd.= " --outputdir $mydir && \\\n";
				$cmd.= "cp $mydir/output/output.mutect2.vcf.gz $mydir/$sampid.somatic.vcf.gz && \\\n";
				$cmd.= "cp $mydir/output/output.mutect2.vcf.gz.tbi $mydir/$sampid.somatic.vcf.gz.tbi && \\\n";
				$cmd.= "rm -rf $mydir/output";
			}
			if($configs{metaphlan} =~ /true/i || $configs{kraken} =~ /true/i || $configs{hlaminer} =~ /true/i){
				$cmd.= " && \\\n";
				$cmd.= "$configs{samtools} fastq -@ 10 -N -f 4 -o $unmapfqs{$sampid}{$dtype} $bqsrbam{$sampid}{$dtype}";
			}
			if($configs{hlaminer} =~ /true/i){
				$cmd.= " && \\\n";
				$cmd.= "$configs{samtools} view -h -@ 10 -L $configs{TASRbed} $bqsrbam{$sampid}{$dtype} |";
				$cmd.= " $configs{samtools} fastq -@ 10 -N -o $outdir/$sampid/2.megaBOLT/$sampid-$dtype/HLA.fq.gz";
			}
			generateShell($shell, $cmd);
			$scripts{MegaBOLT}{$sampid}{$dtype} = $shell;

			# germline mutation annotation
			$shell = "$outdir/$sampid/0.shell/a02.megaBOLT/a02.annosnp_$dtype.sh";
			$cmd = "source /hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/RNApipeline/.bashrc_guanxiangyu\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz.tbi\n";
			$cmd.= "gzip -d $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf.gz\n";
			$cmd.= "$configs{GATK4} VariantFiltration -R $configs{reference} $configs{SNPfilterParameter}";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf";
			$cmd.= " --verbosity ERROR -O $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf && \\\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.$Suffix.sortdup.bqsr.hc.vcf && \\\n";
			$cmd.= "awk '\$1~/^#/ || \$7==\"PASS\"' $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf >";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf && \\\n";
			$cmd.= "$configs{bgzip} -c $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf >";
			$cmd.= " $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz && \\\n";
			$cmd.= "$configs{tabix} -p vcf $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz && \\\n";
			$cmd.= "$configs{bgzip} -@ 2 -f $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf && \\\n";
			$cmd.= "$configs{tabix} -p vcf $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.filter.vcf.gz && \\\n";
			$cmd.= "rm $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf && \\\n";
			$cmd.= "perl $configs{ANNOVAR} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid.$dtype.hc.final.vcf.gz";
			$cmd.= " $configs{ANNOVAR_refdir} -out $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot $configs{annovar_par} && \\\n";
			$cmd.= "#less $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot.homo_multianno.txt | cut -f1,2,4,5,6 |";
			$cmd.= " grep 'rs' > $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.txt && \\\n";
			$cmd.= "ls $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.germline.variants.annot* | while read file;";
			$cmd.= "do $configs{bgzip} -@ 2 -f \$file; done";
			generateShell($shell,$cmd);
			$scripts{gatksnp_ann}{$sampid}{$dtype}=$shell;

			##snp chemicalDrug
			$shell = "$outdir/$sampid/0.shell/a16.medicine/a16.chemicalDrug_$dtype.sh";
			$cmd = "perl $configs{chemicaldrug_tool} $outdir/$sampid/2.megaBOLT/$sampid-$dtype/$sampid-$dtype.final.variants.txt >";
			$cmd.= " $outdir/$sampid/9.medicine/chemicalDrug/chemicalDrug.out";
			generateShell($shell,$cmd);
			$scripts{chemicalDurg_ann}{$sampid}{$dtype}=$shell;
			# }

			if($configs{tumor_only} =~ /true/i && $dtype =~ /Tumor/i){
				# tumor only result filter
				my $mydir = "$outdir/$sampid/5.somatic/tumor_only";
				system("mkdir -p $mydir/contamination/tmp") unless (-e "$mydir/contamination/tmp");
				system("mkdir -p $mydir/artifact/tmp") unless (-e "$mydir/artifact/tmp");
				my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/tumor_only.bqsr.contamination.sh";
				my $tbam = $bqsrbam{$sampid}{$dtype};
				my $cmd = "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
				$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
				$cmd.= " GetPileupSummaries -I $tbam -V $configs{GetPileupSummaries}";
				$cmd.= " -L $configs{reference_bed}";
				$cmd.= " -O $mydir/contamination/$sampid-Tumor.pileups.table && \\\n";
				$cmd.= "$configs{GATK4} --java-options \"-Xmx6g -XX:ParallelGCThreads=1";
				$cmd.= " -Dsamjdk.compression_level=5 -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/contamination/tmp\"";
				$cmd.= " CalculateContamination -I $mydir/contamination/$sampid-Tumor.pileups.table";
				$cmd.= " -O $mydir/contamination/$sampid.contamination.table";
				generateShell($shell,$cmd);
				$scripts{tumor_only}{$sampid}{contamination} = $shell;
				
				## Filter somatic mutation Calls
				$shell = "$outdir/$sampid/0.shell/a02.megaBOLT/tumor_only.bqsr.filter.sh";
				$cmd = "mkdir -p $mydir/tmp && \\\n";
				$cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
				$cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterMutectCalls";
				$cmd.= " -V $mydir/$sampid.somatic.vcf.gz -contamination-table $mydir/contamination/$sampid.contamination.table";
				$cmd.= " -O $mydir/$sampid.somatic.filter1.vcf.gz && \\\n";
				# $cmd.= "$configs{GATK4} --java-options \"-Xmx25g -XX:ParallelGCThreads=1 -Dsamjdk.compression_level=5";
				# $cmd.= " -XX:-UseGCOverheadLimit -Djava.io.tmpdir=$mydir/tmp\" FilterByOrientationBias";
				# $cmd.= " --artifact-modes \"G/T\" --artifact-modes \"C/T\" -V $mydir/$sampid.somatic.filter1.vcf.gz";
				# $cmd.= " -P $mydir/artifact/$sampid.artifact.pre_adapter_detail_metrics.txt";
				# $cmd.= " -O $mydir/$sampid.somatic.filter2.vcf.gz && \\\n";
				# $cmd.= "gzip -dc $mydir/$sampid.somatic.filter2.vcf.gz |";
				$cmd.= "gzip -dc $mydir/$sampid.somatic.filter1.vcf.gz |";
				$cmd.= " perl -ne \'if(/^#/){print \$_;}else{my \$FILTER=(split /\\t/)[6];if(\$FILTER=~\"PASS\"){print \$_;}}\' |";
				$cmd.= " gzip -c > $mydir/$sampid.somatic.PASS.vcf.gz && \\\n";
				$cmd.= "perl $configs{ANNOVAR} $mydir/$sampid.somatic.PASS.vcf.gz";
				$cmd.= " $configs{ANNOVAR_refdir} -out $mydir/$sampid.somatic.annot $configs{annovar_par} && \\\n";
				$cmd.= "ls $mydir/$sampid.somatic.annot* | while read file;";
				$cmd.= "do $configs{bgzip} -@ 2 -f \$file; done && \\\n";
				$cmd.= "rm -rf $mydir/artifact $mydir/contamination $mydir/$sampid.somatic.vcf.gz* $mydir/$sampid.somatic.filter1* $mydir/$sampid.somatic.filter2*";
				generateShell($shell,$cmd);
				$scripts{tumor_only}{$sampid}{filter} = $shell;
			}
		}
	}
}

sub megaboltqc{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my ($tbam, $nbam);
		for my $dtype(keys %{$hash{$sampid}}){
			$tbam = $hash{$sampid}{$dtype} if($dtype =~ /Tumor/i);
			$nbam = $hash{$sampid}{$dtype} if($dtype =~ /Normal/i);
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/$dtype.coverage.sh";
			my $mydir = "$outdir/$sampid/2.megaBOLT/coverage";
			system("mkdir -p $mydir/$dtype") unless (-e "$mydir/$dtype");
			my $cmd = "$configs{samtools} bedcov $configs{targetRegion} $hash{$sampid}{$dtype} > $mydir/$dtype/bedcov.txt && \\\n";
			$cmd.= "$configs{samtools} coverage $hash{$sampid}{$dtype} > $mydir/$dtype/coverage.txt && \\\n";
			$cmd.= "$configs{samtools} depth -a -o $mydir/$dtype/depth.txt $hash{$sampid}{$dtype} && \\\n";
			$cmd.= "$configs{bgzip} -@ 2 $mydir/$dtype/depth.txt";
			generateShell($shell,$cmd);
			$scripts{coverage}{$sampid}{$dtype} = $shell;          
		}
		unless($configs{tumor_only} =~ /true/i || ! -e $configs{qcvcf}){
			my $mydir = "$outdir/$sampid/2.megaBOLT/QC";
			system("mkdir -p $mydir/tmp") unless (-e "$mydir/tmp");
			my $shell = "$outdir/$sampid/0.shell/a02.megaBOLT/QC.sh";
			my $cmd = "$configs{python} $configs{qctool} --bam1 $nbam --bam2 $tbam";
			$cmd.= " --output $mydir/bam_matcher.report.txt --vcf $configs{qcvcf}";
			$cmd.= " --reference $configs{reference} --scratch-dir $mydir/tmp $configs{qcpar} && \\\n";
			$cmd.= "rm -rf $mydir/tmp";
			generateShell($shell,$cmd);
			$scripts{megaboltqc}{$sampid} = $shell;
		}
	}
}

sub metaphlan{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $shell = "$outdir/$sampid/0.shell/a03.metaphlan/a03.metaphlan_$dtype.sh";
			my $fastq = $hash{$sampid}{$dtype};
			my $cmd = "export PATH=\"/hwfssz1/ST_SUPERCELLS/P21Z10200N0125/guanxiangyu/02.apps/miniconda3/envs/metaphlan/bin:\$PATH\" && \\\n";
			$cmd.= "mkdir -p $outdir/$sampid/3.microbiome/$dtype && \\\n";
			$cmd.= "$configs{metaphlan3} $fastq --nproc $configs{metaphlan_nproc} --input_type fastq";
			$cmd.= " --bowtie2out $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_metagenome.bowtie2.bz2";
			$cmd.= " -o $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_metaphlan_wavg_l.txt";
			$cmd.= " --index $configs{metaphlan_index} --bowtie2db $configs{bowtie2db}";
			$cmd.= " --bowtie2_exe $configs{bowtie2} --bowtie2_build $configs{bowtie2_build}";
			$cmd.= " --sample_id ${sampid}_${dtype} --add_viruses --stat wavg_l && \\\n";
			$cmd.= "$configs{metaphlan3} $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_metagenome.bowtie2.bz2 --input_type bowtie2out";
			$cmd.= " -o $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_metaphlan_wavg_g.txt";
			$cmd.= " --index $configs{metaphlan_index} --bowtie2db $configs{bowtie2db}";
			$cmd.= " --sample_id ${sampid}_${dtype} --add_viruses --stat wavg_g";
			generateShell($shell,$cmd);
			$scripts{metaphlan}{$sampid}{$dtype} = $shell;
		}
	}
}

sub kraken2{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			my $shell = "$outdir/$sampid/0.shell/a04.kraken2/a04.kraken2_$dtype.sh";
			my $fastq = $hash{$sampid}{$dtype};
			my $cmd = "$configs{kraken2} --db $configs{Kraken2DB} --threads $configs{Kraken2_threads} $configs{Kraken2par}";
			$cmd.= " --output $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_kraken2.out";
			$cmd.= " --report $outdir/$sampid/3.microbiome/$dtype/${sampid}_${dtype}_kraken2.report $fastq";
			generateShell($shell, $cmd);
			$scripts{kraken2}{$sampid}{$dtype} = $shell;
		}
	}
}
sub readHLAlist{
	open LIST, "<$configs{hlalist}" or die "$!";
	while(my $line = <LIST>){
		chomp $line;
		my ($sampid,$hlafile) = (split /\t/,$line)[0,1];
		$hlatype{$sampid} = $hlafile;
	}
	close LIST;
}

sub hlasomatic{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		my $mydir = "$outdir/$sampid/5.somatic/hlasomatic";
		for my $dtype(keys %{$hash{$sampid}}){
			if($dtype =~ /RNA/i){ next; }
			system("mkdir -p $mydir/$dtype/Polysolver") unless(-e "$mydir/$dtype/Polysolver");
			my $bam = $hash{$sampid}{$dtype};
			my $shell = "$outdir/$sampid/0.shell/a18.hlasomatic/$sampid.$dtype.Polysolver.sh";
			my $cmd = "source $configs{polysolver_config} && \\\n";
			$cmd.= "$configs{polysolver} $bam $configs{hlapar} $mydir/$dtype/Polysolver && \\\n";
			$cmd.= "perl $configs{combinescript} $mydir/$dtype/Polysolver $mydir/$sampid.$dtype.Polysolver.csv && \\\n";
			$cmd.= "rm -rf $mydir/$dtype/Polysolver/*.bam $mydir/$dtype/Polysolver/*.lik* $mydir/$dtype/Polysolver/*.sam";
			generateShell($shell, $cmd);
			$scripts{polysolver}{$sampid}{$dtype} = $shell;
			$hlatype{$sampid}{$dtype} = "$mydir/$sampid.$dtype.Polysolver.csv";
		}
	}
}

sub HLAminer{
	my (%hash) = @_;
	for my $sampid(keys %hash){
		for my $dtype(keys %{$hash{$sampid}}){
			if ($dtype =~ /RNA/i){ next; }
			my $mydir = "$outdir/$sampid/5.somatic/hlasomatic/$dtype/HLAminer";
			system("mkdir -p $mydir") unless (-e $mydir);
			open FOF, "> $mydir/$sampid.$dtype.fof";
			# print FOF "$outdir/$sampid/2.sentieon/$dtype/$sampid.$dtype.HLA.fq.gz\n";
			print FOF "$outdir/$sampid/2.megaBOLT/$sampid-$dtype/HLA.fq.gz\n";
			print FOF "$unmapfqs{$sampid}{$dtype}\n";
			close FOF;
			my $shell = "$outdir/$sampid/0.shell/a18.hlasomatic/$sampid.$dtype.HLAminer.sh";
			my $cmd = "###TASR\necho \"Running TASR...\"\n";
			$cmd.= "$configs{HLAminer}/TASR -f $mydir/$sampid.$dtype.fof -s $configs{TASRtargets}";
			$cmd.= " -b $mydir/$sampid.$dtype.TASRhla $configs{TASRpar} && \\\n";
			$cmd.= "###Restrict 200nt+ contigs\n";
			$cmd.= "cat $mydir/$sampid.$dtype.TASRhla.contigs |";
			$cmd.= " perl -ne 'if(/size(\\d+)/){if(\$1>=200){\$flag=1;print;}else{\$flag=0;}}else{print if(\$flag);}' >";
			$cmd.= " $mydir/$sampid.$dtype.TASRhla200.contigs && \\\n";
			$cmd.= "###Create a [NCBI] blastable database\n";
			$cmd.= "echo \"Formatting blastable database...\"\n";
			$cmd.= "$configs{HLAminer}/formatdb -p F -i $mydir/$sampid.$dtype.TASRhla200.contigs -l $mydir/$sampid.$dtype.formatdb.log && \\\n";
			$cmd.= "###Align contigs against database\n";
			$cmd.= "echo \"Aligning TASR contigs to HLA references...\"\n";
			$cmd.= "$configs{HLAminer}/parseXMLblast.pl -c $configs{HLAminer}/ncbiBlastConfig.txt -d $configs{TASRtargets}";
			$cmd.= " -i $mydir/$sampid.$dtype.TASRhla200.contigs -o 0 -a 1 > $mydir/$sampid.$dtype.tig_vs_hla-ncbi.coord && \\\n";
			$cmd.= "###Align HLA references to contigs\n";
			$cmd.= "echo \"Aligning HLA references to TASR contigs (go have a coffee, it may take a while)...\"\n";
			$cmd.= "$configs{HLAminer}/parseXMLblast.pl -c $configs{HLAminer}/ncbiBlastConfig.txt -d $mydir/$sampid.$dtype.TASRhla200.contigs";
			$cmd.= " -i $configs{TASRtargets} -o 0 -a 1 > $mydir/$sampid.$dtype.hla_vs_tig-ncbi.coord && \\\n";
			$cmd.= "###Predict HLA alleles\n";
			$cmd.= "echo \"Predicting HLA alleles...\"\n";
			$cmd.= "$configs{HLAminer}/HLAminer.pl -b $mydir/$sampid.$dtype.tig_vs_hla-ncbi.coord -r $mydir/$sampid.$dtype.hla_vs_tig-ncbi.coord";
			$cmd.= " -c $mydir/$sampid.$dtype.TASRhla200.contigs -h $configs{TASRtargets} -m $mydir -p $configs{hlanomp} && \\\n";
			$cmd.= "$configs{HLAminer_Parse} $mydir/HLAminer_HPTASR.csv > $outdir/$sampid/5.somatic/hlasomatic/$sampid.$dtype.HLAminer_HPTASR.csv";
			generateShell($shell, $cmd);
			$scripts{HLAminer}{$sampid}{$dtype} = $shell;
		}
	}
}


sub edgeList{
	my ($day) = @_;
	my $mydir = "$outdir/shell_run";
	open EDGE, ">$outdir/shell_run/edge.$configs{projectname}.$day.list" or die $!;
	%fqinputs = %bqsrbam if ($configs{alignment} =~ /false/i);
	for my $sampid(sort keys %fqinputs){
		for my $dtype(sort keys %{$fqinputs{$sampid}}){
			if (exists $scripts{clean}){        
				for my $lane(sort keys %{$scripts{clean}{$sampid}{$dtype}}){
					foreach my $num(sort keys %{$scripts{clean}{$sampid}{$dtype}{$lane}}){
						if($configs{split_num} > 1 && $splitfqs{$sampid}{$dtype}{$lane}{fqstat_exist} =~ /true/i){
							## split fq
							my $t = $num % $configs{split_num};
							my $zui = ($num - $t) / $configs{split_num}; 
							print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.1}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";
							print EDGE "$scripts{split}{$sampid}{$dtype}{$lane}{$zui.2}:1G:1CPU $scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU\n";   
						}
						## clean megaBOLT
						print EDGE "$scripts{clean}{$sampid}{$dtype}{$lane}{$num}:8G:1CPU $scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q}\n"; 
					}
			   }
			}

			## alignment and qc
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{megaboltqc}{$sampid}:10G:1CPU\n" if($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /true/i && $configs{tumor_only} =~ /false/i && -e $configs{qcvcf});
			print EDGE "$scripts{megaboltqc}{$sampid}:10G:1CPU\n" if ($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /false/i && $dtype =~ /Normal/ && $configs{tumor_only} =~ /false/i && -e $configs{qcvcf});
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{coverage}{$sampid}{$dtype}:10G:1CPU\n" if ($configs{megaboltqc} =~ /true/i && $configs{alignment} =~ /false/i);
			## microbiome analysis
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{metaphlan}{$sampid}{$dtype}:$configs{metaphlan_q}\n" if($configs{metaphlan} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{metaphlan}{$sampid}{$dtype}:$configs{metaphlan_q}\n" if($configs{metaphlan} =~ /true/i && $configs{alignment} =~ /false/i);
			print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{kraken2}{$sampid}{$dtype}:$configs{kraken2_q}\n" if($configs{kraken} =~ /true/i && $configs{alignment} =~ /true/i);
			print EDGE "$scripts{kraken2}{$sampid}{$dtype}:$configs{kraken2_q}\n" if($configs{kraken} =~ /true/i && $configs{alignment} =~ /false/i);                
			##hla dection
			if($configs{alignment} =~ /true/i){

			}
			if($configs{alignment} =~ /false/i){
				unless(-e $configs{hlalist}){
					## polysolver
					print EDGE "$scripts{polysolver}{$sampid}{$dtype}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
					## HLAminer
					print EDGE "$scripts{HLAminer}{$sampid}{$dtype}:$configs{HLAminer_q}\n" if($configs{hlaminer} =~ /true/i);
				}
			}
			else{
				unless(-e $configs{hlalist}){
					## polysolver
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{polysolver}{$sampid}{$dtype}:4G:1CPU\n" if($configs{hlasomatic} =~ /true/i);
					## HLAminer
					print EDGE "$scripts{MegaBOLT}{$sampid}{$dtype}:$configs{megaBOLT_q} $scripts{HLAminer}{$sampid}{$dtype}:$configs{HLAminer_q}\n" if($configs{hlaminer} =~ /true/i);
					
				}
			}
		}
	}
	close EDGE;
}