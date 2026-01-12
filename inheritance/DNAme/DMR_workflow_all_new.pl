#!usr/bin/perl -w


=head1 Version
Author: Xing Datong
Date: 2023-08-07
=head1 Usage
perl  DMR_workflow.pl methy1 methy2 [options]
-win                   <int>       input the window size,default 200
-depth                 <int>       setup the low depth cutoff,default 4
-pvalue1               <float>     setup the p value of the fisher test of DMR,default 0.01
-pvalue2               <float>     setup the p value of the fisher test of DmC,default 0.05
-gap                   <int>       setup the extend length of different DMR,default 100
-method                <int>       call methylation method,1 bsmap ,2 bismark,default 1
-fold_DMR              <float>     setup the foldchang of methylation regions between different samples,default 2
-fold_DMC              <float>     setup the foldchang of methylated C between different samples,default 2
-chr_len               <str>       input the chr length *.fai
-type                  <str>       the dmr found type,default(all), or each type(CG, CHG, CHH); all or each
-help                              output help information to screen
=head1 Exmple
perl DMR_workflow.pl methy1 methy2 id1 id2 -chr_len
=cut
use strict;
use Getopt::Long;
use FindBin qw($Bin $Script);
my ($win, $dep, $pvalue1, $pvalue2, $gap, $help, $fold, $chr, $method, $type1,$Bin_each,$fold_dmc);
GetOptions(
    "win:i"=>\$win,
    "depth:i"=>\$dep,
    "pvalue1:f"=>\$pvalue1,
    "pvalue2:f"=>\$pvalue2,
    "gap:i"=>\$gap,
    "fold_DMR:f"=>\$fold,
    "fold_DMC:f"=>\$fold_dmc,
    "method:i"=>\$method,
    "chr_len:s"=>\$chr,
    "type:s"=>\$type1,
    "help"=>\$help
);

$gap    ||=100;      $win    ||=200;    $dep    ||=4;$type1 ||="all";
$pvalue1    ||=0.01; $pvalue2   ||=0.05;$fold   ||=2;$method    ||=1;

die `pod2text $0` if (@ARGV < 4 || $help);
my $sa1=$ARGV[2];my $sa2=$ARGV[3];
$Bin="/mnt/chenbw/Baihua_Genome_2016/04.NMR/shell";
$Bin_each="/media/newdisks/HDD1/xingdt/bin/Methylation/DMR/each_type-new";
if($type1 eq "all"){
    open OUT,">$sa1\_$sa2\_dep$dep\_all.sh";
    #print OUT "#PBS -l nodes=1:ppn=20\n#PBS -q public01\n#PBS -l mem=10G\ncd \$PBS_O_WORKDIR\n";
    print OUT "perl $Bin/step1_identify_candidate_DMR_4step.pl $ARGV[0] $ARGV[1] $dep $win $chr $method > $sa1\_$sa2.dep$dep.win$win\n";
    print OUT "perl $Bin/step2_BH_correct_and_identify_DmC_4step.pl  $sa1\_$sa2.dep$dep.win$win $ARGV[0] $ARGV[1] $dep $pvalue1 $win $sa1\_$sa2.dep$dep.win$win.dmc  $pvalue2 $fold $method \n";
    print OUT "perl $Bin/step3_identify_hyper_and_hypo_jocob.pl $sa1\_$sa2.dep$dep.win$win.dmc 7 $gap $fold $ARGV[0] $ARGV[1] 4 > $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap\n";
######
	print OUT "cat $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap|sed -e '1d;s\/-\/\\t\/;s\/:\/\\t/' | awk -v OFS=\"\\t\" '{print \$1,\$2,\$3,\$4}'|sort -k1,1 -k2,2n > $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-DMR.bed\n";
	print OUT "grep \"hyper\" $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-DMR.bed > $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-hyper-DMR.bed\n";
	print OUT "grep \"hypo\" $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-DMR.bed > $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-hypo-DMR.bed\n";
	########验证最后合并后的DMR之间存在overlap的情况，是否都是因为该区域同时存在hyper和hypo DMR导致
	print OUT "bedtools intersect -a  $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-DMR.bed -b  $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-hyper-DMR.bed -wa -wb|awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}'|uniq -d|bedtools intersect -a - -b $sa1\_$sa2.dep$dep.win$win.dmc7.gap$gap-hypo-DMR.bed -wa -wb |awk -v OFS='\\t' '{print \$5,\$6,\$7,\$7-\$6,\$8}' >  $sa1\_$sa2-all-tC.hyper-hypo-overlap-DMR\n";
    close OUT;
}
if($type1 eq "each"){
    my @type=qw/CG CHG CHH/;
    open OUT,">$sa1\_$sa2\_dep$dep\_each.sh";
    foreach my $k(@type){
        my $type=$k;
        print OUT "perl $Bin_each/step1_identify_candidate_DMR_4step_each.pl $ARGV[0] $ARGV[1] $dep $win $type $method $chr > $sa1\_$sa2\_$type.dep$dep.win$win\n";
        print OUT "perl $Bin_each/step2_BH_correct_and_identify_DmC_4step_each.pl $sa1\_$sa2\_$type.dep$dep.win$win $ARGV[0] $ARGV[1] $dep $pvalue1 $win $sa1\_$sa2\_$type.dep$dep.win$win.dmc $pvalue2 $fold $type $method $fold_dmc\n";
        print OUT "perl $Bin_each/step3_identify_hyper_and_hypo_jocob_each.pl $sa1\_$sa2\_$type.dep$dep.win$win.dmc 7 $gap $ARGV[0] $ARGV[1] $fold $pvalue2 $type $method $sa1 $sa2 $fold $fold_dmc $dep > $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap\n";
	######
	print OUT "cat $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap|sed -e '1d;s\/-\/\\t\/;s\/:\/\\t/' | awk -v OFS=\"\\t\" '{print \$1,\$2,\$3,\$4}'|sort -k1,1 -k2,2n > $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed\n";
	print OUT "grep \"hyper\" $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed > $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-hyper-DMR.bed\n";
	print OUT "grep \"hypo\" $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed > $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-hypo-DMR.bed\n";
	########验证最后合并后的DMR之间存在overlap的情况，是否都是因为该区域同时存在hyper和hypo DMR导致
	print OUT "bedtools intersect -a $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed -b  $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed -wa -wb|awk '{print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$4}'|uniq -d|bedtools intersect -a - -b $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap-DMR.bed -wa -wb |awk -v OFS='\\t' '{print \$5,\$6,\$7,\$7-\$6,\$8}' > $sa1\_$sa2\_$type.dep$dep.$type.win$win.dmc7.gap$gap.hyper-hypo-overlap-DMR\n";
    }
    close OUT;
}
if($type1 ne "all" && $type1 ne "each"){
    print "error\tplease input the type , all or each\n";
}
