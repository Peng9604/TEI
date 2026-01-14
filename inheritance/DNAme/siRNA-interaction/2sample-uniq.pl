#!/usr/bin/perl -w
use strict;
die"usage:perl $0 <sample1-uniq.bed> <sample2-uniq.bed> <sample1> <sample2> <OUT-prefix>\n"unless @ARGV==5;
#####当大小区域在一个文件中，多个小区域对应相同的大区域。统计两个样本的相同大区域中，各自存在多少个小区域。共两个文件，一个文件是sample1的每个大区域对应的小区域，另一个文件是sample2的每个大区域对应的小区域。每个大区域必须有名字（如CG_NI-1）和样本名（如Col）。文件必须包括：chr1 1000 1400(大区域位置，其实有名字后，可以没有大区域的坐标) CG_NI-1 Col chr1 1020 1044（小区域位置）
open IN1,$ARGV[0];
open IN2,$ARGV[1];
my %tcm;
while(<IN1>){
	chomp;my @a=split/\s+/;
	my $id=join "\t",@a[0..2];
	push @{$tcm{$a[4]}{$a[3]}},$id;  ###key1是NI编号、key2是Col/C24，value是一个数组，包括符合key1、key2的所有24nt siRNA的比对位置
}
close IN1;

while(<IN2>){
	chomp;my @a=split/\s+/;
	my $id=join "\t",@a[0..2];
	push @{$tcm{$a[4]}{$a[3]}},$id;
}
close IN2;


open OUT1,">$ARGV[4]-2sample-uniq_number.txt";
open OUT2,">$ARGV[4]-$ARGV[2]-uniq.bed5";
open OUT3,">$ARGV[4]-$ARGV[3]-uniq.bed5";

print OUT1 "TCM\t$ARGV[2]_$ARGV[3]\t$ARGV[2]\t$ARGV[3]\n";
foreach my $k1(sort keys %tcm){
        my $samp1_num;my $samp2_num;
        #my @sample=sort keys %{$tcm{$k3}};
        if(exists $tcm{$k1}{$ARGV[2]}){
            $samp1_num=scalar(@{$tcm{$k1}{$ARGV[2]}});
        }else{
            $samp1_num=0
        }
        if(exists $tcm{$k1}{$ARGV[3]}){
	     $samp2_num=scalar(@{$tcm{$k1}{$ARGV[3]}});
        }else{
            $samp2_num=0
        }
	next if $samp1_num==0;
	next if $samp2_num==0;
        my $all_num=$samp1_num+$samp2_num;
	print OUT1 "$k1\t$all_num\t$samp1_num\t$samp2_num\n";
	foreach my $k2(@{$tcm{$k1}{$ARGV[2]}}){
		print OUT2 "$k2\t$ARGV[2]-$k1\n";
	}
	foreach my $k3(@{$tcm{$k1}{$ARGV[3]}}){
		print OUT3 "$k3\t$ARGV[3]-$k1\n";
	}
}
