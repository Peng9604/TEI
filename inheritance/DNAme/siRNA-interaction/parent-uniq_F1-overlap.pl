#!/usr/bin/perl -w
use strict;
die"usage:perl $0 <sample1-uniq.bed> <sample2-uniq.bed> <sample1> <sample2> <index> <OUT-prefix>\n"unless @ARGV==6;

open IN1,$ARGV[0];
open IN2,$ARGV[1];
my %tcm;
while(<IN1>){
	chomp;my @a=split/\s+/;my $str;
        if($_=~/.*($ARGV[4].*)$/){
              $str=$1;
#              print "$str\n";
        }
	my $id=join "\t",@a[0..2];
	push @{$tcm{$str}{$ARGV[2]}},$id;
}
close IN1;

while(<IN2>){
	chomp;my @a=split/\s+/;my $str;
	if($_=~/.*($ARGV[4].*)$/){
             $str=$1;
        }
	my $id=join "\t",@a[0..2];
	push @{$tcm{$str}{$ARGV[3]}},$id;
}
close IN2;

open OUT1,">$ARGV[5]-parent-uniq_F1-overlap-number.txt";
open OUT2,">$ARGV[5]-$ARGV[2]-qualified-uniq.bed5";
open OUT3,">$ARGV[5]-$ARGV[3]-qualified-uniq.bed5";

print OUT1 "TCM\t$ARGV[2]_$ARGV[3]\t$ARGV[2]\t$ARGV[3]\n";
foreach my $k1(sort keys %tcm){
        my $samp1_num;my $samp2_num;
        #my @sample=sort keys %{$tcm{$k3}};
        if(exists $tcm{$k1}{$ARGV[2]}){
            $samp1_num=scalar(@{$tcm{$k1}{$ARGV[2]}});
        }else{
            $samp1_num=0;
        }
        if(exists $tcm{$k1}{$ARGV[3]}){
	     $samp2_num=scalar(@{$tcm{$k1}{$ARGV[3]}});
        }else{
            $samp2_num=0;
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
