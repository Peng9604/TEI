#!usr/bin/perl -w
####1.统计meth1中每个C的tC、mC、mC/tC，过滤掉tC小于4（命令行可以设置）的C
####2.符合第一项条件的C，将这个C的tC、mC、mC/tC累加给该C所在bin为起始的window和上游的三个window（每个window的起始bin之差均为50bp，meth2同理）。这样，便实现了200bp/window，50bp/步移。但是如果该C在一个样本中有覆盖，而在另一个样本中没有覆盖，那么在这个脚本中会被去除。
####3.过滤掉meth1和meth2都没有合格覆盖度C的window 。接着利用fishier检验，计算每个window的meth1和meth2的甲基化差异程度，返回p值和FDR；并且利用这个window的meth1和meth2的平均甲基化水平（每个C的甲基化水平累加除以C的个数）计算FC（meth1和meth2哪个大，哪个当分子）
####4.输出文件（只要meth1和meth2至少有一个合格C（tC>4）的window都被保留，无论是否鉴定为潜在dmr）：染色质编号，window的编号，meth1在这个window的tC、mC、C的数量，meth2在这个window的tC、mC、C的数量，meth1和meth2在这个window中甲基化的FC，这个window中meth1和meth2甲基化差异的p值，矫正过的FDR
####5.因为该脚本没有卡每个dmr中dmc的数量，所以需要通过step2脚本
######输入：第一个和第二个参数分别是比较的两个样本的methyltable
    #第三个：C的覆盖度阈值，一般为4
    #第四个：window大小，一般为200
    #第五个：CG/CHG/CHH
    #第六个：1
    #第七个：基因组每条染色质的总长度
use strict;
use Statistics::Descriptive;
use Statistics::Distributions;
use Text::NSP::Measures::2D::Fisher::twotailed;
die "usage:perl $0 <methy1.txt> <methy2.txt> <low depth> <window> <methy_type> <type:1-bsmap/2-bismark> <chr length>\n" unless @ARGV == 7;

open A,$ARGV[0];open B,$ARGV[1];open LEN,$ARGV[6];
my (%methy1, %methy2, %type ,$seq, $methy_seq, %num_methy1, %num_methy2);my $step=$ARGV[3]/4;  my %chr_len;   #$step是步移，是window大小的1/4

my %tr=(
    'CGA'=>"CG",'CGT'=>"CG",'CGC'=>"CG",'CGG'=>"CG",'CAG'=>"CHG",'CTG'=>"CHG",'CCG'=>"CHG"
);
while(<LEN>){
    chomp;my @a=split;    $chr_len{$a[0]}=$a[1];    #key是染色质编号，value是染色质长度
}
close LEN;
if($ARGV[5] == 1){
    <A>;<B>;
}
while(<A>){
    chomp;my @a=split;
    my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
    if($ARGV[5] == 2){
        $mdepth=$a[3];$tdepth=$a[3]+$a[4]; 
        next if($tdepth == 0);		
		$rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
    }
	
    if($ARGV[5] == 1){
		next if $tdepth < $ARGV[2];   #排除该位点测序深度不达标的（reads覆盖度小于4）
        $methy_seq=&identify_methy_type($a[2],$a[3]);
    }
	if($ARGV[4] ne "C"){
	    next if($methy_seq ne $ARGV[4]);
	}
    my $sub_window=int($a[1]/$step);    ###以C所在的bin为首的window的编号：进入本次循环的C的基因组位置除以步移（50bp）的结果，取整。
    $methy1{$a[0]}{$sub_window}[0]+=$tdepth;   $methy1{$a[0]}{$sub_window-1}[0]+=$tdepth;    #key1是进入这次循环的C的染色质编号，key2是这个C的位置所在的bin的编号（全基因组范围给bin编号），除在该bin中累加这个C的信息（即将该C的信息累加到以该bin为首的window）外，还将这个C累加到分别以上游3个bin为首的window中（这就是200bp/window，50bp/步移）。value1是以这个bin为首的window的totalC，value2是以这个bin为首的window的mC，value3是每进入循环一个C，就增加1（代表以这个bin为首的window的总C个数），value4是以这个bin为首的mC/tC
    $methy1{$a[0]}{$sub_window-2}[0]+=$tdepth; $methy1{$a[0]}{$sub_window-3}[0]+=$tdepth;

    $methy1{$a[0]}{$sub_window}[1]+=$mdepth;   $methy1{$a[0]}{$sub_window-1}[1]+=$mdepth;
    $methy1{$a[0]}{$sub_window-2}[1]+=$mdepth; $methy1{$a[0]}{$sub_window-3}[1]+=$mdepth;

    $methy1{$a[0]}{$sub_window}[3]+=$rate;   $methy1{$a[0]}{$sub_window-1}[3]+=$rate;
    $methy1{$a[0]}{$sub_window-2}[3]+=$rate; $methy1{$a[0]}{$sub_window-3}[3]+=$rate;

    $methy1{$a[0]}{$sub_window}[2]++;         $methy1{$a[0]}{$sub_window-1}[2]++;
    $methy1{$a[0]}{$sub_window-2}[2]++;       $methy1{$a[0]}{$sub_window-3}[2]++;
}
close A;
while(<B>){
    chomp;my @a=split;
    my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
    if($ARGV[5] == 2){
        $mdepth=$a[3];$tdepth=$a[3]+$a[4];  
next if($tdepth == 0);        $rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
    }
	
    if($ARGV[5] == 1){
		next if $tdepth < $ARGV[2];   #排除该位点测序深度不达标的（reads覆盖度小于4）
        $methy_seq=&identify_methy_type($a[2],$a[3]);
    }
	if($ARGV[4] ne "C"){
	    next if($methy_seq ne $ARGV[4]);
	}
	my $sub_window=int ($a[1]/$step);
    $methy2{$a[0]}{$sub_window}[0]+=$tdepth;   $methy2{$a[0]}{$sub_window-1}[0]+=$tdepth;
    $methy2{$a[0]}{$sub_window-2}[0]+=$tdepth; $methy2{$a[0]}{$sub_window-3}[0]+=$tdepth;

    $methy2{$a[0]}{$sub_window}[1]+=$mdepth;   $methy2{$a[0]}{$sub_window-1}[1]+=$mdepth;
    $methy2{$a[0]}{$sub_window-2}[1]+=$mdepth; $methy2{$a[0]}{$sub_window-3}[1]+=$mdepth;

    $methy2{$a[0]}{$sub_window}[3]+=$rate;   $methy2{$a[0]}{$sub_window-1}[3]+=$rate;
    $methy2{$a[0]}{$sub_window-2}[3]+=$rate; $methy2{$a[0]}{$sub_window-3}[3]+=$rate;

    $methy2{$a[0]}{$sub_window}[2]++;         $methy2{$a[0]}{$sub_window-1}[2]++;
    $methy2{$a[0]}{$sub_window-2}[2]++;       $methy2{$a[0]}{$sub_window-3}[2]++;
}
close B;

my (%all_info, %dmr_select, %dmr_info);

foreach my $k1(sort keys %chr_len){
    my $window_len=int ($chr_len{$k1}/$step);   #$window_len代表染色质共被分成多少个window（200bp/window，50bp/步移）
    foreach my $k2(0..$window_len){     #####从第一个window开始计算
        next if(!exists $methy1{$k1}{$k2} || !exists $methy2{$k1}{$k2});   #过滤掉在meth1或meth2中没有合格的C的window
        my $npp = $methy1{$k1}{$k2}[0]+$methy2{$k1}{$k2}[0]; my $n1p = $methy1{$k1}{$k2}[0];      #$npp是这个window中methytable1的tC和methytable2的tC之和，$n1p是这个window中methytable1的tC，$np1这个window中methytable1的mC和methytable2的mC之和，$n11是这个window中methytable1的mC
        my $np1 = $methy1{$k1}{$k2}[1]+$methy2{$k1}{$k2}[1]; my $n11 = $methy1{$k1}{$k2}[1];
        my $fisher_exact_test_twotail=calculateStatistic(n11=>$n11, n1p=>$n1p, np1=>$np1, npp=>$npp);
        if(!defined $fisher_exact_test_twotail){
            $fisher_exact_test_twotail=1;
        }
        my $rate;
		if($methy1{$k1}{$k2}[1]/$methy1{$k1}{$k2}[0] > $methy2{$k1}{$k2}[1]/$methy2{$k1}{$k2}[0]){  #判断methytable1和methytable2在这个window中谁的平均甲基化水平更高：用这个window中所有C的甲基化率之和除以这window中C的总数
            $rate=($methy1{$k1}{$k2}[1]/$methy1{$k1}{$k2}[0]+0.01)/($methy2{$k1}{$k2}[1]/$methy2{$k1}{$k2}[0]+0.01);    #计算methytable1和methytable2在这个window中甲基化水平的FC：利用在这个window中甲基化水平高的methytable的平均mC/tC除以甲基化水平低的methytable的平均mC/tC（都加了0.01是害怕meth2在这个window中所有C甲基化率加和是0）。
        }
        else{
            $rate=($methy2{$k1}{$k2}[1]/$methy2{$k1}{$k2}[0]+0.01)/($methy1{$k1}{$k2}[1]/$methy1{$k1}{$k2}[0]+0.01);
        }
        #if($methy1{$k1}{$k2}[3]/$methy1{$k1}{$k2}[2] > $methy2{$k1}{$k2}[3]/$methy2{$k1}{$k2}[2]){  #判断methytable1和methytable2在这个window中谁的平均甲基化水平更高：用这个window中所有C的甲基化率之和除以这window中C的总数
        #    $rate=(($methy1{$k1}{$k2}[3]+0.01)/$methy1{$k1}{$k2}[2])/(($methy2{$k1}{$k2}[3]+0.01)/$methy2{$k1}{$k2}[2]);    #计算methytable1和methytable2在这个window中甲基化水平的FC：利用在这个window中甲基化水平高的methytable的平均mC/tC除以甲基化水平低的methytable的平均mC/tC（都加了0.01是害怕meth2在这个window中所有C甲基化率加和是0）。
        #}
        #else{
        #    $rate=(($methy2{$k1}{$k2}[3]+0.01)/$methy2{$k1}{$k2}[2])/(($methy1{$k1}{$k2}[3]+0.01)/$methy1{$k1}{$k2}[2]);
        #}
        my $temp1=join "\t",@{$methy1{$k1}{$k2}}[0..2];
        my $temp2=join "\t",@{$methy2{$k1}{$k2}}[0..2];
        $dmr_info{$k1}{$k2}=$fisher_exact_test_twotail;   #%dmr_info的key1是染色质，key2是window的基因组编号，value是meth1和meth2在这个window的甲基化水平差异的p值
        $all_info{$k1}{$k2}=$temp1."\t".$temp2."\t".$rate."\t".$fisher_exact_test_twotail;   #%all_info的key1是染色质，key2是window的基因组编号，value是一个组合出来的字符串：meth1在这个window的tC、mC、C的数量，meth2在这个window的tC、mC、C的数量，meth1和meth2在这个bin中甲基化的FC，这个window中meth1和meth2甲基化差异的p值
    }
}
foreach my $k1(sort keys %dmr_info){
    my %pvalue_correct=&BH_correct(\%{$dmr_info{$k1}});
    foreach my $k2(sort {$a<=>$b} keys %{$dmr_info{$k1}}){
        print "$k1\t$k2\t$all_info{$k1}{$k2}\t$pvalue_correct{$k2}\n";   #染色质编号，window的编号，meth1在这个window的tC、mC、C的数量，meth2在这个window的tC、mC、C的数量，meth1和meth2在这个window中甲基化的FC，这个window中meth1和meth2甲基化差异的p值，矫正过的FDR
    }
}
sub BH_correct{
    my $record=0; my $initial=shift; my %result;
    my @sort_pvalue=sort {$a<=>$b} values %$initial;
    my @sort_site=sort {$$initial{$a} <=> $$initial{$b}} keys %$initial;
    my $num=$#sort_pvalue+1;$record=$sort_pvalue[-1];
    $result{$sort_site[-1]}=$sort_pvalue[-1];
    for(my $i=$#sort_pvalue-1;$i >= 0;$i--){
        my $site=$i+1;
        my $temp=$sort_pvalue[$i]*$num/$site;
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data($record,$temp);
        my $min = $stat->min();
        $result{$sort_site[$i]}=$min;
    }
    return %result;
}

sub identify_methy_type{
    my $direction=shift;
    my $file=shift;
    my ($seq, $methy_seq);
    if($direction eq "-"){
        $seq=substr($file,0,3);
        $seq=reverse($seq);
        $seq=~tr/ATCG/TAGC/;
    }
    else{
        $seq=substr($file,2);
    }
    if(exists $tr{$seq}){
        $methy_seq=$tr{$seq};
    }
    else{
        $methy_seq="CHH";
    }
    return $methy_seq;
}
