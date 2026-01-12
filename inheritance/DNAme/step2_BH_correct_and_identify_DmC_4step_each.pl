#!usr/bin/perl -w
####该脚本主旨：1. 保留符合DMR条件的window（200bp）作为候选DMR；2. 计算候选DMR的DMC数量
#1.该脚本首先对step1输出文件中，只有meth1和meth2在这个window上的甲基化差异FDR小于或等于DMR pvalue cutoff、并且meth1和meth2在这个window上的甲基化差异FC大于2（命令行可以修改）的window被保留下来。
#2.计算DMC的C的基本要求：1.这个C在两个样本中都有大于4个reads覆盖;2. 这个C落在以上过滤后得到的window中；3. 这个C在至少一个样本的mC不为0。
#3.判断这个C是DMC的标准：1. 这个C在其中一个样本的mC为0时，另一个样本的甲基化水平大于该context type的差异阈值（CG 0.4，CHG 0.2，CHH 0.1）时，判定这个C为DMC；2.这个C在meth1和meth2的mC都不是0，需要满足如下三个条件才能被判为DMC：1.fishier检验得到的p值矫正后小于阈值；2. 这个C在两个样本上的甲基化水平之差大于该context type的差异阈值（CG 0.4，CHG 0.2，CHH 0.1）；3.两个样本的甲基化水平之比（foldchange）大于阈值
#3.一旦这个C被判断为dmC，那么会被累加到这个C所在bin为首的window和上游三个bin为首的window上（只累加到没有被过滤掉的DMR中）。并且最后会输出符合DMR条件的window，也就是并不像第一个文件把所有的window都输出出来。
#####输出文件（每一行都是一个符合dmr条件的window）：染色质编号：每个window的start和end（长度200bp），这个window中meth1甲基化水平低的dmC，这个window中meth2甲基化水平低的dmC，meth1在这个window的tC、mC、总C的数量，meth2在这个window的tC、mC、总C的数量，meth1和meth2在这个window中甲基化的FC，这个window中meth1和meth2甲基化差异的p值，矫正过的FDR
#####200bp/window，50bp/步移的思路：我们以window=200bp，步移=50bp。每个bin可以作为四个window的步移，这个bin的甲基化信息会累加到这四个window中，这四个window的起始bin分别是该bin、该bin- 1、该bin-2、该bin-3. 同样，每个bin都是一个window的起始，这个bin（即这个window）都包含了下游三个bin的甲基化状态。所以每个window的编号就是起始bin的编号。
#####输入：第一个：上一步输入的候选DMR文件。第二个和第三个是两个样本的methtable，第四个是C的覆盖度阈值（一般为4），第五个是DMR的FDR阈值，第六个是200bp的窗口大小，第七个是输出文件名，第八个是DMC的FDR阈值，第九个是过滤第一步候选DMR foldchange的阈值，第十个是CG/CHG/CHH，第十一个是1，第十二个是DMC的foldchange阈值
use strict;
use Statistics::Descriptive;
use Statistics::Distributions;
use Text::NSP::Measures::2D::Fisher::twotailed;

die "usage:perl $0 <candidate DMR info> <methy1> <methy2> <low depth cutoff> <DMR pvalue cutoff> <window size> <output> <DmC pvalue cutoff> <第一步计算的DMR的foldchange cutoff> <methy type> <type 1-bsmap 2-bismark> <DMC foldchange cutoff>\n" unless @ARGV == 12;

#第一个输入文件（上一步的输出）：染色质编号，window的编号，meth1在这个window的tC、mC、C的数量，meth2在这个window的tC、mC、C的数量，meth1和meth2在这个window中甲基化的FC，这个window中meth1和meth2甲基化差异的p值，矫正过的FDR



open A,$ARGV[0];open B,$ARGV[1];open C,$ARGV[2];
open OUT1,">$ARGV[6]";
my %tr=(
	'CGA'=>"CG",'CGT'=>"CG",'CGC'=>"CG",'CGG'=>"CG",'CAG'=>"CHG",'CTG'=>"CHG",'CCG'=>"CHG"
);
my (%dmr_info, %methy1, %methy2, %dmr_select, %dmc_num, %each_win_num);
my $step=$ARGV[5]/4;   #200/4（50bp/步移）

while(<A>){
	chomp;    my @a=split;
	if($a[-1] <= $ARGV[4] && $a[8] > $ARGV[8]){   #meth1和meth2在这个window上的甲基化差异FDR小于或等于DMR pvalue cutoff，并且meth1和meth2在这个window上的甲基化差异FC大于2（命令行可以修改）的window被保留下来
		$dmr_select{$a[0]}{$a[1]}=join "\t",@a[2..$#a]; 
		$dmc_num{$a[0]}{$a[1]}[0]=0;$dmc_num{$a[0]}{$a[1]}[1]=0;   #%dmr_select是通过DMR的FRD阈值和FC阈值过滤后留下的window
	}
}
close A;
if($ARGV[10] == 1){
	<C>;<B>;
}
while(<B>){
	chomp;my @a=split;my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
	if($ARGV[10] == 2){
		$mdepth=$a[3]; $tdepth=$a[3]+$a[4];        
		next if($tdepth == 0);        
		$rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
	}
	if($ARGV[10] == 1){
		$methy_seq=&identify_methy_type($a[2],$a[3]);
	}
	next if($tdepth < $ARGV[3]);     #####过滤掉覆盖度不足4（命令行可修改）的C
	if($ARGV[9] ne "C"){
		next if($methy_seq ne $ARGV[9]);
	}
	my $win=int($a[1]/$step);
	if(exists $dmr_select{$a[0]}{$win}){        $methy1{$a[0]}{$a[1]}[0]=$tdepth;    $methy1{$a[0]}{$a[1]}[1]=$mdepth;    }  #重新计算过滤后每个window的C的甲基化状态：这个C的信息会同时被归于所在bin为首的window以及上游三个bin为首的window中（如果window被过滤掉了，那么不会计算这个window的情况）
	if(exists $dmr_select{$a[0]}{$win-1}){      $methy1{$a[0]}{$a[1]}[0]=$tdepth;    $methy1{$a[0]}{$a[1]}[1]=$mdepth;    }
	if(exists $dmr_select{$a[0]}{$win-2}){      $methy1{$a[0]}{$a[1]}[0]=$tdepth;    $methy1{$a[0]}{$a[1]}[1]=$mdepth;    }
	if(exists $dmr_select{$a[0]}{$win-3}){      $methy1{$a[0]}{$a[1]}[0]=$tdepth;    $methy1{$a[0]}{$a[1]}[1]=$mdepth;    }
}
close B;
my %type=(
	"CG"=>'0.4',
	"CHG"=>'0.2',
	"CHH"=>'0.1',
);
while(<C>){
	chomp;my @a=split;  my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
	if($ARGV[10] == 2){
		$mdepth=$a[3];$tdepth=$a[3]+$a[4];
		next if($tdepth == 0);
		$rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
	}
	if($ARGV[10] == 1){
		$methy_seq=&identify_methy_type($a[2],$a[3]);
	}
	next if($tdepth < $ARGV[3]);
	if($ARGV[9] ne "C"){
		next if($methy_seq ne $ARGV[9]);
	}
	my $win=int ($a[1]/$step);
	if(exists $methy1{$a[0]}{$a[1]}){      ####只有这个C在两个样本中都有大于4个reads覆盖，才会往下走
		next if($methy1{$a[0]}{$a[1]}[1] == 0 && $mdepth == 0);  #如果这个C的meth1和meth2的mC都为0，则过滤掉这个C
		my $record=0;
		if($methy1{$a[0]}{$a[1]}[1] == 0 && $mdepth/$tdepth > $type{$methy_seq}){  #如果这个C上，meth1的mC为0，但是meth2上的甲基化水平大于指定context的差异水平阈值，那么这个C就是符合条件的dmC，则这个C所在的bin和上游的三个bin上的dmC数量都加1
			$dmc_num{$a[0]}{$win}[0]++;   $dmc_num{$a[0]}{$win-1}[0]++;	$record=1;
			$dmc_num{$a[0]}{$win-2}[0]++; $dmc_num{$a[0]}{$win-3}[0]++;
		}
		if($mdepth == 0 && $methy1{$a[0]}{$a[1]}[1]/$methy1{$a[0]}{$a[1]}[0] > $type{$methy_seq}){
			$dmc_num{$a[0]}{$win}[1]++;   $dmc_num{$a[0]}{$win-1}[1]++;	$record=1;
			$dmc_num{$a[0]}{$win-2}[1]++; $dmc_num{$a[0]}{$win-3}[1]++;
		}
		next if($record == 1);  #已经被判定为dmC的C不会接着运行
		if($mdepth != 0 && $methy1{$a[0]}{$a[1]}[1] != 0){
			my $npp = $methy1{$a[0]}{$a[1]}[0]+$tdepth; my $n1p =$methy1{$a[0]}{$a[1]}[0];
			my $np1 = $methy1{$a[0]}{$a[1]}[1]+$mdepth; my $n11 =$methy1{$a[0]}{$a[1]}[1];
			my $fisher_exact_test_twotail=calculateStatistic( n11=>$n11, n1p=>$n1p, np1=>$np1, npp=>$npp);  #如果这个C上，meth1和meth2的mC都不是0，需要通过fishier检验来判断这个dmC是否有差异
			if(defined $fisher_exact_test_twotail && $fisher_exact_test_twotail <= $ARGV[7]){     ###判断dmC的p值是否存在，且判断这个C的差异p值是否小于或等于dmC的阈值，$ARGV[7]是dmC的p值的阈值
				my $res1=abs($methy1{$a[0]}{$a[1]}[1]/$methy1{$a[0]}{$a[1]}[0]-$mdepth/$tdepth);    ###$res1是两个样本甲基化水平之差
				if($methy_seq eq "CG"){					next if($res1 < 0.4);			     }   ###只保留两个样本甲基化水平差异符合指定context阈值的候选dmc（通过fishier检验p值的阈值的dmc）
				if($methy_seq eq "CHG"){				next if($res1 < 0.2);			     }
				if($methy_seq eq "CHH"){				next if($res1 < 0.1);			     }
				if($methy1{$a[0]}{$a[1]}[1]/$methy1{$a[0]}{$a[1]}[0] <= $mdepth/$tdepth){        ##
					my $rate=$methy1{$a[0]}{$a[1]}[1]/$methy1{$a[0]}{$a[1]}[0];
					my $foldchange=($mdepth/$tdepth)/$rate;
					next if($foldchange < $ARGV[11]);     ####在前面过滤掉不符合dmc FDR和差异阈值之后的候选dmc，还需要经过FC小于2的过滤才能最终留下。
					$dmc_num{$a[0]}{$win}[0]++;$dmc_num{$a[0]}{$win-1}[0]++;
					$dmc_num{$a[0]}{$win-2}[0]++;$dmc_num{$a[0]}{$win-3}[0]++;
				}
				else{
					my $rate=$mdepth/$tdepth;
					my $foldchange=($methy1{$a[0]}{$a[1]}[1]/$methy1{$a[0]}{$a[1]}[0])/$rate;
					next if($foldchange < $ARGV[11]);
					$dmc_num{$a[0]}{$win}[1]++;$dmc_num{$a[0]}{$win-1}[1]++;
					$dmc_num{$a[0]}{$win-2}[1]++;$dmc_num{$a[0]}{$win-3}[1]++;
				}
			}
		}
	}
}
close C;
foreach my $k1(sort keys %dmr_select){
	foreach my $k2(sort {$a<=>$b} keys %{$dmr_select{$k1}}){
		my $start=$k2*$step;my $end=$start+$ARGV[5]-1;
		print OUT1 "$k1:$start-$end\t$dmc_num{$k1}{$k2}[0]\t$dmc_num{$k1}{$k2}[1]\t$dmr_select{$k1}{$k2}\n";        #####染色质编号:这个bin往下游200bp的start和end，这个bin中meth1甲基化水平低的dmC个数，这个bin中meth2甲基化水平低的dmC个数，（从第四列到最后一列是第一个输出文件筛选DMR时获得的）meth1在这个bin的tC、mC、C的数量，meth2在这个bin的tC、mC、C的数量，meth1和meth2在这个bin中甲基化的FC，这个bin中meth1和meth2甲基化差异的p值，矫正过的FDR
	}
}
close OUT1;

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
