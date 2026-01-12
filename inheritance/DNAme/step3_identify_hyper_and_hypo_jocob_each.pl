#!usr/bin/perl -w

use strict;
use Text::NSP::Measures::2D::Fisher::twotailed;
die "usage:perl $0 <dmc> <DMR hyper或hypo dmc个数的cutoff> <merge length> <methy1> <methy2> <这个window（没合并前）被判定为DMR的foldchange cutoff（第一步fishier检验得到的）> <该脚本重算DmC pvalue cutoff> <methy_type> <type 1-bsmap 2-bismark> <id1> <id2> <合并后DMR-fold change阈值> <该脚本重算DMC foldchange cutoff> <C的覆盖度阈值，exp: 4>\n" unless @ARGV == 14;

#1.DMR（200bp）得到保留的条件：1.foldchange高于阈值 2.hyper或hypo DMC个数高于阈值
#2.将合格的DMR（200bp）进行合并，gap一般设置为100（输入命令可调整）
#3. 计算合并后DMR中每个C是否为DMC，符合要求的C具有以下特征：1.这C在两个样本中的覆盖度都高于4；2.这个C在候选DMR中；3.这个C至少在一个样本中的mC不等于0
#第一个输入文件（每一行都是一个符合dmr条件的window）：染色质编号:每个window的start-end（长度200bp），这个window中meth1甲基化水平低的dmC，这个window中meth2甲基化水平低的dmC，meth1在这个window的tC、mC、总C的数量，meth2在这个window的tC、mC、总C的数量，meth1和meth2在这个window中甲基化的FC，这个window中meth1和meth2甲基化差异的p值，矫正过的FDR
#####输出文件：染色质编号:每个DMR的start-end，hyper/hypo,这个DMR中meth2大于meth1甲基化水平(hyper)的dmC个数，这个DMR中meth2小于meth1甲基化水平（hypo）的dmC个数，这个DMR上的总C的个数（这些C上至少有一个样本的甲基化水平大于0），这个DMR上meth1的甲基化水平，这个DMR上meth2的甲基化水平，这个DMR的log2FC
#####输入文件：1.上一步输出的每个合格DMR的DMC信息，2.DMR hyper或hypo dmc个数的阈值，3.两个DMR合并的gap阈值，一般为100bp，4.样本1的methytable，5.样本2的methytable，6.window（没合并前）被判定为DMR的foldchange 阈值（第一步fishier检验得到的），6.该脚本重算DmC时FDR的阈值（重算原因是合并后会有原来没被计算的gap区域），7.CG/CHG/CHH，8.设置为1，10.样本1的名字，11.样本2的名字，12.合并后DMR foldchange的阈值，13.该脚本重算DMC foldchange的阈值，14.C的覆盖度阈值，一般为4
my (%info, %direction, %all_info, $column_num, %merge_result, %res);

my %tr=(
	'CGA'=>"CG",'CGT'=>"CG",'CGC'=>"CG",'CGG'=>"CG",'CAG'=>"CHG",'CTG'=>"CHG",'CCG'=>"CHG"
);
open A,$ARGV[0];open IN1,$ARGV[3];open IN2,$ARGV[4];
while(<A>){
	chomp;my @a=split;	my $inf=join "\t",@a[1..$#a];	my @b=split /[:-]/,$a[0];next if($a[9] <= $ARGV[5]); 		next if($a[1] < $ARGV[1] && $a[2] < $ARGV[1]);  ###DMR的FC大于阈值且至少保证hyper DMC或hypo DMC个数达到阈值才能被留下
	if($a[7]/$a[6] > $a[4]/$a[3]){
		$info{"hyper"}{$b[0]}{$b[1]}=$b[2];    #####key1:hyper/hypo，key2时chr，key3是这个window的start，value是这个window的end（比start大200）
	}
	if($a[7]/$a[6] < $a[4]/$a[3]){
		$info{"hypo"}{$b[0]}{$b[1]}=$b[2];
	}
}
close A;
##%info：储存未合并的DMR的信息。key1是hyper/hypo，key2是chr，key3是start，value是end
foreach my $k3(sort keys %info){     ##$k3是hyper/hypo
	foreach my $k1(sort keys %{$info{$k3}}){  ###$k1是chr
		my @site=sort {$a<=>$b} keys %{$info{$k3}{$k1}};	##@site储存了所有DMR的start（按基因组顺序）
		my $record=0;	my (%tar, %merge_info);	$tar{$record}[0]=$site[0];      ####%tar：key是合并后DMR的序号；value1是合并后新DMR的start，value2是新DMR的end
		foreach my $k2(1..$#site){
			if($info{$k3}{$k1}{$site[$k2-1]} >= $site[$k2]-$ARGV[2]){
				$tar{$record}[1]=$info{$k3}{$k1}{$site[$k2]};   ###重叠的最后一个DMR的end作为合并后DMR的新end
			}
			else{
				$tar{$record}[1]=$info{$k3}{$k1}{$site[$k2-1]};  
				$record++;
				$tar{$record}[0]=$site[$k2];   ###刷新下一个start
			}
		}
		$tar{$record}[1]=$info{$k3}{$k1}{$site[$#site]};   ###封住该染色质最后一个end
		foreach my $k2(sort {$a<=>$b} keys %tar){
			my $id=$k1.":".$tar{$k2}[0]."-".$tar{$k2}[1];
			foreach my $kk($tar{$k2}[0]..$tar{$k2}[1]){
				$res{$k1}{$kk}=$id;   ####key1染色质号，key2是合并后的候选DMR的每个碱基位置，value是这个候选DMR的位置信息（chr:start-end）
			}
			$merge_result{$id}[0]=0;$merge_result{$id}[1]=0;$merge_result{$id}[2]=0;$merge_result{$id}[3]=$k3;
		}
	}
}
my (%methy1, %methy2, %methy_info1);
if($ARGV[8] == 1){
	<IN1>;<IN2>;
}
while(<IN1>){
	chomp;my @a=split;my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
	if($ARGV[8] == 2){
		$mdepth=$a[3];$tdepth=$a[3]+$a[4];		
		next if($tdepth == 0);$rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
	}
	if($ARGV[8] == 1){
		$methy_seq=&identify_methy_type($a[2],$a[3]);
	}
	next if($tdepth <$ARGV[13]);   #####next if($tdepth < $ARGV[3]);
	if($ARGV[7] ne "C"){
		next if($methy_seq ne $ARGV[7]);
	}
	if(exists $res{$a[0]}{$a[1]}){    ####如果这个C在合并候选的DMR中，则储存这个C的mC和tC（%methy_info1）且将这个C的mC和tC储存在这个合并候选DMR中（methy1）
		$methy_info1{$a[0]}{$a[1]}[0]=$tdepth;     $methy_info1{$a[0]}{$a[1]}[1]=$mdepth;   ###储存这个C的信息
		$methy1{$res{$a[0]}{$a[1]}}[0]+=$tdepth;   $methy1{$res{$a[0]}{$a[1]}}[1]+=$mdepth;   ####储存这个C属于的新DMR的信息
	}
}
close IN1;
my %type=(
	"CG"=>'0.4',
	"CHG"=>'0.2',
	"CHH"=>'0.1',
);
while(<IN2>){
	chomp;my @a=split;    my ($mdepth, $tdepth, $rate, $methy_seq)=($a[-1],$a[-2],$a[-3],$a[-4]);
	if($ARGV[8] == 2){
		$mdepth=$a[3];$tdepth=$a[3]+$a[4];		
		next if($tdepth == 0);$rate=$mdepth/$tdepth;
		$methy_seq=$a[-2];
	}
	if($ARGV[8] == 1){
		$methy_seq=&identify_methy_type($a[2],$a[3]);
	}
	next if($tdepth <$ARGV[13]);   ####next if($tdepth < $ARGV[3]);
	if($ARGV[7] ne "C"){
		next if($methy_seq ne $ARGV[7]);
	}
	if(exists $res{$a[0]}{$a[1]}){
		$methy2{$res{$a[0]}{$a[1]}}[0]+=$tdepth;$methy2{$res{$a[0]}{$a[1]}}[1]+=$mdepth;  ###将这个C的mC和tC储存在这个合并候选DMR中（methy2）
		if(exists $methy_info1{$a[0]}{$a[1]}){   #####这个C必须在两个样本都有覆盖
			next if($methy_info1{$a[0]}{$a[1]}[1] == 0 && $mdepth == 0);  ####这个C至少在一个样本中的mC不等于0
			$merge_result{$res{$a[0]}{$a[1]}}[2]++;   #####这个C所在的合并候选DMR的C的个数加1
			if($methy_info1{$a[0]}{$a[1]}[1] == 0 && $mdepth/$tdepth > $type{$methy_seq}){   #####
				$merge_result{$res{$a[0]}{$a[1]}}[0]++;  #####这个C所在的合并候选DMR的hyper DMC加1
			}
			if($mdepth == 0 && $methy_info1{$a[0]}{$a[1]}[1]/$methy_info1{$a[0]}{$a[1]}[0] > $type{$methy_seq}){
				$merge_result{$res{$a[0]}{$a[1]}}[1]++;   #####这个C所在的合并候选DMR的hypo DMC加1
			}
			if($mdepth != 0 && $methy_info1{$a[0]}{$a[1]}[1] != 0){
				my $npp = $methy_info1{$a[0]}{$a[1]}[0]+$tdepth; my $n1p =$methy_info1{$a[0]}{$a[1]}[0];
				my $np1 = $methy_info1{$a[0]}{$a[1]}[1]+$mdepth; my $n11 =$methy_info1{$a[0]}{$a[1]}[1];
				my $fisher_exact_test_twotail=calculateStatistic( n11=>$n11, n1p=>$n1p, np1=>$np1, npp=>$npp);
				if(defined $fisher_exact_test_twotail && $fisher_exact_test_twotail <= $ARGV[6]){
					my $res1=abs($methy_info1{$a[0]}{$a[1]}[1]/$methy_info1{$a[0]}{$a[1]}[0]-$mdepth/$tdepth);
					if($methy_seq eq "CG"){                                 next if($res1 < 0.4);           }
					if($methy_seq eq "CHG"){                                next if($res1 < 0.2);           }
					if($methy_seq eq "CHH"){                                next if($res1 < 0.1);           }
					if($methy_info1{$a[0]}{$a[1]}[1]/$methy_info1{$a[0]}{$a[1]}[0] < $mdepth/$tdepth){
						my $rate=$methy_info1{$a[0]}{$a[1]}[1]/$methy_info1{$a[0]}{$a[1]}[0];
						my $foldchange=($mdepth/$tdepth)/$rate;
						next if($foldchange < $ARGV[12]);
						$merge_result{$res{$a[0]}{$a[1]}}[0]++;  ####这个C所在的合并候选DMR的hyper DMC加1
					}
					else{
						my $rate=$mdepth/$tdepth;
						my $foldchange=($methy_info1{$a[0]}{$a[1]}[1]/$methy_info1{$a[0]}{$a[1]}[0])/$rate;
						next if($foldchange < $ARGV[12]);
						$merge_result{$res{$a[0]}{$a[1]}}[1]++;   ####这个C所在的合并候选DMR的hypo DMC加1
					}
				}
			}
		}
	}
}
close IN2;
print "Region\ttype\tDmC-hyper\tDmC-hypo\ttotal-C\t$ARGV[9]-methy\t$ARGV[10]-methy\tfoldchange\n";
foreach my $k(sort keys %merge_result){
	next if($merge_result{$k}[3] eq "hyper" && $merge_result{$k}[0] < $ARGV[1]);     #####只保留DMR类型（hyper或hypo）的DMC个数达到阈值的DMR
	next if($merge_result{$k}[3] eq "hypo"  && $merge_result{$k}[1] < $ARGV[1]);
	my ($rate1, $rate2)=(0, 0);	my $fc=1;
	if(exists $methy1{$k}){
		$rate1=$methy1{$k}[1]/$methy1{$k}[0];
	}
	if(exists $methy2{$k}){
		$rate2=$methy2{$k}[1]/$methy2{$k}[0];     #####计算合并候选得到DMR的FC
	}
	if($rate1 == 0 && $rate2 != 0){
		$fc=100*$rate2;   ####相当于让$rate1等于0.01，$fc=$rate2/0.01
	}
	if($rate1 != 0 && $rate2 == 0){
		$fc=100*$rate1;
	}
	if($rate1 != 0 && $rate2 != 0){   
		$fc=$rate1/$rate2;
		if($rate1 < $rate2){
			$fc=$rate2/$rate1;
		}
	}
	next if($fc <= $ARGV[11]);
	$fc=log($fc)/log(2);
	print "$k\t$merge_result{$k}[3]\t$merge_result{$k}[0]\t$merge_result{$k}[1]\t$merge_result{$k}[2]\t$rate1\t$rate2\t$fc\n";
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
