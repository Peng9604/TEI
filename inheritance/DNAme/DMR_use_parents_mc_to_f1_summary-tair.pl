#!/usr/local/bin/perl -w
use strict;
use Statistics::Descriptive;
use Statistics::Distributions;
use Text::NSP::Measures::2D::Fisher::twotailed;
die"usage:perl $0 <dmr> <S1_meth_table> <S2_meth_table> <F1_meth_table> <meth_type：C/CG/CHG/CHH> <min_read，一般是4> <OUT>\n"unless @ARGV == 7;
######功能：统计两个亲本在两个亲本DMR上的mC总数和tC总数(命名为mp和tp)，统计一种杂交（正交或反交）在DMR上的mC总数和tC总数(命名为mf和tf)，过滤掉F1少于7个合格的C（合格的C指F1中这个C上有大于4个reads覆盖度）的这个dmr
##########接下来对于合格的DMR，F1的mC和tC与两个亲本总mC和总tC做fishier检验，及F1与PAV（mp与tp之比）比较是否具有显著差异，保留FDR输出，该脚本不会过滤FDR。对F1高于PAV、低于PAV的DMR命名为H和L，也一并输出。
####输入：
#1. DMR文件，三部曲直出
#2. 第二、三、四文件依次是两个亲本和子代的methytable
#3. 统计的甲基化类型：C/CG/CHG/CHH
#4. 在加载methytable时过滤掉少于n个reads覆盖的C

####输出文件（最后一个参数）：LXC_COL_LUC-CHG-MI.interaction.out
####格式：dmr的chr、start、end、dmr的chr，dmr在这个染色质中排第几个，PAV比F1高：L、PAV比F1低：H，两者相等：C)，hyper/hypo、PAV的甲基化水平、F1的甲基化水平、p值、FDR（fishier检验得到的）


open IN1,$ARGV[0];open IN2,$ARGV[1];open IN3,$ARGV[2];open IN4,$ARGV[3];
open OA,">$ARGV[6]";
my %tr=(
    'CGA'=>"CG",'CGT'=>"CG",'CGC'=>"CG",'CGG'=>"CG",'CAG'=>"CHG",'CTG'=>"CHG",'CCG'=>"CHG"
);

my (%dmr_info,%order,%info,%meth1,%meth2,%meth3);

while(<IN2>){
    chomp;my @a=split/\t/;next if $a[0]=~/chr$/;next if $a[5] < $ARGV[5];
    if($ARGV[4] eq 'C'){
        $meth1{$a[0]}{$a[1]}=[$a[5],$a[6]];    #将亲本methytable中的染色质和位置作为key，这个位置的totalC和mC为value
    }
    else{
        my $type=&identify_methy_type($a[2],$a[3]);
        next if $type ne $ARGV[4];
        $meth1{$a[0]}{$a[1]}=[$a[5],$a[6]];
    }
}

close IN2;

while(<IN3>){
    chomp;my @a=split/\t/;next if $a[0]=~/chr$/;next if $a[5] < $ARGV[5];
    if($ARGV[4] eq 'C'){
        $meth2{$a[0]}{$a[1]}=[$a[5],$a[6]];
    }
    else{
        my $type=&identify_methy_type($a[2],$a[3]);
        next if $type ne $ARGV[4];
        $meth2{$a[0]}{$a[1]}=[$a[5],$a[6]];
    }
}

close IN3;

while(<IN4>){
    chomp;my @a=split/\t/;next if $a[0]=~/chr$/;next if $a[5] < $ARGV[5];
    if($ARGV[4] eq 'C'){
        $meth3{$a[0]}{$a[1]}=[$a[5],$a[6]];
    }
    else{
        my $type=&identify_methy_type($a[2],$a[3]);
        next if $type ne $ARGV[4];
        $meth3{$a[0]}{$a[1]}=[$a[5],$a[6]];
    }
}

close IN4;

my $n = 0;
my @chr = ();

while(<IN1>){
    chomp;my @a=split/\t/;next if $a[0]=~/^R/;
    my @b=split/[:-]/,$a[0];

    ####每个染色质的DMR按照顺序依次对应一个n值。每次换到新的chr，$n就重新变为1
	if(!@chr){
		push @chr,$b[0];
		$n++;
	}else{
		if($chr[-1] eq "$b[0]"){
			$n++;
		}else{
			push @chr,$b[0];
			$n = 1;
		}
	}

	my $count_C_in_f1 = 0;

	my ($mp,$tp,$mf,$tf)=(0,0,0,0);
    foreach my $k1($b[1]..$b[2]){
        if(exists $meth1{$b[0]}{$k1}){
            $mp+=$meth1{$b[0]}{$k1}[1];$tp+=$meth1{$b[0]}{$k1}[0];  #$mp是父母本在这个DMR中的甲基化C数量之和，$tp是父母本在这个DMR中的所有C的数量之和
        }
        if(exists $meth2{$b[0]}{$k1}){
            $mp+=$meth2{$b[0]}{$k1}[1];$tp+=$meth2{$b[0]}{$k1}[0];
        }
        if(exists $meth3{$b[0]}{$k1}){
            $mf+=$meth3{$b[0]}{$k1}[1];$tf+=$meth3{$b[0]}{$k1}[0];   #$mf是杂交在这个DMR中的甲基化C数量之和，$tf是杂交在这个DMR中的所有C的数量之和
			$count_C_in_f1++;
        }
    }
    
	next if $tf == 0; 
	next if $count_C_in_f1 < 7; ####在这个dmr上，F1合格的C的个数是$count_C_in_f1。过滤掉F1少于7个合格的C的这个dmr

	my $npp = $tp + $tf; my $n1p = $tf;
    my $np1 = $mp + $mf; my $n11 = $mf;
    my $fisher_exact_test_twotail=calculateStatistic( n11=>$n11, n1p=>$n1p, np1=>$np1, npp=>$npp);
    if(!$fisher_exact_test_twotail){ $fisher_exact_test_twotail = 1;}

	$dmr_info{$b[0]}{$n}=$fisher_exact_test_twotail;
    my $PAV = $mp/$tp; my $F1 = $mf/$tf;    #PAV算法：一个dmr上两个亲本的总甲基化C的数量除以tC的数量
    my $id = "";
    if( $F1 > $PAV ){
        $id = "H";
    }
    elsif ( $F1 < $PAV ){
        $id = "L";
    }
    else{
        $id = "C";
    }
    $order{$b[0]}{$n}=$a[0];
    $info{$b[0]}{$n}=[$a[1],$id,$PAV,$F1,$fisher_exact_test_twotail];
}

my @out = ("chr1","chr2","chr3","chr4","chr5");

foreach my $k1(@out){
    my %Re_BH = &BH_correct( \%{$dmr_info{$k1}} );
    foreach my $k2(sort {$a <=> $b} keys %{$dmr_info{$k1}}){
    	print OA "$order{$k1}{$k2}\t$k1\t$k2";
    	foreach my $k3(@{$info{$k1}{$k2}}){
        	print OA "\t$k3";
    	}
    	print OA "\t$Re_BH{$k2}\n";   #输出依次是：dmr的位置信息、$id(PAV比F1高：L、PAV比F1低：H，两者相等：C)、hyper/hypo、PAV的甲基化水平、F1的甲基化水平、p值、FDR
    }
}

sub identify_methy_type{
    my $direction=shift;
    my $file=shift;
    my ($seq, $methy_seq);
    if($direction eq "-"){
        $seq=substr($file,0,3);        $seq=reverse($seq);        $seq=~tr/ATCG/TAGC/;
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
