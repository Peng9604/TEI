awk -v OFS='\t' '{if($2-1000<0){start=0}else{start=$2-1000};print $1,start,$3+1000}' CXL-CG-TCMLXC-CG-TCM.bed > CXL-LXC-CG-SMR-TCM_1k.txt
awk '{print "CG-SMR-TCM_"NR"\t"$1"\t"$2"\t"$3}' CXL-LXC-CG-SMR-TCM_1k.txt > CG-SMR-TCM.bed
ls /media/newdisks/HDD6/xingdt/nrpenrpd-WT_smallRNA/bowtie-mapped-nolimit-3mis/*.sort.bam|perl -e 'while(<>){chomp;@a=split/\//;@b=split/\./,$a[-1];print"$b[0]\t$_\n";}' > sortbam.list
python /media/newdisks/HDD1/xingdt/bin/smallRNA/Cal-smallRNA-readcounts-step1-1.py sortbam.list CG-SMR-TCM.bed
###C24XCol
awk -v OFS='\t' '{if($8 == 24){print $6,$7,$10,$1,$2,$3,$4}}' C24XCol_no-tRNA_rRNA_snRNA_tair-CG-SMR-TCM.readcounts.txt > C24XCol-CG-SMR-TCM.readcounts.txt
 sort -k1,1 -k2,2n C24XCol-CG-SMR-TCM.readcounts.txt > C24XCol-CG-SMR-TCM.readcounts-sorted.txt
python /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/fenlei-fold.py C24XCol-CG-SMR-TCM.readcounts-sorted.txt C24XCol-CG-SMR-TCM.readcounts-sorted-fold.txt
###ColXC24
awk -v OFS='\t' '{if($8 == 24){print $6,$7,$10,$1,$2,$3,$4}}' ColXC24_no-tRNA_rRNA_snRNA_tair-CG-SMR-TCM.readcounts.txt > ColXC24-CG-SMR-TCM.readcounts.txt
sort -k1,1 -k2,2n ColXC24-CG-SMR-TCM.readcounts.txt > ColXC24-CG-SMR-TCM.readcounts-sorted.txt
python /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/fenlei-fold.py ColXC24-CG-SMR-TCM.readcounts-sorted.txt ColXC24-CG-SMR-TCM.readcounts-sorted-fold.txt
###Col
awk -v OFS='\t' '{if($8 == 24){print $6,$7,$10,$1,$2,$3,$4}}' Col_no-tRNA_rRNA_snRNA_tair-CG-SMR-TCM.readcounts.txt > Col-CG-SMR-TCM.readcounts.txt
sort -k1,1 -k2,2n Col-CG-SMR-TCM.readcounts.txt > Col-CG-SMR-TCM.readcounts-sorted.txt
python /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/fenlei-fold.py Col-CG-SMR-TCM.readcounts-sorted.txt Col-CG-SMR-TCM.readcounts-sorted-fold.txt
###C24
awk -v OFS='\t' '{if($8 == 24){print $6,$7,$10,$1,$2,$3,$4}}' C24_no-tRNA_rRNA_snRNA_tair-CG-SMR-TCM.readcounts.txt > C24-CG-SMR-TCM.readcounts.txt
 sort -k1,1 -k2,2n C24-CG-SMR-TCM.readcounts.txt > C24-CG-SMR-TCM.readcounts-sorted.txt
python /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/fenlei-fold.py C24-CG-SMR-TCM.readcounts-sorted.txt C24-CG-SMR-TCM.readcounts-sorted-fold.txt

awk -v OFS='\t' '{print $1,$2,$3,"Col-"$4}' Col-CG-SMR-TCM.readcounts-sorted-fold.txt > Col-CG-SMR-TCM.readcounts-sorted-fold-forC24.txt
awk -v OFS='\t' '{print $1,$2,$3,"C24-"$4}' C24-CG-SMR-TCM.readcounts-sorted-fold.txt > C24-CG-SMR-TCM.readcounts-sorted-fold-forCol.txt
cat Col-CG-SMR-TCM.readcounts-sorted-fold-forC24.txt C24-CG-SMR-TCM.readcounts-sorted-fold-forCol.txt |sort -k1,1 -k2,2n > Col-C24-CG-SMR-TCM-overlap-sort.txt
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py Col-C24-CG-SMR-TCM-overlap-sort.txt Col-C24-CG-SMR-TCM-overlap-sort-venn.txt
grep "Col" Col-C24-CG-SMR-TCM-overlap-sort-venn.txt|grep -v "C24"|sed -e '1,$s/Col-/Col\t/g' > Col_no-C24-CG-SMR-TCM.txt
grep "C24" Col-C24-CG-SMR-TCM-overlap-sort-venn.txt|grep -v "Col"|sed -e '1,$s/C24-/C24\t/g' > C24_no-Col-CG-SMR-TCM.txt
perl /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/2sample-uniq.pl Col_no-C24-CG-SMR-TCM.txt C24_no-Col-CG-SMR-TCM.txt Col C24 Col-C24-CG-SMR-TCM

awk -v OFS='\t' '{print $1,$2,$3,"ColXC24-"$4}' ColXC24-CG-SMR-TCM.readcounts-sorted-fold.txt > ColXC24-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt
awk -v OFS='\t' '{print $1,$2,$3,"C24XCol-"$4}' C24XCol-CG-SMR-TCM.readcounts-sorted-fold.txt > C24XCol-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt
 cat ColXC24-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt C24XCol-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt Col-C24-CG-SMR-TCM-Col-uniq.bed5 |sort -k1,1 -k2,2n > CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-sorted.txt
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-sorted.txt CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-sorted-venn.txt
grep -E 'C24XCol.*Col.*ColXC24' CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-sorted-venn.txt > CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-last.txt

cat ColXC24-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt C24XCol-CG-SMR-TCM.readcounts-sorted-fold-forparents.txt Col-C24-CG-SMR-TCM-C24-uniq.bed5 |sort -k1,1 -k2,2n > CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-sorted.txt
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-sorted.txt CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-sorted-venn.txt
grep -E 'C24.*C24XCol.*ColXC24' CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-sorted-venn.txt > CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-last.txt

perl /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC/new_DMR_SMR/PAV/TCM-siRNA_interaction/parent-uniq_F1-overlap.pl CG-SMR-TCM-Col-uniq_ColXC24_C24XCol-last.txt CG-SMR-TCM-C24-uniq_ColXC24_C24XCol-last.txt Col C24 CG-SMR-TCM CG-SMR-TCM-parent-uniq_F1-overlap

