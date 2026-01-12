perl /media/newdisks/HDD1/xingdt/bin/Methylation/DMR_use_parents_mc_to_f1_summary-tair.pl /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC-new/COL_LUC_CG.dep4.CG.win200.dmc7.gap100 /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/COL_methy.txt /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/LUC_methy.txt /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/LXC_methy.txt CG 4 LXC_COL_LUC-CG-MI.interaction.out
perl /media/newdisks/HDD1/xingdt/bin/Methylation/MI_split.pl LXC_COL_LUC-CG-MI.interaction.out LXC_COL_LUC-CG

perl /media/newdisks/HDD1/xingdt/bin/Methylation/DMR_use_parents_mc_to_f1_summary-tair.pl /media/newdisks/HDD6/xingdt/AT_meth/01.dmr/COL_LUC-new/COL_LUC_CG.dep4.CG.win200.dmc7.gap100 /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/COL_methy.txt /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/LUC_methy.txt /media/newdisks/HDD6/xingdt/AT_meth/methytable_noMtPt/CXL_methy.txt CG 4 CXL_COL_LUC-CG-MI.interaction.out
perl /media/newdisks/HDD1/xingdt/bin/Methylation/MI_split.pl CXL_COL_LUC-CG-MI.interaction.out CXL_COL_LUC-CG



ls *.out|while read dd;do sed -i '1,$s/:/\t/g' $dd;done
ls *.out|while read dd;do sed -i '1,$s/\-/\t/g' $dd;done
ls *.out|while read dd;do a=`echo $dd|cut -d '_' -f1`;b=`echo $dd|cut -d '.' -f2`;awk -v z=$a -v OFS="\t" -v k=$b '{print $1,$2,$3,z"-"k}' $dd > $dd-1;done
cat CXL_COL_LUC-CG.NI.out-1 LXC_COL_LUC-CG.NI.out-1 > CXL-LXC-CG.NI.out
cat CXL_COL_LUC-CG.TCM.out-1 LXC_COL_LUC-CG.TCM.out-1 > CXL-LXC-CG.TCM.out
cat CXL_COL_LUC-CG.TCdM.out-1 LXC_COL_LUC-CG.TCdM.out-1 > CXL-LXC-CG.TCdM.out
sort -k1,1 -k2,2n CXL-LXC-CG.NI.out > CXL-LXC-CG.NI.out-sorted
sort -k1,1 -k2,2n CXL-LXC-CG.TCM.out > CXL-LXC-CG.TCM.out-sorted
sort -k1,1 -k2,2n CXL-LXC-CG.TCdM.out > CXL-LXC-CG.TCdM.out-sorted

python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py CXL-LXC-CG.NI.out-sorted CXL-LXC-CG.NI.out-sorted-venn
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py CXL-LXC-CG.TCM.out-sorted CXL-LXC-CG.TCM.out-sorted-venn
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py CXL-LXC-CG.TCdM.out-sorted CXL-LXC-CG.TCdM.out-sorted-venn
grep 'CXL.*LXC\|LXC.*CXL' CXL-LXC-CG.TCdM.out-sorted-venn > CXL-LXC-CG.TCdM.out-sorted-venn-last.txt
grep 'CXL.*LXC\|LXC.*CXL' CXL-LXC-CG.TCM.out-sorted-venn > CXL-LXC-CG.TCM.out-sorted-venn-last.txt
grep 'CXL.*LXC\|LXC.*CXL' CXL-LXC-CG.NI.out-sorted-venn > CXL-LXC-CG.NI.out-sorted-venn-last.txt

