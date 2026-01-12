grep K27 /mnt/wangyw/other/xingdt/bivo_zhangpeng/make_fig/fig1/upsetR-WT-ddm1-venn/sample-K4-K27-label-data/8sam-WT-all-venn/WT-all-K4K27.txt > WT-all-K27.txt
sort -k1,1 -k2,2n WT-all-K27.txt > WT-all-K27-sorted.txt
python /mnt/wangyw/soybean/K4K27-Bivalent/03/Bivalent/fenlei.py WT-all-K27-sorted.txt WT-all-K27-sorted-venn.txt
cut -f1-3 WT-all-K27-sorted-venn.txt > WT-all-K27-sorted-venn-nocol4.txt
