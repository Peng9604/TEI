python /mnt/wangyw/other/xingdt/bivo_zhangpeng/make_fig/fig4-interaction/WT-ddm1-diff-type/WT-K4-venn/others-xifen.py WT-all-K27-sorted-venn-nocol4-others WT-K27-fHv_rHn WT-K27-fHv_rTCdM WT-K27-fHv_rTCM WT-K27-fTCM_rHn WT-K27-fTCM_rTCdM  WT-K27-fTCM_rHv WT-K27-fTCdM_rHn WT-K27-fTCdM_rHn WT-K27-fTCdM_rHv WT-K27-fHn_rTCdM WT-K27-fHn_rTCM WT-K27-fHn_rHv 2>log
wl WT-K27-* > others-xifen-tongji.xls
cat others-xifen-tongji.xls xifen-tongji.xls jiben-tongji.xls |grep -v "total" |sed '1,$s/^\s*//g' |sed '1,$s/ /\t/g'|sort -nr|perl -ne '($a,$b)=split;print "$b\t$a\n"' > all-type-tongji.xls
