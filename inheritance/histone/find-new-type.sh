#bash/bin

bed=$1
w_j=$2
wj_pa=$3
jw_pa=$4
python /mnt/wangyw/other/xingdt/bivo_zhangpeng/make_fig/fig4-interaction/WT-ddm1-diff-type/WT-K27-venn/find-new-type.py $bed ${w_j} ${wj_pa} ${jw_pa} P-F1-FDR-signal ${bed%.*}-NI ${bed%.*}-up ${bed%.*}-down ${bed%.*}-Hv ${bed%.*}-TCM ${bed%.*}-Hn ${bed%.*}-TCdM ${bed%.*}-others ${bed%.*}-fTCM ${bed%.*}-rTCM ${bed%.*}-fTCdM ${bed%.*}-rTCdM 0.01
