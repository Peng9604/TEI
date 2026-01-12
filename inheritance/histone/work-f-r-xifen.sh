fTCM_K27=$1
rTCM=$2
fTCdM_K27=$3
rTCdM_K27=$4

#fHv
awk '{if($7>=$8&&$9>$7)print$0}{if($8>=$7&&$9>$8)print$0}' $1 >${1}-fHv
#fTCM
awk '{if($7>=$8&&$9<$7)print$0}{if($8>=$7&&$9<$8)print$0}' $1 >${1}-fTCM
#fHn
awk '{if($7>=$8&&$9<$8)print$0}{if($8>=$7&&$9<$7)print$0}' $3>${3}-fHn
#fTCdM
awk '{if($7>=$8&&$9>$8)print$0}{if($8>=$7&&$9>$7)print$0}' $3>${3}-fTCdM

#rHv
awk '{if($7>=$8&&$10>$7)print$0}{if($8>=$7&&$10>$8)print$0}' $2 >${2}-rHv
#rTCM
awk '{if($7>=$8&&$10<$7)print$0}{if($8>=$7&&$10<$8)print$0}' $2 >${2}-rTCM
#rHn
awk '{if($7>=$8&&$10<$8)print$0}{if($8>=$7&&$10<$7)print$0}' $4>${4}-rHn
#rTCdM
awk '{if($7>=$8&&$10>$8)print$0}{if($8>=$7&&$10>$7)print$0}' $4>${4}-rTCdM
