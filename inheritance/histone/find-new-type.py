#!/usr/bin/python
import sys 
import collections

def get_peak(peak_file):
    peak_dic=collections.OrderedDict()
    with open(peak_file,'r') as f:
        for line in f:
            a=line.rstrip().split('\t')
            if a[0] in peak_dic:
                peak_dic[a[0]][a[1]]=a[2]
            else:
                peak_dic[a[0]]=collections.OrderedDict()
                peak_dic[a[0]][a[1]]=a[2]
    return peak_dic
def get_FDR(inputfile,peak_dic):
    FDR_dic=collections.OrderedDict()
    with open(inputfile,'r') as f:
        head=f.readline()
        for line in f:
            a=line.rstrip().split('\t')
            if a[1] in peak_dic[a[0]]:
                if float(a[9])!=1:   #######应改成：if float(a[7])>1:   用FC是否大于1来判断A与B谁大
                    if a[0] in FDR_dic:
                        FDR_dic[a[0]][a[1]]=(a[2],'+'+str(a[9]),a[4],a[6],a[7])
                    else:
                        FDR_dic[a[0]]=collections.OrderedDict()
                        FDR_dic[a[0]][a[1]]=(a[2],'+'+str(a[9]),a[4],a[6],a[7])
                else:
                    if a[0] in FDR_dic:
                        FDR_dic[a[0]][a[1]]=(a[2],'-'+str(a[12]),a[4],a[6],a[10])
                    else:
                        FDR_dic[a[0]]=collections.OrderedDict()
                        FDR_dic[a[0]][a[1]]=(a[2],'-'+str(a[12]),a[4],a[6],a[10])
    return FDR_dic

def merge_FDR(parent_FDR,F1_FDR,rF1_FDR,outfile):
    ot=open(outfile,'w')
    merge_dic=collections.OrderedDict()
    ot.write('#chr'+'\t'+'star'+'\t'+'end'+'\t'+'FDR_LUC_vs_COL'+'\t'+'FDR_CXL_vs_pamix'+'\t'+'FDR_LXC_vs_pamix'+'\t'+'LUC'+'\t'+'COL'+'\t'+'CXL'+'\t'+'LXC'+'\t'+'PAV'+'\t'+'FC_LUC_vs_COL'+'\t'+'FC_CXL_vs_PAV'+'\t'+'FC_LXC_vs_PAV'+'\n')
    for key1 in parent_FDR:
        for key2 in parent_FDR[key1]:
             ot.write(key1+'\t'+key2+'\t'+parent_FDR[key1][key2][0]+'\t'+parent_FDR[key1][key2][1]+'\t'+F1_FDR[key1][key2][1]+'\t'+rF1_FDR[key1][key2][1]+'\t'+parent_FDR[key1][key2][2]+'\t'+parent_FDR[key1][key2][3]+'\t'+F1_FDR[key1][key2][2]+'\t'+rF1_FDR[key1][key2][2]+'\t'+rF1_FDR[key1][key2][3]+'\t'+parent_FDR[key1][key2][4]+'\t'+F1_FDR[key1][key2][4]+'\t'+rF1_FDR[key1][key2][4]+'\n')
             if key1 in merge_dic:
                 merge_dic[key1][key2]=(parent_FDR[key1][key2][0],parent_FDR[key1][key2][1],F1_FDR[key1][key2][1],rF1_FDR[key1][key2][1],parent_FDR[key1][key2][2],parent_FDR[key1][key2][3],F1_FDR[key1][key2][2],rF1_FDR[key1][key2][2],rF1_FDR[key1][key2][3])
             else:
                 merge_dic[key1]=collections.OrderedDict()
                 merge_dic[key1][key2]=(parent_FDR[key1][key2][0],parent_FDR[key1][key2][1],F1_FDR[key1][key2][1],rF1_FDR[key1][key2][1],parent_FDR[key1][key2][2],parent_FDR[key1][key2][3],F1_FDR[key1][key2][2],rF1_FDR[key1][key2][2],rF1_FDR[key1][key2][3])
    return merge_dic

def diff_type(merge_dic,NI,up,down,Hv,TCM,Hn,TCdM,others,fTCM,rTCM,fTCdM,rTCdM,FDR):
    NI = open(NI,'w')
    up = open(up,'w')
    down = open(down,'w')
    Hv = open(Hv,'w')
    TCM = open(TCM,'w')
    Hn = open(Hn,'w')
    TCdM = open(TCdM,'w')
    others = open(others,'w')
    fTCM = open(fTCM,'w')
    rTCM = open(rTCM,'w')
    fTCdM = open(fTCdM,'w')
    rTCdM = open(rTCdM,'w')
    High=0
    Low=0
    for key1 in merge_dic:
        for key2 in merge_dic[key1]:
            if abs(float(merge_dic[key1][key2][2]))>float(FDR) and abs(float(merge_dic[key1][key2][3]))>float(FDR):
                NI.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif merge_dic[key1][key2][2].startswith('+') and merge_dic[key1][key2][3].startswith('+'):
                up.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                High=max(float(merge_dic[key1][key2][4]),float(merge_dic[key1][key2][5]))
                if float(merge_dic[key1][key2][2])>float(FDR) and float(merge_dic[key1][key2][3])<float(FDR):
                    rTCM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                elif float(merge_dic[key1][key2][2])<float(FDR) and float(merge_dic[key1][key2][3])>float(FDR):
                    fTCM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                elif float(merge_dic[key1][key2][6])>float(High) and float(merge_dic[key1][key2][7])>float(High):
                    Hv.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                else:
                    TCM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif merge_dic[key1][key2][2].startswith('-') and merge_dic[key1][key2][3].startswith('-'):
                down.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                Low=min(float(merge_dic[key1][key2][4]),float(merge_dic[key1][key2][5]))
                if float(merge_dic[key1][key2][2])>-float(FDR) and float(merge_dic[key1][key2][3])<-float(FDR):
                    fTCdM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                elif float(merge_dic[key1][key2][2])<-float(FDR) and float(merge_dic[key1][key2][3])>-float(FDR):
                    rTCdM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                elif float(merge_dic[key1][key2][6])<float(Low) and float(merge_dic[key1][key2][7])<float(Low):
                    Hn.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
                else:
                    TCdM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif float(merge_dic[key1][key2][2])<float(FDR) and float(merge_dic[key1][key2][3])<-float(FDR):
                    fTCM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif float(merge_dic[key1][key2][2])>-float(FDR) and float(merge_dic[key1][key2][3])>float(FDR):
                    fTCdM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif float(merge_dic[key1][key2][2])<-float(FDR) and float(merge_dic[key1][key2][3])<float(FDR):
                    rTCM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            elif float(merge_dic[key1][key2][2])>float(FDR) and float(merge_dic[key1][key2][3])>-float(FDR):
                    rTCdM.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            # elif float(merge_dic[key1][key2][2])==0 or float(merge_dic[key1][key2][3])==0:
            #     if merge_dic[key1][key2][6].startswith('+') and merge_dic[key1][key2][7].startswith('+'):
            #         Hv.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            #     else:
            #         Hn.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
            else:
                others.write(key1+'\t'+key2+'\t'+merge_dic[key1][key2][0]+'\t'+merge_dic[key1][key2][1]+'\t'+merge_dic[key1][key2][2]+'\t'+merge_dic[key1][key2][3]+'\t'+merge_dic[key1][key2][4]+'\t'+merge_dic[key1][key2][5]+'\t'+merge_dic[key1][key2][6]+'\t'+merge_dic[key1][key2][7]+'\t'+merge_dic[key1][key2][8]+'\n')
    NI.close()
    up.close()
    down.close()
    Hv.close()
    TCM.close()
    Hn.close()
    TCdM.close()
    others.close()
    fTCM.close()
    rTCM.close()
    fTCdM.close()
    rTCdM.close()

    
def main(argv):
    peak=get_peak(argv[1])
    parents=get_FDR(argv[2],peak)
    F1=get_FDR(argv[3],peak)
    rF1=get_FDR(argv[4],peak)
    all_FDR=merge_FDR(parents,F1,rF1,argv[5])
    diff_type(all_FDR,argv[6],argv[7],argv[8],argv[9],argv[10],argv[11],argv[12],argv[13],argv[14],argv[15],argv[16],argv[17],argv[18])

# def diff_type(merge_dic,NI,up,down,Hv,TCM,Hn,TCdM,others,FDR):
if __name__=='__main__':
    main(sys.argv)

