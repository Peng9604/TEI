import pysam
import sys
import collections

def readsnp(inputf):
    snp = collections.OrderedDict()
    with open(inputf,'r') as f:
        for line in f:
            a=line.rstrip().split('\t')
            if a[0] in snp:
                snp[a[0]][a[1]]=(a[2],a[3])
            else:
                snp[a[0]]=collections.OrderedDict()
                snp[a[0]][a[1]]=(a[2],a[3])
    return snp

def Import(inputfile,snp,outf1,outf2):
    P1 = open(outf1,'w')
    P2 = open(outf2,'w')
    bf = pysam.AlignmentFile(inputfile, 'rb')
    for key1 in snp:
    #    print key1
        for key2 in snp[key1]:
  #          print key1,key2
            #A = '\''+str(key1)+'\''
            for pileupcolumn in bf.fetch(key1,int(key2)-1,int(key2)):
#                print pileupcolumn
#                print type(pileupcolumn)
                W = str(pileupcolumn)
                w = W.rstrip('\n').split('\t')
#                print w[9]
#                print W
#                print key2,w[3],int(key2)-int(w[3])-1,snp[key1][key2][0],snp[key1][key2][1]
#                print w[9][int(key2)-int(w[3])-1]
                try:
                    if w[9][int(key2)-int(w[3])-1] == snp[key1][key2][0]:
                        P1.write(W+'\n')
                    if w[9][int(key2)-int(w[3])-1] == snp[key1][key2][1]:
                        P2.write(W+'\n')
                except IndexError:
                    pass
                continue
    P1.close()
    P2.close()
    bf.close()

def main(argv):
    snp = readsnp(argv[1])
    Import(argv[2],snp,argv[3],argv[4])
    
          
if __name__=='__main__':
    main(sys.argv)

