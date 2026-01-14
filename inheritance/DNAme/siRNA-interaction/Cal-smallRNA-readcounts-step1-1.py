#coding:utf-8
import pysam
import sys
#usage:python Cal-smallRNA-readcounts-step1-1.py sortbam.list region.bed

dic = {}
f2 = open(sys.argv[2],'r')
for line in f2:
    a = line.rstrip().split("\t")
    dic[a[0]]=[a[1],a[2],a[3]]
f2.close()

f2name=sys.argv[2].rstrip().split(".")

f1 = open(sys.argv[1],'r')
for line in f1.readlines():
    k = line.rstrip().split("\t")
    bamname = k[1]
    bam = pysam.AlignmentFile(bamname,"rb")
    with open(k[0]+'-'+f2name[0]+'.readcounts.txt','w') as OUT:
        for name in dic.keys():
            chrom=dic[name][0]
            start=int(dic[name][1])
            end=int(dic[name][2])
            reads = bam.fetch(contig=chrom, start=start, end=end)
            for line in reads:
                line1=str(line)
                l = line1.split("\t")
                OUT.write(name+"\t"+dic[name][0]+"\t"+dic[name][1]+"\t"+dic[name][2]+"\t"+l[0]+"\t"+str(line.reference_name)+"\t"+str(line.reference_start)+"\t"+str(line.query_alignment_length)+"\t"+l[9]+"\t"+str(line.reference_end)+"\n")   ####输出格式：TE或基因或bin的名称，region-chr，region-start，region-end，smallRNA的名称，smallRNA比对上的chr，smallRNA比对上的start，smallRNA比对上基因组的长度，smallRNA序列
f1.close()
