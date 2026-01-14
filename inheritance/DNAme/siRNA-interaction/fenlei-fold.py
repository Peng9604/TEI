#!/usr/bin/python
# -*-coding:utf-8 -*-
####输入：前三列是siRNA mapping的区域坐标，常常有多个siRNA mapping到相同位置，所以是冗余的。后四列是siRNA mapping位置所属的大区域名称和坐标。
####功能：去除前三列的坐标冗余
import math
import sys
import os
import collections
import re

def SnowMan(inputfile,outfile):
    target = open(inputfile,'r')
    result = open(outfile,'w')
    firstline = target.readline()
    q = firstline.rstrip().split('\t')
    Chr = q[0]
    start = int(q[1])
    end = int(q[2])
    index = '\t'.join(q[3:])
    for line in target:
        a=line.rstrip().split('\t')
        if a[0] == Chr:
            if end >= int(a[1]) and int(a[2]) >= start:
                start = min(end,int(a[1]),int(a[2]),start)
                end = max(end,int(a[1]),int(a[2]),start)
                index = '\t'.join(a[3:])
            else:
                result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t'+index+'\n')
                start = int(a[1])
                end = int(a[2])
                index = '\t'.join(a[3:])
        else:
            result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t'+index+'\n')
            start = int(a[1])
            end = int(a[2])
            Chr = a[0]
            index = '\t'.join(a[3:])
            
    result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t'+index+'\n')
    result.close()

def main(argv):
    SnowMan(argv[1],argv[2])

if __name__=='__main__':
    main(sys.argv)
