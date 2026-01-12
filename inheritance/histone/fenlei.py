#!/usr/bin/python
# -*-coding:utf-8 -*-
import math
import sys
import os
import collections
import re

def SnowMan(inputfile,outfile):
    target = open(inputfile,'r')
    result = open(outfile,'w')
    name = []
    parentel = ["p1","p2","p"]
    parents = []
    firstline = target.readline()
    q = firstline.rstrip().split('\t')
    Chr = q[0]
    start = int(q[1])
    end = int(q[2])
    Type = q[3]
    name.append(Type)
    if q[3] in parentel:
        parents.append(int(q[1]))
        parents.append(int(q[2]))
    for line in target:
#        print name
        a=line.rstrip().split('\t')
        if a[0] == Chr:
            if end >= int(a[1]) and int(a[2]) >= start:
                start = min(end,int(a[1]),int(a[2]),start)
                end = max(end,int(a[1]),int(a[2]),start)
                if a[3] in parentel:
                    parents.append(int(a[1]))
                    parents.append(int(a[2]))
#                    parents.append(int(a[2]))
                if a[3] in name:
                    continue    
                else:
                    name.append(a[3])
            else:
                if "p1" in name or "p2" in name or "p" in name:
                    kk1 = min(parents)
                    kk2 = max(parents)
                    result.write(Chr+'\t'+str(kk1)+'\t'+str(kk2)+'\t')
                else:
                    result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t')
                nameSort = sorted(name)
                for i in nameSort:
                    result.write(i)
                result.write('\n')
                start = int(a[1])
                end = int(a[2])
                name = []
                parents = []
                name.append(a[3])
                if a[3] in parentel:
                    parents.append(int(a[1]))
                    parents.append(int(a[2]))
#                    parents.append(int(a[2]))
        else:
            if "p1" in name or "p2" in name or "p" in name:
                kk1 = min(parents)
                kk2 = max(parents)
                result.write(Chr+'\t'+str(kk1)+'\t'+str(kk2)+'\t')
            else:
                result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t')
            nameSort = sorted(name)
            for i in nameSort:
                result.write(i)
            result.write('\n')
            start = int(a[1])
            end = int(a[2])
            Chr = a[0]
            name = []
            parents = []
            name.append(a[3])
            if a[3] in parentel:
                parents.append(int(a[1]))
                parents.append(int(a[2]))

    if "p1" in name or "p2" in name or "p" in name:
        kk1 = min(parents)
        kk2 = max(parents)
        result.write(Chr+'\t'+str(kk1)+'\t'+str(kk2)+'\t')
    else:
        result.write(Chr+'\t'+str(start)+'\t'+str(end)+'\t')
    nameSort = sorted(name)
    for i in nameSort:
        result.write(i)
    result.write('\n')
    result.close()

def main(argv):
    SnowMan(argv[1],argv[2])

if __name__=='__main__':
    main(sys.argv)
