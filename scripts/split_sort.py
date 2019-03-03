#!/usr/bin/env python2.7
# python split_sort.py Aligned.out.bed

import sys
import collections

f=open(sys.argv[1]+'/Aligned.out.bed','r')
ff=open(sys.argv[1]+'/line.txt','a')
split_line=collections.defaultdict(dict)
split_start=float("inf")
before_name=''
before_chr=''
line_split=''	

for line in f.readlines():
	line = line.strip()
	items = line.split('\t')
	if items[0]==before_chr or before_chr=='':
		if items[3]==before_name:
			line_split+=line+"\n"			
		else:
			print line
	else:
		print line 
		
	before_name=items[3]
	before_chr=items[0]
	
ff.write(line_split)
line_split=''
					
ff.close()
f.close()
	
		