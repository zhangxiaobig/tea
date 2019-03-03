#!/usr/bin/env python2.7
# python TE_transcript_grow.py RMSK_uniq_reads_overlaped_filtered_merge.bed

import sys

f1=open(sys.argv[1],'r')
chr={}
strand={}
start={}
end={}
uniq_count={}


for line in f1.readlines():
	line = line.strip()
	items = line.split('\t')
	name = items[7]
	
	if name not in chr:
		chr[name]=items[0]
		strand[name]=items[9]
		
		
	if name not in start:
		start[name]=int(items[1])
		end[name]=int(items[2])
		uniq_count[name]=int(items[3])
	else:
		uniq_count[name]+=int(items[3])
		if items[1]<start[name]:
			start[name]=int(items[1])
		if items[2]>end[name]:
			end[name]=int(items[2])
			
for key in chr:
	print chr[key]+"\t"+str(start[key])+"\t"+str(end[key])+"\t"+key+"\t"+strand[key]+"\t"+str(uniq_count[key])
		
f1.close()