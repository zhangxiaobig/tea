#!/usr/bin/env python2.7

import sys

f_bed=open(sys.argv[1],'r')
mp_score={}
tlen={}
tlen_map={}
	
for line in f_bed.readlines():
	line=line.strip()
	items=line.split('\t')
	st=int(items[1])
	ed=int(items[2])
	start=int(items[7])
	end=int(items[8])
	score=float(items[9])
	name=items[3]
	
	tlen[name]=float(ed-st+1)
	if name not in mp_score:
		mp_score[name]=0
		tlen_map[name]=0
	
	if st>=start+1 and end <=ed:
		mp_score[name]+=score*(end-st+1)
		tlen_map[name]+=end-st+1
	elif ed<=end and start+1>=st:
		mp_score[name]+=score*(ed-start)
		tlen_map[name]+=ed-start
	elif start+1<=st and ed<=end:
		mp_score[name]+=score*(ed-st+1)
		tlen_map[name]+=ed-st+1
	else:
		mp_score[name]+=score*(end-start)
		tlen_map[name]+=end-start

for key in sorted(mp_score.keys()):
	print key+"\t"+str(mp_score[key]/tlen[key])

f_bed.close()	
	
	
