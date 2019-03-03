#!/usr/bin/env python2.7

import sys

f1=open(sys.argv[1],'r')
f2=open(sys.argv[2],'r')

mulit_count={}

for line in f1.readlines():
	line=line.strip()
	items=line.split("\t")
	issue1=items[0].rstrip()
	mulit_count[issue1]=int(items[1])
		
for line in f2.readlines():
	line=line.strip()
	items=line.split("\t")
	issue1=items[0].rstrip()
	mulit_count[issue1]+=int(items[1])
	
for key in mulit_count:
	print key+"\t"+str(mulit_count[key])

f1.close()
f2.close()	