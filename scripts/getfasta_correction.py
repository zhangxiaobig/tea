#!/usr/bin/env python2.7

import sys

f=open(sys.argv[1],'r')

for line in f.readlines():
	line = line.strip()
	items = line.split(':')
	if len(items)>0:
		print items[0]
	else:
		print line

f.close()