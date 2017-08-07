#! /usr/bin/env python
from collections import namedtuple
def intlist(k):
	k=k.split(",")
	li=[int(n) for n in k if k!=""]
	return li
def load_table(lines,cast="n",header="y"):
	lines=iter(lines)
	
	if header=="y":
		 
		hl = lines.next()
		hl=hl.replace("#","")
		if cast=="y":
			types=[n.strip().split(":")[1] for n in hl.split("\t")]
		lt=[n.strip().split(":")[0] for n in hl.split("\t")]
		Table=namedtuple('Table',lt)

		
	for l in lines:
		l=l.rstrip()
		lt=l.split("\t")
		l1=lt
		tbl=Table(*l1)
		yield tbl

if __name__=="__main__":
	
	from optparse import OptionParser
	parser=OptionParser()
	parser.add_option("-c","--cast",dest="cast",default="n")
	(opts,args)=parser.parse_args()
	import sys
	f=open(args[0])
	lines=f.readlines()
	for t in load_table(lines,opts.cast):
		print t.name
