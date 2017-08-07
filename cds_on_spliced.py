#!/usr/bin/env python
# coding: utf-8
from table_loader import load_table
import sys
import numpy
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from man_seq import *
from collections import defaultdict
from collections import namedtuple
from optparse import OptionParser
from scipy import stats
from sequence_data import hg19
from aux_py.genes.gene_models import *
def cds(a_lines,g_lines,lines,system_path,opts): #use as input gene models. Return name and start end of CDS on spliced sequence
	circ2ref=defaultdict(str)
	g_frame=defaultdict(list)
	g_exonStarts=defaultdict(list)
	g_exonEnds=defaultdict(list)
	g_cons=defaultdict(numpy.array)

	for a in a_lines:
		a=a.rstrip().split("\t")

		if len(a)>23:
			ref=a[26]
			ref=ref.replace("best=","").replace('"',"").replace(";","")
			circ2ref[a[3]]=ref
	
	for l in get_transcripts(g_lines,system_path):
		g_frame[l.name]=l.exonFrames
		g_exonStarts[l.name]= l.exonStarts
		g_exonEnds[l.name]= l.exonEnds

	for l in get_transcripts(lines,system_path):
		i=0

		if circ2ref[l.name] and circ2ref[l.name] in g_frame:				
			name=circ2ref[l.name]
			cons=numpy.concatenate(map(lambda s,e:hg19.placentalmammals_conservation.get_oriented(l.chrom,s,e,l.strand),g_exonStarts[name],g_exonEnds[name]))
			for s,e in zip(g_exonStarts[name],g_exonEnds[name]):

				if s<=l.cdsStart<e:
					frame=g_frame[name][i]
					break

				i+=1	 	
			cdsC=l.cds_on_spliced
#decide whether the name is an id or coordinate

			if opts=="id":
				name=l.name

			elif opts=="coor":

				name="%s(%s):%s-%s"%(l.chrom,l.strand,l.txStart,l.txEnd)
			yield (name,cdsC,frame,cons)


def third_codon(name,cdsC,frame,cons):
		pos3=list()
		pos12=list()
		pos3=[n for i,n in enumerate(cons[cdsC.start:cdsC.end]) if (i+frame)%3==0]
		pos12=[n for i,n in enumerate(cons[cdsC.start:cdsC.end]) if (i+frame)%3!=0]
		if len(pos3)>0 and len(pos12)>0:
			statistic = '%6.3f\t%6.4f' %  stats.ttest_ind(pos3,pos12)
		else: statistic = 'NA\tNA' 
		print "%s\t%s" %(name,statistic)

	#print n, stats.ttest_ind(pos3,pos12)

	

if __name__=="__main__":
	parser=OptionParser()
	parse=parser.add_option("-n","--name-format",dest="name_format",default="id")
	parse=parser.add_option("-g","--gene_models",dest="gene_models",default="/data/BIO2/annae/uscs/hg19_ucsc_genes.csv")
	parse=parser.add_option("-a","--annotation",dest="annotation",default="/data/BIO2/sponge/code/hum_anna_annotated.bed")
	parse=parser.add_option("-S","--system_path",dest="system_path",default="/media/Backup/indexes/hg19/genome.fa")
	(opts,args)=parser.parse_args()

	if args:
		ff=open(args[0],"r")
	else:
		ff=sys.stdin
	c_lines=ff.readlines()
	ff3=open(opts.gene_models)
	g_lines=ff3.readlines()
	ff4=open(opts.annotation)
	a_lines=ff4.readlines()
	for c in cds(a_lines,g_lines,c_lines,opts.system_path,opts.name_format):
		third_codon(*c)

