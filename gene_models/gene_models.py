#! /usr/bin/env python
from table_loader import load_table
from twoBit_get_slice import fetch_subseq
from twobitreader import *
from bisect import *
from collections import namedtuple
class ExonChains(object):
	
	#initialize Transcript class with all fields from gene_models table followed by the path where the files for the given model system are stored
	def __init__(self,bin_,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonS,exonE,score,name2,cdsStartStat,cdsEndStat,exonFrames,system_path):
		import os
		self.bin_=bin_
		self.genome=system_path
		self.name=name
		self.cdsStartStat=cdsStartStat
		self.cdsEndStat=cdsEndStat
		self.gene_id=name2
		self.score=score
		self.chrom=chrom
		self.strand=strand
		self.txStart=int(txStart)
		self.txEnd=int(txEnd)
		self.cdsStart=int(cdsStart)
		self.cdsEnd=int(cdsEnd)
		self.exonCount=int(exonCount)
		self.exonFrames_=exonFrames
		self.dir_=int("%s1"%strand)
 
		self.exonStarts=[int(n) for n in exonS.split(",") if n!=""]
		self.exonEnds=[int(n) for n in exonE.split(",") if n!=""]

	@property
	def exonFrames(self):
		exonFrames=[int(n) for n in self.exonFrames_.split(",") if n!=""]
		return exonFrames


	@property
	def exons(self):

		exons=list()

		for es,ee in zip(self.exonStarts,self.exonEnds):
			exons.append(fetch_subseq(self.genome,self.chrom,es,ee,self.strand)) #fetch subseq from genome
		return exons

	@property 
	def spliced_seq(self):
	
		sq=""
		exons_d=self.exons[::self.dir_]
		for e in exons_d:
			sq+=e
		return sq

	@property
	def spliced_len(self):
		return len(self.spliced_seq)


	@property
	def exon_lengths(self):
		lengths=[len(e) for e in self.exons]
		return lengths
	
	@property
	def spliced_exon_starts(self):
		k=[0]
		for l in self.exon_lengths[:-1]:
			k.append(k[-1]+l)
		return k
	
	@property
	def spliced_cds_coors(self):
		CdsC=namedtuple("CdsC",["start","end"])	
		cdsS=self.cdsStart-self.txStart
		cdsE=self.cdsEnd-self.txStart
		length=self.txEnd-self.txStart
		for i in range(1,len(self.exonStarts)):
			if self.cdsStart >= self.exonStarts[i]:
				cdsS=cdsS-(self.exonStarts[i]-self.exonEnds[i-1])
			if self.cdsEnd>=self.exonStarts[i]:
				cdsE=cdsE-(self.exonStarts[i]-self.exonEnds[i-1])
			length=length-(self.exonStarts[i]-self.exonEnds[i-1])
#format start,e according to +/- strand
		s,e=[cdsS,cdsE][::self.dir_]
		le=[0,length][::self.dir_]
		s=(s-le[0])*self.dir_
		e=(e-le[0])*self.dir_
		return CdsC(s,e)




	def slice_seqs(self,start,end):
		return  self.spliced_seq[start:end]


	def sub_exon_seq(self,**args):
		seq_l=list()
		try:
			for a in args:
				seq_l.append(self.exons[a])
				
				
		except IndexError:
			print "exon number out of bounds"
		else:	
			seq="".join(seq_l[::self.dir_])
			return seq

	def intersect(self,start,end):
		fl=bisect_left(self.exonStarts,start)
		fr=bisect_right(self.exonEnds,start)
		#el=bisect_left(self.exonStarts,end)
		#er=bisect_right(self.exonStarts,end)
		el=bisect_left(self.exonStarts,end)
		er=bisect_right(self.exonEnds,end)
	#transcript before coors
		try:
			b_exon_starts=self.exonStarts[:fl]
			b_exon_ends=self.exonEnds[:fl]
			b_exon_count=len(b_exon_starts)
			b_exon_ends[-1]=min(start,b_exon_ends[-1])
			b_tx_s=b_exon_starts[0]
			b_tx_e=b_exon_ends[-1]
			if b_tx_s==b_tx_e:
				raise IndexError
			b_cds_s=max(self.cdsStart,b_exon_starts[0])
			b_cds_e=min(self.cdsEnd,b_exon_ends[-1])
			before= ExonChains(self.bin_,self.name,self.chrom,self.strand,b_tx_s,b_tx_e,b_cds_s,b_cds_e,b_exon_count,",".join(map(str,b_exon_starts)),",".join(map(str,b_exon_ends)),self.score,self.gene_id,self.cdsStartStat,self.cdsEndStat,self.exonFrames_,self.genome)
		except IndexError:
			before= None

		#	print "before",self.name,start,self.exonStarts[0]
	#transcipt between		
		try:
			#print fl,el,fr,er,self.exonStarts,self.exonEnds
			c_exon_starts=self.exonStarts[max(fr,0):el]
			c_exon_ends=self.exonEnds[max(fr,0):el]
			c_exon_count=len(c_exon_starts)
			c_exon_starts[0]=max(start,c_exon_starts[0])
			c_exon_ends[-1]=min(end,c_exon_ends[-1])

			c_tx_s=c_exon_starts[0]
			c_tx_e=c_exon_ends[-1]
			c_cds_s=max(self.cdsStart,c_exon_starts[0])
			if c_tx_s==c_tx_e:
				raise IndexError
			c_cds_e=min(self.cdsEnd,c_exon_ends[-1])
			chain= ExonChains(self.bin_,self.name,self.chrom,self.strand,c_tx_s,c_tx_e,c_cds_s,c_cds_e,c_exon_count,",".join(map(str,c_exon_starts)),",".join(map(str,c_exon_ends)),self.score,self.gene_id,self.cdsStartStat,self.cdsEndStat,self.exonFrames_,self.genome)
		except IndexError as exc:
			chain= None

	#transcript after
		try:
			#a_exon_starts=self.exonStarts[er-1:]
			a_exon_starts=self.exonStarts[er:]
			a_exon_ends=self.exonEnds[er:]
			#a_exon_ends=self.exonEnds[er-1:]
			a_exon_count=len(a_exon_starts)
			a_exon_starts[0]=max(end,a_exon_starts[0])
			#print self.name,end,a_exon_starts[0]
			a_tx_s=a_exon_starts[0]
			a_tx_e=a_exon_ends[-1]
			a_cds_s=max(self.cdsStart,a_exon_starts[0])
			if a_tx_s==a_tx_e:
				raise IndexError
			a_cds_e=min(self.cdsEnd,a_exon_ends[-1])
			after= ExonChains(self.bin_,self.name,self.chrom,self.strand,a_tx_s,a_tx_e,a_cds_s,a_cds_e,a_exon_count,",".join(map(str,a_exon_starts)),",".join(map(str,a_exon_ends)),self.score,self.gene_id,self.cdsStartStat,self.cdsEndStat,self.exonFrames_,self.genome)
		except IndexError:
			after= None
		#	print "after",self.name,start,end,self.exonStarts,self.exonEnds,a_exon_ends,f,e
		return before,chain,after


	def map_from_spliced(self,pos):

		if self.strand=="+":
			p=pos

		elif self.strand=="-":
			p=self.spliced_len-1-pos
		k=0	
		lengths=self.exon_lengths

		for i in range(0,len(lengths)):
			if k<=p<k+lengths[i]:
				m=p-k+self.exonStarts[i]
			k=k+lengths[i]
		try:
			return m

		except UnboundLocalError:
				return "NA"

	def map_to_spliced(self,pos):

		p=pos
		k=0
		lengths=self.exon_lengths

		for i in range(0,len(lengths)):
			if self.exonStarts[i]<=p<self.exonEnds[i]:
				m=p+k-self.exonStarts[i]
			k=k+lengths[i]
		try:	
			if self.strand=="+":
				m=m
			elif self.strand=="-":
				m=self.spliced_len-1-m
		
			return m
		except UnboundLocalError:
				return "NA"
	@property
	def as_bed(self):
		rgb={"+":"255,0,0","-":"0,0,255",}	
		es=",".join([str(l-self.txStart) for l in self.exonStarts])
		return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s"%(self.chrom,self.txStart,self.txEnd,self.name,self.score,self.strand,self.txStart,self.txEnd,rgb[self.strand],self.exonCount,",".join(map(str,self.exon_lengths)),es)

class Transcript(ExonChains):

	def __init__(self,bin_,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonS,exonE,score="0",name2=None,cdsStartStat=None,cdsEndStat=None,exonFrames=None,chr_pos=None,system_path=None):
		super(Transcript,self).__init__(bin_,name,chrom,strand,txStart,txEnd,cdsStart,cdsEnd,exonCount,exonS,exonE,score,name2,cdsStartStat,cdsEndStat,exonFrames,system_path)
		if cdsStart==cdsEnd:
			
			self.UTR3=None
			self.CDS=None
			self.UTR5=None
			self.descr="ncRNA"
		else:
			self.descr="coding"
			try:
				if self.strand=="+":
					self.UTR5,self.CDS,self.UTR3=self.intersect(self.cdsStart,self.cdsEnd)
				else:
					self.UTR3,self.CDS,self.UTR5=self.intersect(self.cdsStart,self.cdsEnd)
			except IndexError:
				print self.name,self.cdsStart,self.exonStarts[:2]
#generator for transcripts
def get_transcripts(file_,twoBit_path,cast="n"):
#choose model
	system_path=twoBit_path
	 
	if file=="-":
		fa=sys.stdin
	else:				
		fa=open(file_,"r")
	lines=fa.readlines()

 
 
	for tbl in load_table(lines,cast):
		tbl=list(tbl)
		
		yield Transcript(*tbl,system_path=system_path)


if __name__=="__main__":
	from optparse import OptionParser
	parser=OptionParser()
	parser.add_option("-p","--system_path",dest="system_path")
	parser.add_option("-c","--cast",default="n",dest="cast")
	
	(opts,args)=parser.parse_args()
	print 'track name="name" description="description" colorByStrand="255,0,0 0,0,255"'
	for tx in get_transcripts(args[0],opts.system_path,opts.cast):
		print tx.as_bed
				
