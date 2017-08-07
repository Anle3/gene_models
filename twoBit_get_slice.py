import twobitreader 
import string

def fetch_subseq(path,chrom,start,end,strand):
	genome = twobitreader.TwoBitFile(path)
	subseq= genome[chrom].get_slice(start,end)

	if strand=="-":
			subseq=reverse(subseq)
	 
	return subseq

#Reverse sequence if it comes from minus strand
def reverse(seq):
	 
	seq=seq[::-1].translate(string.maketrans('AUTGCautgc','TAACGTAACG'))
	return seq

 
			
if __name__=="__main__":
	from optparse import OptionParser
	parser=OptionParser()
	parser.add_option("-p","--system_path",dest="system_path")
	parser.add_option("-c","--chrom",dest="chrom")
	parser.add_option("-s","--start",type=int,dest="start")
	parser.add_option("-e","--end",type=int,dest="end")
	parser.add_option("-r","--strand",dest="strand",default="+")
		
	(opts,args)=parser.parse_args()
	  
	print  fetch_subseq(opts.system_path,opts.chrom,opts.start,opts.end,opts.strand)
						

