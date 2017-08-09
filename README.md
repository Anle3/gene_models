 
## Gene Models
 
### Description  
Gene models provide a toolbox for mapping and extracting transript level information to and from genome.

### Dependencies  
twobitreader

### Input files  
A path to a 2bit formated genome and a tab delimited file with gene_models downloaded from UCSCS Table Browser (RefSEq, GENCODE)

### Install  
```pip install git+https://github.com/Anle3/gene_models.git```  

### Examples
```
#Import gene models module
from gene_models import *
#define paths
genecode_models= "/gpfs01/home/gecoo/omics/hg19/hg19_gene_models.txt"
genome=" /gpfs01/home/gecoo/omics/hg19/hg19.2bit"

#add Trascript object in a dictionary where names are the keys
transcripts={}
tx=get_transcripts(genecode_models,genome)
 for t in tx:
  transcripts[t.name]=t

##Operate in one specific transcript
t1=t["ENST00000237247.6"]
## get gene ID
 t1.gene_id
Out: 'SGIP1'
##get spliced sequence
t1.spliced_seq
##map trascript position to genome
  t1.map_from_spliced(100)
Out: 67000003 ##100 position on the transcript corresponds to 67000003 
```
