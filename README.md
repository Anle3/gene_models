 
## Gene Models
 
### Description  
Gene models provide a toolbox for mapping and extracting transript level information to and from genome.

### Dependencies  
twobitreader

### Input files  
A path to a 2bit formated genome and a tab delimited file with gene_models downloaded from UCSCS Table Browser (RefSEq, GENCODE)

### Download input files
####TwoBit Genome
Form terminal install 
```
pip install ucscgenome
#for local install on servers use --user if you don't have rights to install globaly
pip install --user ucscgenome
```
```
#!/usr/bin/env python
#from ucscgenome import Genome 
hg19 = Genome('hg19')
```
This will install 2bit in ~/.ucscgenome
```
ls ~/.ucscgenome

hg19.2bit
```
 
 ####Gene models
 * From UCSC table browser select assembly (e.g hg19)
 * Select gene model track (e.g. GENCODE Genes V19)
 * Select output format: all fields from selected table
 * Savbe in output file
 
### Install gene_models
```pip install git+https://github.com/Anle3/gene_models.git```  

### Examples
```
#Import gene models module
from gene_models import *

#define paths
genecode_models= "path_to_gene_models/hg19_gene_models.txt"
genome="path_to_2bit_genome/hg19.2bit"

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
