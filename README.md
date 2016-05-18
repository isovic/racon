# consise
Consensus module for raw de novo DNA assembly of long uncorrected reads.  

## Quick start
Clone and make Consise:
```  
git clone https://github.com/isovic/consise.git  && cd consise && make modules && make tools && make -j  
```
Run an example script:  
```  
./example1-lambda.sh  
```  

## Dependencies
1. gcc >= 4.8  

## Installation
```  
git clone https://github.com/isovic/consise.git  
cd consise  
make modules  
make tools  
make -j  
```  

## Usage
```  
bin/consise [options] <raw_contigs.fasta> <alignments.sam> <out_consensus.fasta>  
```  
For detailed info on various options, run ```bin/consise``` without arguments.  
