# Racon
Consensus module for raw de novo DNA assembly of long uncorrected reads.  

## Quick start
Clone and make Racon:
```  
git clone https://github.com/isovic/racon.git  && cd racon && make modules && make tools && make -j  
```
Run an example script:  
```  
./example1-lambda.sh  
```  
Tip: Running Racon iterativelly will produce better consensus sequences. (But don't forget to align your reads to the consensus sequence from the previous iteration.)  

## Dependencies
1. gcc >= 4.8  

## Installation
```  
git clone https://github.com/isovic/racon.git  
cd racon  
make modules  
make tools  
make -j  
```  

## Usage
```  
bin/racon [options] <raw_contigs.fasta> <alignments.sam> <out_consensus.fasta>  
```  
For detailed info on various options, run ```bin/racon``` without arguments.  
