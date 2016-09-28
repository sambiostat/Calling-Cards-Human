#!/usr/bin/python
#Used to clean and sort the refGene.bed file.Input is the bed file download from UCSC genome browser

import numpy as np
import pandas as pd

infile="/home/comp/rmlab/ssankararaman/Samantha/Calling_Cards/HUMAN/Code/human.ref.bed"
outfile="/home/comp/rmlab/ssankararaman/Samantha/Calling_Cards/HUMAN/Code/refGene.Sorted.bed"

old=pd.read_csv(infile,sep='\t')
out=pd.DataFrame(columns=["chr","start","end","name","name2","strand"])
out["chr"]=old["chrom"]
out["start"]=old["txStart"]
out["end"]=old["txEnd"]
out["name"]=old["name"]
out["name2"]=old["name2"]
out["strand"]=old["strand"]

block=out[out["chr"].isin (["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"])]

result=block.sort(["chr","start","end"])
result.to_csv(outfile,sep='\t',header=False, index=False)
