#!/usr/bin/python

##########################################################
#Porject: Calling_Cards for HUMAN
#Data:2015/03/17
#Lastest update: 2016/01/16
#Author:JINCHUN ZHANG
#
#Input: txt data with format'chromosome	position count'
#Output: 24 .txt format data for 24 chromosomes, sorted by position, 3rd column is their number of independent insertion with 5^N bonus for N>1
#	Header is:
#		Choromsome
#		Position 
#       Independent_Insertions_withBonus
#		Independent_Insertions
#Cutoff is used for filter out the independent insertions based on their own count which is the 3rd column of the input dataset
#Usage is python Precut.py --e <.././experiment.gnashy> --w <../../wild.gnashy> --o <../../OUTFOLDER/> --C <integer>
#########################################################

import argparse
import numpy as np
import os
import subprocess

parser = argparse.ArgumentParser(prog='Precut.py', description='Count the Independent Insertions for each position and cut the file into 24 files.')
parser.add_argument("--e",dest="exp",type=str,nargs=1,help="Path and name for the experimental file")
parser.add_argument("--c",dest="cutoff", type=int, nargs=1, default=0,help="Set the cutoff for the read.counts, only keep Insertion when it has count bigger than the cutoff")
parser.add_argument("--w",dest="wt", type=str, nargs=1,help="Path and name for wild type dataset")
parser.add_argument("--o",dest="out", type=str, nargs=1,default="Output",help="Path and name for output folder, <../../female_Output>")

args=parser.parse_args()
cutoff=args.cutoff[0]
exp=args.exp[0]
wt=args.wt[0]
out=args.out[0]

print "Precut.py start!"
print "Cutoff used in insertion filtering is: ", cutoff

for x in ("EXP","WT") : 
	if not os.path.exists(out+'/PrecutOut/'+x):
        	os.makedirs(out+'/PrecutOut/'+x)
	if x=="EXP":
		infile=exp
	elif x=="WT":
		infile=wt
	sdin=open(infile,'r')
	inf=sdin.readlines()
	sdin.close()
	chr={
            1:[],2:[],3:[],4:[],5:[],6:[],7:[],8:[],9:[],10:[],11:[],12:[],
            13:[],14:[],15:[],16:[],17:[],18:[],19:[],20:[],21:[],22:[],23:[],24:[]
  	  }
	#filter out based on cutoff
	for i in inf:
        	a=i.split('\t',2)
        	a[0]=int(a[0])
        	a[1]=int(a[1])
        	a[2]=int(a[2].replace('\n',''))
        	if a[2]>cutoff:
            		chr[a[0]].append([a[1],a[2]])
	#count insertion events 
   	for i in chr:
		if len(chr[i])>0:
	        	chr[i]=np.array(chr[i])
        		uniq=np.unique(chr[i][:,0])
        		new=np.zeros((len(uniq),4),dtype='int32')
	        	new[:,1]=uniq
			new[:,0]=i
        		for n in range(len(uniq)):
            			freq=chr[i][chr[i][:,0]==uniq[n]]
	            		freq=freq[:,0]
				new[n,3]=len(freq)
				if len(freq)==1:
                    			new[n,2]=1
	                	else:
        	            		new[n,2]=len(freq)+5**len(freq)		
        		outf=out+'/PrecutOut/'+x+'/chr'+str(i)+'.txt'
			with open(outf,'w') as of:
				of.write('Chromosome\tPosition\tIndependent_Insertions_withBonus\tIndependent_Insertions\n')
        			np.savetxt(of, new,fmt="%d")			

subprocess.call("tail -n +2 -q "+out+"/PrecutOut/WT/chr*.txt> "+ out+"/PrecutOut/WT/combine.txt", shell=True)
print "Precutchr.py completed."
