#!usr/bin/python
##########################################################
#Porject: Calling_Cards
#Data:2015/03/17
#Author:JINCHUN ZHANG
#
#    QTcluster.py used to partition experiment data into clusters based on the position.
#    Input: .txt format data, i.e. the output file of Precutchr.py.
#    Output: 24 csv files for 24 chromosomes with following header: 
#	Position
#	Cluster
#	Independent_Insertions_withBonus
#	Independent_Insertions
#	Chromosome
#
#Usage: python QTcluster.py --o <../../OUTFOLDER/> --d <Distance_integer>
#########################################################

print " QTcluster.py start!"
   
import numpy as np
import pandas as pd
import os
import sys
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser = argparse.ArgumentParser(prog='QTcluster.py', description='Use Q_T clustering method to cluster the insertion events based on position')
parser.add_argument("--o",dest="out", type=str, nargs=1, help="Type the directory of Output, <../../Output/>")
parser.add_argument("--d",dest="dist", type=int, nargs=1, help="Set the maximum within-cluster distance, <Integer>")

args=parser.parse_args()
output=args.out[0]
d=args.dist[0]

print "The maximum within-cluster distance is:", d

if not os.path.exists(output+'PrecutOut/EXP'):
	print "Wrong Directory, please put the one contains output from PrecutOut"
	sys.exit("Ending Script")	    

if not os.path.exists(output+'QTout'):
    os.makedirs(output+'QTout')

#Function _SUBCLUST, find subclusters of input points and within distance less than d
def _SUBCLUST_ (index,pos,d):
	aa=0
	pos=sorted(pos)
	while aa==0:
		left_index=min(index)
		right_index=max(index)
		sub=[pos[x] for x in index] 
		if left_index==0 and right_index<len(pos)-1 and pos[right_index+1]-pos[left_index]<= d:
			sub.append(pos[right_index+1])	
			index.append(right_index+1)
		elif right_index==len(pos)-1 and left_index>0 and pos[right_index]-pos[left_index-1]<=d:
			index.append(left_index-1)
			sub.append(pos[left_index-1])
		elif left_index>0 and right_index<len(pos)-1:
			left=pos[left_index-1] 
			right=pos[right_index+1]
			center=(pos[right_index]+pos[left_index])/2
			if center-left>right-center and right-pos[left_index]<=d:
					sub.append(right)
					index.append(right_index+1)
			elif center-left<right-center and pos[right_index]-left<=d:
					sub.append(left)
                        	        index.append(left_index-1)
			elif center-left==right-center and right-left<=d:
                                	sub.append(right)
					sub.append(left)
					index.append(left_index-1)
					index.append(right_index+1)
			else: aa=1
		else: aa=2
	return sub

#Function _MAX_, used to find the subclusters with maximum within-cluster distance, i.e. the biggest subcluster
def _MAX_ (pos):
	MAX=0
	MAX_sub=[]
	for p in range(len(pos)):
		index=[p]
		subp=_SUBCLUST_(index,pos,d)
		if max(subp)-min(subp)>=MAX:
			MAX=max(subp)-min(subp)
                        MAX_sub=subp
	return MAX_sub

#Function _QTCLUST_	
#	INPUT: pos:position set
#	OUTPUT: subclustered data with format 'position, number_of_cluster. p.s. 99999 refers to singleton 
def _QTCLUST_(pos):
        subnum=0
        pp=0
        qtclustered=[]
        while pp==0:
           	pos=sorted(pos)
                length=len(pos)
		if length==0: 
			pp=1
                if length==1:
                        pp=2
                        qtclustered.append([pos[0],99999])
                elif length>1:
			sub=_MAX_(pos)
			if len(sub)==1:
				qtclustered.append([sub[0],99999])
			else:
				subnum=subnum+1 
				for b in sub:
					qtclustered.append([b,subnum])
			pos=list(set(pos)-set(sub))
			pos=sorted(pos)
	return qtclustered

'''Read in one Chromosome once and process as following,
	Input: Chr sorted by pos with format 'chr position count'
	Output: save QTclustered dataframe data for each chr with format
		Chr
		Position
		Window
		Count
'''
	
for i in range(1,25):
        inf=str(output)+'PrecutOut/EXP'+'/chr'+str(i)+'.txt'
	if os.path.isfile(inf):
	        chr=np.loadtxt(inf,skiprows=1,dtype=int)
        	positions=chr[:,1]
		print 'chr'+str(i)	
		qtc=np.array(_QTCLUST_(positions))
		df=pd.DataFrame(qtc,index=qtc[:,0],columns=['Position','Cluster'])
		s1=pd.Series(chr[:,2],index=chr[:,1])
		s2=pd.Series(chr[:,3],index=chr[:,1])
		df['Independent_Insertion_withBonus']=s1
		df['Independent_Insertion']=s2
		df['Chromosome']=i
		df=df.sort(columns='Position',axis=0)	
		outf=output+'QTout/chr'+str(i)+'clustered.csv'
		df.to_csv(outf, index=False,header="Position,Cluster,Independent_Insertions_withBonus,Independent_Insertions,Chromosome")
print "QTcluster completed!"
