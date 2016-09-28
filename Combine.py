#!/bin/python

##########################################################
#Porject: Calling_Cards for HUMAN
#Date:2015/03/17
#Author:JINCHUN ZHANG
##########################################################
'''
Used to combine the Q-T clustered experimental data and the corresponding WILD data. 
Input: the directory of QTcluster.py output
Output: One csv file,each row is one cluster, with format:
	Chromosome
	Cluster
	Start
	Stop
	BkgHop
	BkgHop_withBonus
	ExpHop
	ExpHop_withBonus
'''
 
import pandas as pd
import numpy as np
import os
import sys
import argparse

parser = argparse.ArgumentParser(add_help=True)
parser = argparse.ArgumentParser(prog='Combine.py', description="Reformat the QToutput data to Chr, Start, End")
parser.add_argument("--o",dest="out", type=str, nargs=1, help="Type the directory of Output, <../../Output/>")

args=parser.parse_args()
output=args.out[0]

if not os.path.exists(output):
	print "Wrong Directory, please put the one ended with Output/"
	sys.exit("Ending Script")

#Outer join two dataframe
def _COMBINE_ (exp,bkg):
        exppos=exp['Position'].tolist()
        window=exp['Cluster'].tolist()
        window=pd.Series(window,index=exppos)
        expcount=exp['Independent_Insertion_withBonus'].tolist()
	expcount2=exp['Independent_Insertion'].tolist()
        expcount=pd.Series(expcount,index=exppos)
	expcount2=pd.Series(expcount2,index=exppos)
        bkgcount=pd.Series(list(bkg[:,2]),index=bkg[:,1])
	bkgcount2=pd.Series(list(bkg[:,3]),index=bkg[:,1])
        pos=np.append(np.array(exppos),bkg[:,1])
        pos=np.unique(pos)
        pos=pd.Series(list(pos),index=list(pos))
        new=pd.DataFrame({'Position':pos,'Cluster':window,'ExpHop_withBonus':expcount,'BkgHop_withBonus':bkgcount,'ExpHop':expcount2,'BkgHop':bkgcount2})
        new.index=range(len(new))
        return new

#Fill out the cluster number for each row
def _FILL_ (comb,clust):
        for p in range(1,clust+1):
		ind=comb[comb['Cluster']==p].index
                if len(ind)>2 and (ind[-1]-ind[0]-len(ind)) != -1:
			for a in range(ind[0],ind[-1]+1):
				comb.loc[a,'Cluster']=p
                comb=comb.fillna(0)
	return 	comb

#Transfer the dataset into the output format
def _TRANS_ (comb,clust):
        data=comb[comb['Cluster']==99999]
        data.insert(0,'Chromosome',c)
	s=data['Cluster']
	del data['Cluster']
	data.insert(1,'Cluster',s)
        data.insert(2,'Start',data['Position'])
        data.insert(3,'Stop',data['Position'])
        del data['Position'] 
        idx=data.index[-1]
	data=data.astype('int64')
        for grp in range(1,clust+1):
                block=comb[comb['Cluster']==grp]
                N=np.array(block['Position'])
                data.loc[idx+grp,'Chromosome']=c
                data.loc[idx+grp,'Start']=np.amin(N)
                data.loc[idx+grp,'Stop']=np.amax(N)
                data.loc[idx+grp,'BkgHop']=sum(block['BkgHop'])
                data.loc[idx+grp,'ExpHop']=sum(block['ExpHop'])
                data.loc[idx+grp,'BkgHop_withBonus']=sum(block['BkgHop_withBonus'])
                data.loc[idx+grp,'ExpHop_withBonus']=sum(block['ExpHop_withBonus'])
		data.loc[idx+grp,'Cluster']=grp
	data.index=range(len(data))
	return data

#Translate format when there is only one cluster exist. 
def _TRANS1_(comb):
        comb=comb.fillna(0)
        data=comb[comb['Cluster']==99999]
        data.insert(0,'Chromosome',c)
        s=data['Cluster']
        del data['Cluster']
        data.insert(1,'Cluster',s)
        data.insert(2,'Start',data['Position'])
        data.insert(3,'Stop',data['Position'])
        del data['Position']
        idx=data.index[-1]
        data=data.astype('int64')
        return data

#Process for each chromosome
print 'Combine.py Started!'
outf=output+'QTout/CCTest_Combine.csv'
counter=0

for c in range(1,25):
	expf=str(output)+'QTout/chr'+str(c)+'clustered.csv'
	bkgf=str(output)+'PrecutOut/WT/chr'+str(c)+'.txt'
	if os.path.isfile(expf):
		counter=counter+1
		exp=pd.read_csv(expf)
        	clust=np.array(exp['Cluster']).tolist()
	        clust=np.unique(clust)[-2]
        	bkg=np.loadtxt(bkgf,skiprows=1,dtype=int)
                comb=_COMBINE_(exp,bkg)
                if len(np.unique(clust))==1:
                        comb=_TRANS1_(comb)
                if len(np.unique(clust))>1:
                        clust=np.unique(clust)[-2]
                        comb=_FILL_(comb,clust)
                        comb=_TRANS_(comb,clust)
		comb=comb.astype('int64')
		if counter==1:
        	        comb.to_csv(outf,header=True,index=False)
	        else:
        		comb.to_csv(outf,header=False,index=False,mode='a')
		counter=counter+1
print 'Combine.py Completed!'
