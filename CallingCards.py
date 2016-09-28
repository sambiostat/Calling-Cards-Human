#usr/bin/python

##########################################################
#Porject: Calling_Cards
#Description: This is the pipeline used to call peaks in human calling cards experiment 
#Data:2015/03/17
#Last Update: 2016/01/05
#Author:JINCHUN ZHANG
##########################################################
import subprocess
import argparse

parser = argparse.ArgumentParser(prog='CallingCards.py', description="Calling Cards Wrap Script. Usage is python ../CallingCards.py --p <../../CODE_PATH/> --e <../../exp.gnashy> --w <../../wild.gnashy> --c <integer> --d <integer> --f <float> --a <alpha_value> --r <ref_option> --o <../../OUTFOLDER/>. ")
parser.add_argument("--p",dest="path", type=str, default=["./"],nargs=1, help="Type the Calling_Cards Code directory,<../../CODE/> ")
parser.add_argument("--e",dest="exp", type=str, nargs=1, help="Type the path and name experiment data, <../../expriment.gnashy")
parser.add_argument("--w",dest="wt", type=str, nargs=1,help="Type the path and name of the wild type data, <../../wild.gnashy")
parser.add_argument("--c",dest="cutoff", type=str,default=["0"], nargs=1, help="Set the cutoff(Default value is 0) for the read.counts, only keep Insertion when it's count>=cutoff")
parser.add_argument("--d",dest="dist", type=str, default=["2500"], nargs=1, help="Set the maximum within cluster distance(Default value is 2500bp)")
parser.add_argument("--f",dest="pseudo", type=str, default=["0.05"], nargs=1, help="Set the pseudocount_value(Default value is 0.05) for the statistic testing")
parser.add_argument("--a",dest="alpha", type=str, default=['0.05'], nargs=1, help="Set the alpha value (Default value is 0.05, range from 0 to 1) for statistical testing, when P_value smaller than this alpha value we say its significant.") 
parser.add_argument("--r",dest="ref", type=str, nargs=1, help="Choose the reference genome used to annotate output clusters, it has to be one of the following <mm9; mm10; hg18; hg19>")
parser.add_argument("--o",dest="out", type=str, nargs=1, default=["Output"],help="Set the name of output folder, such as <../../Female_Output/>")

args=parser.parse_args()
code_pwd=args.path[0]
cutoff=args.cutoff[0]
d=args.dist[0]
f=args.pseudo[0]
a=args.alpha[0]
exp=args.exp[0]
wt=args.wt[0]
output=args.out[0]
ref=args.ref[0]

ref_pwd="/home/comp/rmlab/ssankararaman/Samantha/Calling_Cards/Calling_Cards_Ref/"
ref_genome=ref_pwd+"refGene."+ref+".Sorted.bed"

print 'Calling_Cards Pipeline for Human Started!'

subprocess.call("python "+code_pwd+"Precutchr.py  --c " + cutoff+" --e "+exp+" --w "+wt+" --o "+output, shell=True)
subprocess.call("python "+code_pwd+"QTcluster.py --o "+output+" --d "+ d, shell=True )
subprocess.call("python "+code_pwd+"Combine.py --o "+output,shell=True)
subprocess.call("python "+code_pwd+"Simpletest.py --f " +f+" --a "+a+" --o "+output, shell=True)
subprocess.call("python "+code_pwd+"Annotated.py --r "+ref_genome+" --o "+output,shell=True)
#subprocess.call("mv "+output+"/Annotation/Total_Insertions_P_value_Full_Annotated.txt "+output+"/Annotation/"+exp+"_Insertions_Annotated.txt",shell=True)
print 'Calling_Cards Pipeline Completed!'
print 'Please check the User-specified Output Folder, the final annotated data with P_value is saved at the Annotation folder:)'
