import os
from multiprocessing import Pool
import pandas as pd
import numpy as np
import os,random, string,sys

class paralfunJuicer(object):
    def __init__(self,myfun,num_processer):
        self.myfun = myfun
        self.num_processer = num_processer

    def run(self,data,normalize,resolution,gt,outname,juicer,chrolist):
        #chrlist= list(range(1,maxchr+1)); chrlist.append("X")
        #chrolist = ["chr"+str(a) for a in chrlist]

        resultlist=[]
        p = Pool(self.num_processer)
        for r in chrolist:
            p.apply_async(self.myfun, args=(r,data,normalize,resolution,gt,outname,juicer))
        p.close()
        p.join()

def oneJuicer(chrom,data,normalize,resolution,gt,outname,juicer):
    #print("Input data: ", data)
    print("Dumping contact " +chrom+ " matrix from .hic file ......")
    randomindex=''.join(random.sample(string.ascii_letters + string.digits, 8))

    codepath = os.path.dirname(os.path.realpath(__file__))
    makeIntra = codepath+"/bashcode/makeMatrixIntra.sh"
    if not juicer:
        juicer = codepath+"/jc/jctool_1.11.04.jar"
    foldername = outname
    os.system("bash "+makeIntra+" "+normalize+" "+"."+" "+data+" "+
            str(resolution)+" "+gt+" "+juicer+" "+chrom+" "+foldername + "> info.txt"+randomindex)
    try: os.system("rm info.txt"+randomindex)
    except: pass
    print("Dump finished, output is in ./"+outname)

def manyJuicer(data,normalize,resolution,gt,outname,juicer,chrolist,threads=24): #chrolist 为类似["chr21","chr22"]
    paralfunJuicer(oneJuicer,threads).run(data,normalize,resolution,gt,outname,juicer,chrolist)
