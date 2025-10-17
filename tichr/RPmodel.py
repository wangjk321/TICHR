import numpy as np
import math
import pandas as pd
import os
import subprocess
import time
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

from .preprocess_hic import *


codepath = os.path.dirname(os.path.realpath(__file__))


def wrapped_calculateRP_wang(args):
    index, calculateRP_wang_func, pbar = args
    result = calculateRP_wang_func(index)
    pbar.update(1)
    return result


def extract_values_from_bedgraph_pandas(df, chromosome, start, end):
    # df 是给定bedgraph文件，
    start=start-1
    #end=end+1

    # 过滤特定染色体的数据
    #！！！！！！！！！！！！！！这一步特别特别慢
    df['chrom'] = pd.Categorical(df['chrom'])
    df_chrom = df[df['chrom'] == chromosome]
    
    # 过滤特定区间的数据
    df_filtered = df_chrom[(df_chrom['end'] >= start) & (df_chrom['start'] <= end)]
    
    # 提取每个bp的值

    result_length = end - start
    values = np.zeros(result_length)

    for _, row in df_filtered.iterrows():
        overlap_start = max(row['start'], start)
        overlap_end = min(row['end'], end)

        # Calculate the relative positions in the result array
        rel_start = overlap_start - start
        rel_end = overlap_end - start

        values[rel_start:rel_end] = row['value']
    
    return values

def calcu_weight(z,decay_distance,rpweight): #z是正数
    if rpweight == 'marge':
        lamda = -math.log(1/3)*1e5/decay_distance
        weight = 2*math.exp(-lamda*z/1e5)/(1.0+math.exp(-lamda*z/1e5)) 
    elif rpweight == 'classic':
        weight = 2 ** -(z/decay_distance)
    elif rpweight == 'powerlaw':
        lamda = math.log(0.5)/math.log(decay_distance)
        weight = (z+0.1)**lamda
    return(weight)

def obtain_RP_bdgdf(bamfilelist,codepath,gtfile,species):
    bdgdf_values = []
    for bamfile_i in bamfilelist:
        subprocess.run(["bash", codepath + "/bashcode/preprocessRP.sh", bamfile_i, gtfile,species], 
                stdout=open("info.txt", "w"), stderr=subprocess.STDOUT)
        bdgdf = pd.read_csv("midfiletmp/inputbam_treat_pileup.bdg", sep='\t', 
                                        header=None, names=['chrom', 'start', 'end', 'value'])
        bdgdf_values.append(list(bdgdf['value']))
    
    bdgdf['value'] = np.array(bdgdf_values).mean(axis=0)

    return(bdgdf)

def obtain_RPEP_candidatebeddf(bamfilelist,codepath,gtfile,species,candidatebed,generef,
                               quantileref,signaltype,separatepromoter,quantile_method):
    if not generef: generef=''
    peakdf_values = []
    for bamfile_i in bamfilelist:
        print("computing "+str(bamfile_i))
        subprocess.run(["bash", codepath + "/bashcode/preprocessRP.sh", bamfile_i, gtfile, species,candidatebed,generef], 
                       stdout=open("info.txt", "w"), stderr=subprocess.STDOUT)
        candidatebeddf = pd.read_csv("midfiletmp/candidatepeak.read.final",sep='\t',header=None)
        peakdf_values.append(list(candidatebeddf[3]))
    candidatebeddf[3] = np.array(peakdf_values).mean(axis=0)

    if not quantile_method:
        candidatebeddf = candidatebeddf
    elif quantileref and signaltype:
        print("Using reference file to make quantiles of ChIP-seq signals")
        refdf = pd.read_csv(quantileref, sep='\t')
        if separatepromoter:
            print("separate promoter")
            candidatebeddf = quantile_normalize_separate(candidatebeddf,refdf,signaltype,method=quantile_method)
        else:
            print("not separate promoter")
            candidatebeddf = quantile_normalize_any(candidatebeddf,refdf,signaltype,method=quantile_method)
    else:
        print("Please specify singaltype as DNase or H3K27ac")
        exit(1)

    return(candidatebeddf)


class myRP:
    def __init__(self,genebedfile,bamfilelist,gtfile,padding = int(1e5),decay_distance=10000,rpweight="classic",species='hs',
                 predictEP=False,candidatebed=None,rpresolution=1,tssrange=500,bamfilelist2=None,
                 quantileref=None,signaltype=None,separatepromoter=False,quantile_method=None,signaltype2=None,generef=None):
        self.geneDF= pd.read_csv(genebedfile,sep='\t',header=None) #chr,start,end,symbol,geneid,strands 第五列实在不行换其他的也行

        print("Extracting read from raw read files...")
        #os.system("bash "+codepath+"/bashcode/preprocessRP.sh "+bamfile+" "+gtfile+" "+species+" &> info.txt") 
        if not candidatebed:
            bdgdf = obtain_RP_bdgdf(bamfilelist,codepath,gtfile,species)
            if bamfilelist2:
                bdgdf2 = obtain_RP_bdgdf(bamfilelist2,codepath,gtfile,species)
                bdgdf['value'] = (bdgdf['value'] * bdgdf2['value']) ** 0.5
            self.bdgdf = bdgdf
        else:
            print("Extracting one read file over enhancer sites...")
            candidatebeddf = obtain_RPEP_candidatebeddf(bamfilelist,codepath,gtfile,species,candidatebed,generef,
                                                        quantileref,signaltype,separatepromoter,quantile_method)
            if bamfilelist2:
                candidatebeddf2 = obtain_RPEP_candidatebeddf(bamfilelist2,codepath,gtfile,species,candidatebed,generef,
                                                             quantileref,signaltype2,separatepromoter,quantile_method)
                candidatebeddf[3] = (candidatebeddf[3] * candidatebeddf2[3]) ** 0.5
            candidatebeddf[0] = pd.Categorical(candidatebeddf[0])
            self.candidatebeddf = candidatebeddf

        print("Extracting read finished")

        #os.system("rm -rf midfiletmp")
        self.gtfile = gtfile
        self.gtDF = pd.read_csv(gtfile,sep='\t',header=None)
        self.padding = padding
        self.rpresolution = rpresolution
        self.tssrange = tssrange
        self.rpweight = rpweight
        self.decay_distance = decay_distance

        self.weight = [calcu_weight(abs(z),decay_distance,rpweight) for z in range(-padding,padding+1,rpresolution)]
    
    def calcuRPEP(self,i,logrp=False,setpromoter1=True):
        if i % 500 == 0: print(str(i)+" genes processed.")

        genei = self.geneDF.iloc[i]
        chri = genei[0]
        strand = genei[5]
        tssi = genei[1] if strand == "+" else genei[2]

        geneipeaklist = self.candidatebeddf[(self.candidatebeddf[0]==chri) & (self.candidatebeddf[2]> tssi-self.padding) & (self.candidatebeddf[1]< tssi+self.padding)]
        geneipeakcenterlist = (geneipeaklist[1]+geneipeaklist[2])/2
        ifpromoter = (geneipeaklist[2]> tssi-self.tssrange) & (geneipeaklist[1]< tssi+self.tssrange)

        geneiweightlist = [calcu_weight(abs(i-tssi),self.decay_distance,self.rpweight) for i in geneipeakcenterlist]
        geneirp = np.array(geneiweightlist)*np.array(geneipeaklist[3])
        if logrp:
            geneirp_iflog=np.log1p(geneirp)
        else:
            geneirp_iflog = geneirp
        
        if geneirp.sum()>0: #判断是否全为0
            geneirp_percent = geneirp_iflog / (geneirp_iflog.sum())
        else:
            geneirp_percent = geneirp_iflog * 0

        if setpromoter1:
            geneirp_percent[ifpromoter] = 1 #启动子设置为1

        geneirp_df=geneipeaklist.iloc[:,0:4].copy()
        geneirp_df.columns = ["peakchr",'peakstart','peakend','activity']
        geneirp_df["genesymbol"] = genei[3]
        geneirp_df["genechr"] = genei[0]
        geneirp_df["genestart"] = genei[1]
        geneirp_df["geneend"] = genei[2]
        geneirp_df["genestrand"] = genei[5]
        geneirp_df["geneid"] = genei[4]
        geneirp_df["rpweight"] = list(geneiweightlist)
        geneirp_df["rp_rawvalue"] = list(geneirp)
        geneirp_df["rp_percent"] = list(geneirp_percent)

        return(geneirp_df)
    
    def calcuRPEP_allgene(self,logrp=False,setpromoter1=True):
        rpepdf = pd.DataFrame()
        for genei in range(self.geneDF.shape[0]):
            rpepdf_genei = self.calcuRPEP(genei,logrp=logrp,setpromoter1=setpromoter1)
            rpepdf = pd.concat([rpepdf,rpepdf_genei], ignore_index=True)
        return(rpepdf)

    def calculateRP_wang(self,i): #计算第i个基因的RP
        if i % 500 == 0:
            print(str(i)+" genes processed.")

        genei = self.geneDF.iloc[i]
        chri = genei[0]
        strand = genei[5]
        tssi = genei[1] if strand == "+" else genei[2]

        around_starti = max(tssi - self.padding,0)
        chrisize = int(self.gtDF[self.gtDF[0] == chri][1])
        around_endi = min(tssi + self.padding,chrisize)
        fix_bins = 2 * self.padding +1
        
        content = extract_values_from_bedgraph_pandas(self.bdgdf, chri, around_starti, around_endi)[::self.rpresolution]

        if len(content) > 0:
            signallist = np.array(content,dtype=np.float64)
            singallist_final =np.hstack((
                    np.zeros(self.padding-(tssi-around_starti))[::self.rpresolution],
                    signallist,
                    np.zeros(self.padding-(around_endi-tssi))[::self.rpresolution]
                    ))
        else:
            singallist_final = np.zeros(fix_bins)

        if np.sum(singallist_final) == 0:
            return(0)
            #return(np.sum(singallist_final))
        else:
            return(np.dot(singallist_final,self.weight))
            #return(np.sum(singallist_final))

    
    def calculateRP_wang_allgene(self,threads=8):
        numgene = len(self.geneDF)
        if threads > 1:
            pool = Pool(threads)
            rplist=pool.map(self.calculateRP_wang, range(numgene))
            # 关闭进程池
            pool.close()
            pool.join()
        elif threads==1:
            print("using 1 threads")
            rplist=[]
            for x in range(numgene):
                rplist.append(self.calculateRP_wang(x))
        else:
            print("please use positive number of threads")
            exit(1)


        outdf = self.geneDF.copy()
        outdf['rplist']=rplist
        
        return(outdf)

    def calculateRP(self,i): #计算第i个基因的RP
        if i % 500 == 0:
            print(str(i)+" genes processed.")

        genei = self.geneDF.iloc[i]
        chri = genei[0]
        strand = genei[5]
        if strand == "+":
            tssi = genei[1]
        else:
            tssi = genei[2]

        around_starti = tssi - self.padding
        around_endi = tssi + self.padding
        chrisize = int(self.gtDF[self.gtDF[0] == chri][1])
        if around_starti<0:
            around_starti=0
        if around_endi>chrisize:
            around_endi = chrisize
        
        fix_bins = self.padding + self.padding +1
        around_bins = around_endi - around_starti + 1
        
        cmdcode = codepath+"/bashcode/bigWigSummary"+" "+self.bwfile+" "+chri+" "+str(around_starti)+ \
                    " "+str(around_endi)+" "+str(around_bins) +" 2>> info.txt"    
        
        try:
            content = os.popen(cmdcode).readlines()
            if content:
                signallist = content[0].strip().replace('n/a','0').split('\t')
                signallist = np.array(signallist,dtype=np.float64)
                singallist_final =np.hstack((np.zeros(self.padding-(tssi-around_starti)),
                        signallist,
                        np.zeros(self.padding-(around_endi-tssi))))
            else:
                singallist_final = np.zeros(fix_bins)
                
            return(np.dot(singallist_final,self.weight))
        except:
            singallist_final = np.zeros(fix_bins)
            return(np.dot(singallist_final,self.weight))
    
    def calculateRP_allgene(self,threads=8):
        numgene = len(self.geneDF)

        if threads > 1:
            pool = Pool(threads)
            # 使用 map 函数并行处理任务，并获取结果
            rplist=pool.map(self.calculateRP, range(numgene))
            # 关闭进程池
            pool.close()
            pool.join()
        elif threads==1:
            print("using 1 threads")
            rplist=[]
            for x in range(numgene):
                rplist.append(self.calculateRP(x))
        else:
            print("please use positive number of threads")
            exit(1)

        outdf = self.geneDF.copy()
        outdf['rplist']=rplist

        os.system("rm -rf midfiletmp")
        return(outdf)