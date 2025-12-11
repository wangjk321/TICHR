import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, time
import subprocess
import random,string
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from multiprocessing import Pool
from joblib import Parallel, delayed
from tqdm import tqdm
import shutil

from .preprocess_hic import *
from .tichr_function import *
from .highOrderStructure import *


def is_file_path_list(readFileList):
    if not isinstance(readFileList, list):
        return False
    
    for item in readFileList:
        if not isinstance(item, str):
            return False
        if not os.path.isfile(item):   # 必须是实际存在的文件
            return False
    
    return True


class Tichr:
    def __init__(self,candidatesite,readFileList,gtfile,candidateGeneFile,refGeneFile=None,
                 ifTSSrange=500,S2Gmax=100000,
                 hicfilepath=None,readFileList2=None):
        
        print("***Checking file availablity")
        if os.path.exists(candidatesite) or candidatesite in {"denovo_peak", "surronding_bin", "onlypromoter"}:
            self.candidatesite = candidatesite
        else:
            raise ValueError(
                f"Invalid candidatesite: {candidatesite}. "
                "Must be an existing file path or one of denovo_peak, surronding_bin, onlypromoter.")

        if not is_file_path_list(readFileList):
            raise ValueError("readFileList should be a list of file path")

        self.readFileList = readFileList
        self.readFileList2 = readFileList2

        self.gtfile = gtfile
        self.refgene_file = refGeneFile
        self.ifTSSrange = ifTSSrange
        self.candidateGeneDF = pd.read_csv(candidateGeneFile,sep='\t',header=None) #chr,start,end,symbol,geneid,strands 第五列实在不行换其他的也行
        self.candidateGeneFile = candidateGeneFile
        self.candidateGeneChrList = self.candidateGeneDF[0].unique()
        self.hicfilepath = hicfilepath
        self.peakToGeneMaxDistance = S2Gmax

        self.nomhicdf = None
        self.hicRes = None
        self.hicProcessedDataType = None
        self.ifUseHiCRef = None

        self.structureDF = None
        self.structureType = None
        self.structureWeight = None
    
    def makeSiteBed(self,macs2species='hs',binResolution=100,
                    blackregion=None,tmpdir=None,fixPeakWidth=True, only_promoter_area=100):
        self.macs2species = macs2species
        self.binResolution = binResolution
        if not tmpdir:
            self.tmpdir = "tichr_tmp_"+''.join(random.choices(string.ascii_letters, k=10))
        else:
            self.tmpdir = tmpdir
        self.candidatesite_file = makeSiteBedFunction(self.candidatesite,self.candidateGeneFile,self.readFileList,self.gtfile,
                                                      species=macs2species,binResolution=binResolution,
                                                      peakToGeneMaxDistance=self.peakToGeneMaxDistance,
                                                      blackregion=blackregion, refgene_file=self.refgene_file,tmpdir=self.tmpdir,
                                                      fixPeakWidth=fixPeakWidth, only_promoter_area=only_promoter_area)


    def makeSiteBdg(self, coverageMethod="coverageBed",spmr = False,multiBAMmerge='mean',file_type="bam",
                    quantileref=None,quantile_method=None,
                    separatepromoter=False,signaltype=None,signal2type=None,):
        #coverageMethod in ["macs2RP","coverageBed"]
        self.candidatesite_coverage = makeSiteBdgFunction(self.candidatesite_file,self.readFileList,self.gtfile,
                                                          coverageMethod,spmr = spmr,species=self.macs2species,
                                                          refgene_file=self.refgene_file,tmpdir=self.tmpdir,
                                                          ifTSSrange=self.ifTSSrange,
                                                          quantileref=quantileref,signaltype=signaltype,
                                                          separatepromoter=separatepromoter,quantile_method=quantile_method,
                                                          multiBAMmerge=multiBAMmerge,file_type=file_type)
        
        if self.readFileList2:
            candidatesite_coverage2 = makeSiteBdgFunction(self.candidatesite_file,self.readFileList2,self.gtfile,
                                                          coverageMethod,spmr = spmr,species=self.macs2species,
                                                          refgene_file=self.refgene_file,tmpdir=self.tmpdir,
                                                          ifTSSrange=self.ifTSSrange,
                                                          quantileref=quantileref,signaltype=signal2type,
                                                          separatepromoter=separatepromoter,quantile_method=quantile_method,
                                                          multiBAMmerge=multiBAMmerge,file_type=file_type)
            
            self.candidatesite_coverage[3] = (self.candidatesite_coverage[3] * candidatesite_coverage2[3]) ** 0.5
    
# hic 的标准化分为三个部分
# 1. hicNormType: Juicer本身的KR、VC_SQRT等标准化方法
# 2. contactNorm：标准化矩阵，如使用abc，还是oe方法，还是不进一步标准化
# 3. ifUseHiCRef：是否处以以reference hic

    def proceessHiC(self,hicRes,hicDataType,hicNormType,juicertool=None,
                    threads=8,contactNorm='default',
                    hicprefix="observed.VC_SQRT.",hicsuffix=".matrix.gz"):
        print("***Processing hic ...")
        self.hicRes=hicRes
        self.hicProcessedDataType = hicDataType
        self.contactNorm = contactNorm
        self.nomhicdf = gethicfile(self.hicfilepath,hicRes,hicDataType,self.candidateGeneChrList,
                              hicnorm=hicNormType,gt=self.gtfile,juicertool=juicertool,threads=threads,
                              contactNorm=contactNorm,hicprefix=hicprefix,hicsuffix=hicsuffix)        
    
    def weightStructure(self,structureType,structureFile,structureWeight):
        
        if all(isinstance(var, list) for var in [structureFile, structureType, structureWeight]):
            print("structureDF, structureType and structureWeight are all list")
            structureDF = []
            for file in structureFile:
                df = pd.read_csv(file, header=None, sep="\t")
                df[0] = pd.Categorical(df[0]) 
                structureDF.append(df)
            self.structureDF = structureDF
            self.structureType = structureType
            self.structureWeight = structureWeight
        elif all(isinstance(var, str) for var in [structureFile, structureType]):
            print("structureDF, structureType and structureWeight are all str")
            structureDF = pd.read_csv(structureFile,header=None,sep="\t")
            structureDF[0] = pd.Categorical(structureDF[0]) #make chromosome name to category
            self.structureDF = structureDF
            self.structureType = structureType
            self.structureWeight = structureWeight
        else:
            print("structureDF, structureType and structureWeight are not consistent style")
            exit(1)
        
    
    def computeGenei(self,i,weightType,rpDecayDistance=10000,fixedFunctionType='rp-classic',
                 given_gamma=1.024238616787792, given_scale = 5.9594510043736655,
                 ref_gamma = 0.87, ref_scale = -4.80 + 11.63 * 0.87, hicmindistance=5000,
                 logRgX=False,setpromoter1=False,ifUseHiCRef=True, goldWeightDf=None, noise_ratio=0,noise_quantile=0):
        
        
        t1 = time.time()
        # if i % 1000 == 0:
        #     print(str(i)+" genes processed.")
        
        genei = self.candidateGeneDF.iloc[i]
        geneichr = genei[0]
        geneistrand = genei[5]
        geneitss = genei[1] if geneistrand == "+" else genei[2]

        if self.candidatesite == 'surronding_bin':
            geneitss = (geneitss // self.binResolution) * self.binResolution

        geneipeaklist = self.candidatesite_coverage[(self.candidatesite_coverage[0]==geneichr) & 
                                                    (self.candidatesite_coverage[2]> geneitss-self.peakToGeneMaxDistance) & 
                                                    (self.candidatesite_coverage[1]< geneitss+self.peakToGeneMaxDistance)]
        geneiPeakxEpigenome = np.array(geneipeaklist[3])
        geneipeakcenterlist = (geneipeaklist[1]+geneipeaklist[2])/2
        ifpromoter = (geneipeaklist[2]> geneitss-self.ifTSSrange) & (geneipeaklist[1]< geneitss+self.ifTSSrange)

        # avoiding using hicdata wrongly
        if weightType == "hic": 
            hicProcessedData = self.nomhicdf
        else:
            hicProcessedData = None

        geneiPeakxWeightList = []
        for x in geneipeakcenterlist:
            geneiPeakxWeight = makeWeightFunction(weightType,x,geneitss,rpDecayDistance=rpDecayDistance,fixedFunctionType=fixedFunctionType,
                                    given_gamma=given_gamma, given_scale = given_scale,
                                    ref_gamma = ref_gamma, ref_scale = ref_scale, 
                                    hicmindistance=hicmindistance,hicProcessedData=hicProcessedData,hicRes=self.hicRes,
                                    geneChr=geneichr,hicProcessedDataType=self.hicProcessedDataType,ifUseHiCRef=ifUseHiCRef,
                                    peakToGeneMaxDistance=self.peakToGeneMaxDistance, goldWeightDf=goldWeightDf)
            
            sweight = weightStructureFunc(self.structureType,self.structureDF,self.structureWeight,x,geneitss,geneichr)
            geneiPeakxWeightList.append(geneiPeakxWeight * sweight)

        if fixedFunctionType == 'closest':
            geneiPeakxWeightList = (np.arange(len(geneiPeakxWeightList)) == np.argmin(geneiPeakxWeightList)).astype(int)

        #print(np.nansum(geneiPeakxWeightList))

        if weightType=="hic" and self.contactNorm=='abc':
            #qc gene
            badgene_threshold = 0.01 #排除没有链接的基因
            if np.nansum(geneiPeakxWeightList) < badgene_threshold:
                geneiPeakxWeightList = [getpowerlaw(abs(x-geneitss),given_gamma,given_scale,hicmindistance) for x in geneipeakcenterlist]
        
        #fillna
        geneiPeakxWeightList = np.nan_to_num(geneiPeakxWeightList)

        #t2 = time.time()

        #epigenomes times weight
        geneiWEscore = np.array(geneiPeakxWeightList)*np.array(geneiPeakxEpigenome)
        if logRgX:
            geneiWEscore=np.log1p(geneiWEscore)
        
        if noise_quantile > 0:
            quantile_threshold = np.percentile(geneiWEscore, noise_quantile * 100)
            geneiWEscore[geneiWEscore < quantile_threshold] = 0

        if noise_ratio > 0:
            total_sum = np.sum(geneiWEscore)
            if total_sum > 0:  
                ratio_threshold = noise_ratio * total_sum
                cumulative_sum = np.cumsum(np.sort(geneiWEscore))
                # 找到比例对应的阈值位置
                geneiWEscore[np.searchsorted(cumulative_sum, ratio_threshold)] = 0

        
        #这个是Rg的值
        #geneiWEscoreSum = np.array(geneiPeakxEpigenome).sum()
        geneiWEscoreSum = geneiWEscore.sum()

        
        #Rgx
        if geneiWEscoreSum > 0: #判断是否全为0
            geneiWEscore_percent = geneiWEscore / geneiWEscoreSum
        else:
            geneiWEscore_percent = geneiWEscore * 0

        if setpromoter1:
            geneiWEscore_percent[ifpromoter] = 1 #启动子设置为1

        genei_Rgx_df = geneipeaklist.iloc[:,0:4].copy()
        genei_Rgx_df.columns = ["peakChr",'peakStart','peakEnd','epigenomeActivity']
        genei_Rgx_df["geneSymbol"] = self.candidateGeneDF.iloc[i,3]
        genei_Rgx_df["geneChr"] = self.candidateGeneDF.iloc[i,0]
        genei_Rgx_df["geneStart"] = self.candidateGeneDF.iloc[i,1]
        genei_Rgx_df["geneEnd"] = self.candidateGeneDF.iloc[i,2]
        genei_Rgx_df["geneStrand"] = self.candidateGeneDF.iloc[i,5]
        genei_Rgx_df["geneID"] = self.candidateGeneDF.iloc[i,4]
        genei_Rgx_df["weight"] = list(geneiPeakxWeightList)
        genei_Rgx_df["Rgx_rawvalue"] = list(geneiWEscore)
        genei_Rgx_df["Rgx_percent"] = list(geneiWEscore_percent)

        #t3 = time.time()

        #print (t3-t1, t2-t1, (t2-t1)/(t3-t1))

        return(geneiWEscoreSum,genei_Rgx_df)
    
    def processGene(self,args):
        i, self, weightType, rpDecayDistance, fixedFunctionType, given_gamma, given_scale, \
        ref_gamma, ref_scale, hicmindistance, logRgX, setpromoter1,ifUseHiCRef,goldWeightDf = args
    
        return self.computeGenei(i, weightType, rpDecayDistance, fixedFunctionType, 
                                 given_gamma, given_scale, ref_gamma, ref_scale, 
                                 hicmindistance, logRgX, setpromoter1,ifUseHiCRef,goldWeightDf)
    
    def computeAllGene(self,weightType,halfDistance=10000,fixedFunctionType='exponential',
                 given_gamma=1.024238616787792, given_scale = 5.9594510043736655,
                 ref_gamma = 0.87, ref_scale = -4.80 + 11.63 * 0.87, hicmindistance=5000,
                 logRgX=False,setpromoter1=False,threads=1,ifUseHiCRef=False, 
                 userWeightFile=None,noise_ratio=0,noise_quantile=0):
        
        print("***Computing RgX and Rg")
        print("***Using "+weightType+" mode")
        RgxDfList = []
        RgList=[]
        numOfGene = self.candidateGeneDF.shape[0]

        if threads <= 1:
            # Wrap the loop with tqdm to show a progress bar
            for i in tqdm(range(numOfGene), desc="***Processing genes", unit="gene"):
                geneiWEscoreSum, genei_Rgx_df = self.computeGenei(
                    i, weightType, rpDecayDistance=halfDistance, 
                    fixedFunctionType=fixedFunctionType, given_gamma=given_gamma, 
                    given_scale=given_scale, ref_gamma=ref_gamma, ref_scale=ref_scale, 
                    hicmindistance=hicmindistance, logRgX=logRgX, setpromoter1=setpromoter1, 
                    ifUseHiCRef=ifUseHiCRef, goldWeightDf=userWeightFile,noise_ratio=noise_ratio,noise_quantile=noise_quantile
                )
                RgxDfList.append(genei_Rgx_df)
                RgList.append(geneiWEscoreSum)
        else:
            with ProcessPoolExecutor() as executor:
                futures = [executor.submit(self.computeGenei, i, weightType, rpDecayDistance=halfDistance,
                                            fixedFunctionType=fixedFunctionType, given_gamma=given_gamma,
                                            given_scale=given_scale, ref_gamma=ref_gamma, ref_scale=ref_scale, 
                                            hicmindistance=hicmindistance, logRgX=logRgX, setpromoter1=setpromoter1,
                                            ifUseHiCRef=ifUseHiCRef, goldWeightDf=userWeightFile) 
                        for i in range(numOfGene)]
                
                # 收集结果
                for future in futures:
                    geneiWEscoreSum, genei_Rgx_df = future.result()
                    RgxDfList.append(genei_Rgx_df)
                    RgList.append(geneiWEscoreSum)

            
        
        RgxDf = pd.concat(RgxDfList, ignore_index=True)
        # RgDF -> RgDf
        RgDf = self.candidateGeneDF.copy()
        # 加了header，和RgxDf保持一致
        RgDf.columns = ["geneChr",'geneStart','geneEnd','geneSymbol','geneID','geneStrand']
        RgDf["Rg"] = RgList

        self.RgxDf = RgxDf
        self.RgDf = RgDf

        print("***Finished")

    
    def clean(self,onlytmp=False):
        if os.path.exists(self.tmpdir):
                shutil.rmtree(self.tmpdir)

        if onlytmp:
            print("clean only tmp directory.")
        else:
            print("clean all attribute to release memory.")
            self.candidatesite_coverage = None
            self.nomhicdf = None
            self.RgxDf = None
            self.RgDf = None

        
        
    def save(self,outname="output",header=False):
        print("***Save RgX and Rg to tsv tables.")
        self.RgxDf.to_csv(outname+"_RgxDf.tsv.gz",header=header,sep="\t",index=None,compression="gzip")
        self.RgDf.to_csv(outname+"_RgDf.tsv.gz",header=header,sep="\t",index=None,compression="gzip")

        print("......Saved to "+outname+"_RgxDf.tsv.gz"+" and "+outname+"_RgDf.tsv.gz")

    
    


