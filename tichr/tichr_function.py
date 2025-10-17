import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, time
import tempfile
import subprocess
import pyBigWig

from .preprocess_hic import *



def makeSiteBedFunction(candidatesite,candidateGeneFile,readFileList,gtfile,
                        species='hs',binResolution=100,peakToGeneMaxDistance=100000,
                        blackregion=None, refgene_file=None,tmpdir='tmp_makeSiteBdg',fixPeakWidth=False):
    
    codepath = os.path.dirname(os.path.realpath(__file__))

    if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    print("Temporary directory created at: " + tmpdir)

    if not blackregion:
        with open(tmpdir+"/blackregion.blank.tmp", 'w') as file: pass  # 不写入任何内容
        blackregion=tmpdir+"/blackregion.blank.tmp"

    if not refgene_file:
        with open(tmpdir+"/refgene_file.blank.tmp", 'w') as file: pass  # 不写入任何内容
        refgene_file=tmpdir+"/refgene_file.blank.tmp"
    
    print("Step1. setting candidate sites")
    # candidatesite 的类型有三种：
    # 一是给定bed文件
    # 二是通过对bamfile call peak
    # 三是启动子周围的所有区间
    if candidatesite == "denovo_peak":
        print("Call peaks from the given bam files")
        subprocess.run(["bash", codepath + "/bashcode/ABC_make_candidate.sh",
                    species,tmpdir, " ".join(readFileList), 
                    gtfile, blackregion,refgene_file], stdout=open(tmpdir+"/denovo_peak.log", "w"), stderr=subprocess.STDOUT)
        candidatesite_file = tmpdir+"/candidatebed.enhancer.bed"
    elif candidatesite == 'surronding_bin':
        print("Generating bins for each "+str(binResolution)+"bp")

        subprocess.run(["bash", codepath + "/bashcode/makeSurrondingBin.sh", 
                        candidateGeneFile,str(peakToGeneMaxDistance),gtfile,str(binResolution),tmpdir], 
                        stdout=open(tmpdir+"/SurrondongBin.log", "w"))
        candidatesite_file = tmpdir+"/surrondingbin.bed"
    elif candidatesite == 'onlypromoter':
        print("Only use promoter")
        awk_command = f'''awk -v OFS="\\t" '$6 == "+" {{print $1,$2-100,$2+100}} $6 == "-" {{print $1,$3-100,$3+100}}' {candidateGeneFile} > {tmpdir}/onlyPromoter.bed'''
        subprocess.run(awk_command, shell=True, check=True)
        candidatesite_file = tmpdir+"/onlyPromoter.bed"
    elif os.path.exists(candidatesite):  # bed3 file 
        print("Use a user-defined candidate bed file")
        if fixPeakWidth:
            command = f"cat {candidatesite} | cut -f 1-3 | sortBed | mergeBed | " \
                      f"awk -v OFS='\\t' '{{mid=int(($2+$3)/2); print $1, mid-250, mid+250}}' | awk '$2>0'|sortBed | " \
                      f"intersectBed -v -wa -a stdin -b {blackregion} > {tmpdir}/candidatesite_final.bed"
            subprocess.run(command, shell=True, check=True)
        else:
            subprocess.run(f"intersectBed -v -wa -a {candidatesite} -b {blackregion} > {tmpdir}/candidatesite_final.bed" , shell=True, check=True)
        
        candidatesite_file = tmpdir+"/candidatesite_final.bed"

    else:
        print("please assign ['denovo_peak','surronding_bin' or the correct <filepath> for the candidate sites]")
        exit(1)
    
    return(candidatesite_file)


def makeSiteBdgFunction(candidatesite_file,readFileList,gtfile,coverageMethod,file_type,
                spmr = False,species='hs',refgene_file=None,tmpdir='tmp_makeSiteBdg',ifTSSrange=500,
                quantileref=None,signaltype=None,separatepromoter=False,quantile_method=None,multiBAMmerge="mean"):
    
    if not refgene_file:
        with open(tmpdir+"/refgene_file.blank.tmp", 'w') as file: pass  # 不写入任何内容
        refgene_file=tmpdir+"/refgene_file.blank.tmp"
    
    codepath = os.path.dirname(os.path.realpath(__file__))

    print("Step2. Extracting read counts from bam files")
    spmr_flag = 'yes' if spmr else "no"

    if spmr_flag == 'yes':
        print("use spmr")

    if file_type =="bam":
        if multiBAMmerge == 'mean':
            peakdf_values = []
            for read_bamfile in readFileList:
                print("Processing "+str(read_bamfile))
                #print(["bash", codepath + "/bashcode/coverageBed-tichr.sh", read_bamfile,
                #                gtfile,species,tmpdir,candidatesite_file,coverageMethod,spmr_flag,refgene_file,str(ifTSSrange)])
                
                subprocess.run(["bash", codepath + "/bashcode/coverageBed-tichr.sh", read_bamfile,
                                gtfile,species,tmpdir,candidatesite_file,coverageMethod,spmr_flag,refgene_file,str(ifTSSrange)],
                                stdout=open(tmpdir+"/ExtractRead.log", "w"), stderr=subprocess.STDOUT)
                candidatesite_coverage = pd.read_csv(tmpdir+"/candidateSiteCoverage.bdg",sep='\t',header=None)
                peakdf_values.append(list(candidatesite_coverage[3]))
                candidatesite_coverage[3] = np.array(peakdf_values).mean(axis=0)
        elif multiBAMmerge == 'sum':
            for read_bamfile in readFileList:
                print("Processing "+str(read_bamfile))
            subprocess.run(["bash", codepath + "/bashcode/coverageBed-tichr.sh", " ".join(readFileList),
                                gtfile,species,tmpdir,candidatesite_file,coverageMethod,spmr_flag,refgene_file,str(ifTSSrange)],
                                stdout=open(tmpdir+"/ExtractRead.log", "w"), stderr=subprocess.STDOUT)
            candidatesite_coverage = pd.read_csv(tmpdir+"/candidateSiteCoverage.bdg",sep='\t',header=None)

    elif file_type == "bw" or file_type == "bigWig" or file_type == "bigwig":  
        candidatesite_coverage = pd.read_csv(candidatesite_file, sep="\t", header=None, names=["chr", "start", "end"])
        candidatesite_coverage['value'] = 0
        spmr_weight = 1
        for read_signal_file in readFileList:
            with pyBigWig.open(read_signal_file) as bw:
                chromo_info = bw.chroms()
                if spmr:
                    total_coverage = 0
                    for chrom in bw.chroms().keys():
                        intervals = bw.intervals(chrom)
                        if intervals:
                            for start, end, value in intervals:
                                total_coverage += (end - start) * value
                    if total_coverage != 0:
                        spmr_weight = 1000000 / total_coverage
                    else:
                        print(f"Error: Can not read total coverage for file {read_signal_file}")
                        exit()
                for row in candidatesite_coverage.itertuples():
                    if row.end > chromo_info[row.chr]:
                        print(f"Warning: Candidate site [{row.chr}, {row.start}, {row.end}] is out of the range of the file [{read_signal_file}]")
                        candidatesite_coverage.at[row.Index, 'value'] = 0
                    else:
                        mean_value = bw.stats(row.chr, row.start, row.end, type="mean")[0]
                        if mean_value is None:
                            mean_value = 0
                        add_value = mean_value * (row.end - row.start) * spmr_weight
                        sum_value = row.value + add_value
                        candidatesite_coverage.at[row.Index, 'value'] = sum_value
        if multiBAMmerge == 'mean':
            candidatesite_coverage["value"] = candidatesite_coverage["value"] / len(readFileList)
        elif multiBAMmerge == "sum":
            pass
        else:
            print("Error: multiBAMmerge can only be specified as [sum] or [mean]")
            exit()
        candidatesite_coverage_without_TSScov_dir = f"{tmpdir}/candidatesite_coverage_without_TSScov.temp"
        candidatesite_coverage.to_csv(candidatesite_coverage_without_TSScov_dir, index=False, header=False, sep="\t")
        subprocess.run(["bash", codepath + "/bashcode/coverageTSS.sh", candidatesite_coverage_without_TSScov_dir,
                            tmpdir,refgene_file,str(ifTSSrange)],
                            stdout=open(tmpdir+"/ExtractRead.log", "w"), stderr=subprocess.STDOUT)
        candidatesite_coverage = pd.read_csv(tmpdir+"/candidateSiteCoverage.bdg",sep='\t',header=None)

    elif file_type == "bedGraph" or file_type == "bedgraph":
        if multiBAMmerge == 'mean':
            peakdf_values = []
            for read_signal_file in readFileList:
                print("Processing "+str(read_signal_file))
                subprocess.run(["bash", codepath + "/bashcode/coverageBed_bedGraph.sh", read_signal_file,
                                gtfile,species,tmpdir,candidatesite_file,coverageMethod,spmr_flag,refgene_file,str(ifTSSrange),codepath+"/bashcode"],
                                stdout=open(tmpdir+"/ExtractRead.log", "w"), stderr=subprocess.STDOUT)
                candidatesite_coverage = pd.read_csv(tmpdir+"/candidateSiteCoverage.bdg",sep='\t',header=None)
                peakdf_values.append(list(candidatesite_coverage[3]))
                candidatesite_coverage[3] = np.array(peakdf_values).mean(axis=0)
        elif multiBAMmerge == 'sum':
            for read_signal_file in readFileList:
                print("Processing "+str(read_signal_file))
            subprocess.run(["bash", codepath + "/bashcode/coverageBed_bedGraph.sh", " ".join(readFileList),
                                gtfile,species,tmpdir,candidatesite_file,coverageMethod,spmr_flag,refgene_file,str(ifTSSrange),codepath+"/bashcode"],
                                stdout=open(tmpdir+"/ExtractRead.log", "w"), stderr=subprocess.STDOUT)
            candidatesite_coverage = pd.read_csv(tmpdir+"/candidateSiteCoverage.bdg",sep='\t',header=None)
        else:
            print("Error: multiBAMmerge can only be specified as [sum] or [mean]")
            exit()

    print("Step3. If execute quantile normalization")
    if quantileref and not signaltype:
        print("Please specify singaltype as DNase or H3K27ac")
        exit(1)
    elif quantileref and signaltype and quantile_method:
        print("Using reference file to make quantiles of ChIP-seq signals")
        refdf = pd.read_csv(quantileref, sep='\t')
        if separatepromoter:
            print("separate promoter")
            candidatesite_coverage = quantile_normalize_separate(candidatesite_coverage,refdf,signaltype,method=quantile_method)
        else:
            print("not separate promoter")
            candidatesite_coverage = quantile_normalize_any(candidatesite_coverage,refdf,signaltype,method=quantile_method)
    else:
        candidatesite_coverage = candidatesite_coverage

    candidatesite_coverage[0] = pd.Categorical(candidatesite_coverage[0])
        
    return(candidatesite_coverage)

def getpowerlaw(distance,gamma,scale,mindistance):
    distance=abs(distance)
    if not mindistance: mindistance = 5000
    if distance<mindistance:
        distance=mindistance
    powerlaw = math.exp(scale + -1 * gamma * math.log(distance))
    return(powerlaw)

def weight3Dgenome(structureDataType,structureDataFile,peakPos,tssPos):
    # structureDataType in ["tad","loop","compartment"]
    pass

def makeWeightFunction(weightType,peakPos,tssPos,rpDecayDistance=10000,fixedFunctionType='exponential',
                 given_gamma=1.024238616787792, given_scale = 5.9594510043736655,
                 ref_gamma = 0.87, ref_scale = -4.80 + 11.63 * 0.87, hicmindistance=5000,
                 hicProcessedData=None,hicRes=50000,geneChr=None,hicProcessedDataType=None,ifUseHiCRef=False,
                 peakToGeneMaxDistance=None, goldWeightDf=None
                 ): 
    if weightType == 'fixed_function':
        #print("using fixed function as the weight for epigenome")
        z = abs(peakPos-tssPos)
        if fixedFunctionType == 'sigmoid':
            lamda = -math.log(1/3)*1e5/rpDecayDistance
            weight = 2*math.exp(-lamda*z/1e5)/(1.0+math.exp(-lamda*z/1e5)) 
        elif fixedFunctionType == 'exponential':
            weight = 2 ** -(z/rpDecayDistance)
        elif fixedFunctionType == 'powerlaw':
            lamda = math.log(0.5)/math.log(rpDecayDistance)
            weight = (z+0.1)**lamda
        elif fixedFunctionType == 'normPowerLaw':
            weight = getpowerlaw(z,given_gamma,given_scale,hicmindistance)
        elif fixedFunctionType == 'constant1':
            weight = 1
        elif fixedFunctionType == 'linear-half':
            weight = max(0, 1 - 0.5 * (z / rpDecayDistance))
        elif fixedFunctionType == 'linear':
            weight = max(0, 1 - (z / peakToGeneMaxDistance))
        elif fixedFunctionType == 'closest':
            weight = abs(peakPos-tssPos)
        else:
            print("please give a correct fixedFunctionType")
        return(weight)
    elif weightType == 'hic' and hicProcessedData and geneChr and hicProcessedDataType:
        tssHicIndex = int (tssPos / hicRes)
        peakHiCIndex = int (peakPos / hicRes)

        if hicProcessedDataType == "rawhic_dense" or hicProcessedDataType == "matrix_dense":
            hicInteractionValue = hicProcessedData[geneChr][tssHicIndex,peakHiCIndex]
        elif hicProcessedDataType == 'rawhic_sparse':
            if tssHicIndex <= peakHiCIndex:
                hicInteractionName = str(tssHicIndex*hicRes)+"to"+str(peakHiCIndex * hicRes)
            else:
                hicInteractionName = str(peakHiCIndex * hicRes)+"to"+str(tssHicIndex*hicRes)
            
            try:
                hicInteractionValue = hicProcessedData[geneChr].loc[hicInteractionName]["counts"]
            except:
                hicInteractionValue = 0

        if ifUseHiCRef:
            pseudo_contact = getpowerlaw(abs(peakPos-tssPos),given_gamma,given_scale,hicRes)
            plref_contact = getpowerlaw(abs(peakPos-tssPos),ref_gamma,ref_scale,hicmindistance)
            plfit_contact = getpowerlaw(abs(peakPos-tssPos),given_gamma,given_scale,hicmindistance)
            hicInteractionValue = hicInteractionValue * (plref_contact/plfit_contact) + pseudo_contact

        return(hicInteractionValue)
    
    elif weightType == 'gold_weight':
        """用户指定weight文件, 格式为 chr1,start1,end1,chr2,start2,end2,weight 七列"""
        if goldWeightDf is None:
            print("Please give the goldWeightDf for gold_weight")
            exit(1)
        else:
            df = pd.read_csv(goldWeightDf, sep='\t', header=None)
            filtered_df = df[df[0] == geneChr]

            # print(f"####!!!!geneChr is {geneChr}"
            #       f"peakPos is {peakPos}"
            #       f"tssPos is {tssPos}")

            for _, row in filtered_df.iterrows():
                if row[1] < peakPos < row[2] and row[4] < tssPos < row[5]:
                    return row[6] 
            
            for _, row in filtered_df.iterrows():
                if row[1] < tssPos < row[2] and row[4] < peakPos < row[5]:
                    return row[6] 
                
            return 0


    else:
        print("Please specify the weightType to [hic, fixed_function, gold_weight]")
        exit(1)

        