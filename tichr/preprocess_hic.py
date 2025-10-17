import pandas as pd
import numpy as np
from scipy import interpolate
from multiprocessing import Pool
import scipy.stats as stats
import math
from functools import partial
from sys import exit
import hicstraw
from collections import defaultdict
import threading
import concurrent.futures

from .dumphic import *



def oeNormalizeSparse(records,outType="OE"):
    distance_sums = defaultdict(float)
    distance_counts = defaultdict(float)
    
    for binX, binY, value in records:
        d = abs(binX-binY)
        distance_sums[d] += value
        distance_counts[d] += 1
        
    average_distances = {d: distance_sums[d] / distance_counts[d] for d in distance_sums}
    
    if outType=="OE":
        out_records = [
            [binX, binY, value / average_distances[abs(binX - binY)]]
            for binX, binY, value in records
            ]
    elif outType=="expected":
        out_records = [
            [binX, binY, average_distances[abs(binX - binY)]]
            for binX, binY, value in records
            ]
    df = pd.DataFrame(out_records, columns=["posX", "posY", "counts"])
    df.index = df.posX.astype(str)+"to"+df.posY.astype(str)
    return(df)

def noFurtherNormalizeSparse(records):
    df = pd.DataFrame(records, columns=["posX", "posY", "counts"])
    df.index = df.posX.astype(str)+"to"+df.posY.astype(str)
    return(df)

def scale0to1Sparse(records):
    df = pd.DataFrame(records, columns=["posX", "posY", "counts"])
    df["counts"] = df["counts"] / df["counts"].quantile(0.95) #本来想用max的，但考虑到有极大值。
    df.index = df.posX.astype(str)+"to"+df.posY.astype(str)
    return(df)

def totalNormSparse(records):
    df = pd.DataFrame(records, columns=["posX", "posY", "counts"])
    df["counts"] = (10000000 * df["counts"]) / np.nansum(df["counts"]) 
    df.index = df.posX.astype(str)+"to"+df.posY.astype(str)
    return(df)

def add_records_to_bin_sums(records, bin_sums, start, end):
    for binX, binY, value in records:
        if binX < start or binY > end:
            value /= 2
        if binX == binY:
            bin_sums[binX] += value
        else:
            bin_sums[binX] += value
            bin_sums[binY] += value

def normDF_sparse(records,hicres,chri_len):
    bin_sums = defaultdict(float)
    add_records_to_bin_sums(records, bin_sums,0,chri_len)
    row_mean = np.mean(list(bin_sums.values()))

    df = pd.DataFrame(records, columns=["posX", "posY", "counts"])
    df.index = df.posX.astype(str)+"to"+df.posY.astype(str)

    #divide row sums
    df['counts'] = df['counts']/row_mean

    #correct diagnol
    diagbool = df["posX"] == df["posY"]
    for diagpos in df[diagbool]["posX"]:
        left_count, right_count = 0, 0
        diagleft_posX =  diagpos - hicres
        diagleft_posY = diagpos
        try:
            diagleft_value = df.loc[str(diagleft_posX)+"to"+str(diagleft_posY),"counts"]
            left_count = diagleft_value
        except:
            left_count = 0
        diagright_posX = diagpos
        diagright_posY = diagpos + hicres
        try:
            diagright_value = df.loc[str(diagright_posX)+"to"+str(diagright_posY),"counts"]
            right_count = diagright_value
        except:
            right_count = 0
            
        maxneighbor = max(left_count,right_count)
        if maxneighbor:
            df.loc[str(diagpos)+"to"+str(diagpos),"counts"] = maxneighbor
    
    return(df)


def paral_sparse(hicfilepath,chri,hicnorm,hicres,gtdf,further_normalize_type='abc'):
    chri_len = int(gtdf[gtdf[0]==chri][1])
    hic = hicstraw.HiCFile(hicfilepath)

    try:
        mtobt= hic.getMatrixZoomData(chri, chri, "observed", hicnorm, "BP", hicres)
        mt = mtobt.getRecords(0,chri_len,0,chri_len)
    except:
        print("try to use 1,2... instead of chr1,chr2...  ")
        chri_num = chri.replace('chr', '')
        mtobt= hic.getMatrixZoomData(chri_num, chri_num, "observed", hicnorm, "BP", hicres)
        mt = mtobt.getRecords(0,chri_len,0,chri_len)

    if len(mt)==0:
        print("(null mt) try to use 1,2... instead of chr1,chr2...  ")
        chri_num = chri.replace('chr', '')
        mtobt= hic.getMatrixZoomData(chri_num, chri_num, "observed", hicnorm, "BP", hicres)
        mt = mtobt.getRecords(0,chri_len,0,chri_len)
    
    records = [[r.binX, r.binY, r.counts] for r in mt]
    if further_normalize_type == 'abc':
        nomhicchri = normDF_sparse(records,hicres,chri_len)
    elif further_normalize_type == 'oe':
        nomhicchri = oeNormalizeSparse(records,outType="OE")
    elif further_normalize_type == 'no_further':
        nomhicchri = noFurtherNormalizeSparse(records)
    elif further_normalize_type == '0to1': #如果把Hi-C根据最大值或总和标准化的前提是不同样本之间总的接触强度是一个恒定值，可用去区分局部基因组位点之间的差异。
                                            # 但如果多个样本本身的contact总和不是恒定的，比如一个样本完全没有接触，另一个样本有很多接触，这样标准化其实不够公平。
        nomhicchri = scale0to1Sparse(records)
    elif further_normalize_type == "total":
        nomhicchri = totalNormSparse(records)
    else:
        print("please use one of ['abc','oe','no_further','0to1','total']")

    return chri,nomhicchri

def gethicfile(hicfilepath,hicres,hictype,genechrlist,
               hicnorm="SCALE",gt=None,juicertool=None,threads=8,further_normalize_type='abc'):
    
    if hictype=="matrix_dense":
        for chri in genechrlist:
            nomhicdf[chri] = np.array(normalizehic_dense(hicfilepath))
    elif hictype == 'rawhic_sparse':
        gtdf = pd.read_csv(gt,sep="\t",header=None)
        with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
            nomhicdf={}
            futures = {executor.submit(paral_sparse,hicfilepath,chri,hicnorm,hicres,gtdf,further_normalize_type): 
                    chri for chri in genechrlist}
            for future in concurrent.futures.as_completed(futures):
                chri, nomhicchri = future.result()
                nomhicdf[chri] = nomhicchri
        '''
        for chri in genechrlist:
            print(chri)
            mtobt= hic.getMatrixZoomData(chri, chri, "observed", hicnorm, "BP", hicres)
            chri_len = int(gtdf[gtdf[0]==chri][1])
            mt = mtobt.getRecords(0,chri_len,0,chri_len)
            records = [[r.binX, r.binY, r.counts] for r in mt]
            nomhicdf[chri] = normDF_sparse(records,hicres,chri_len)
        '''
    elif hictype=='rawhic_dense':
        manyJuicer(hicfilepath,hicnorm,hicres,gt,"juicerdump",juicertool,genechrlist,threads=threads)
        for chri in genechrlist:
            #genedfchri = self.genedf[self.genedf[0]==chri]
            #tsslist = genedfchri.apply(lambda row: row[1] if row[5] == "+" else row[2], axis=1)
            matrixfile = "juicerdump/"+str(hicres)+"/observed."+str(hicnorm)+"."+str(chri)+".matrix.gz"
            nomhicdf[chri]= np.array(normalizehic_dense(matrixfile))
    else:
        print("Please give the hictype parameter ['matrix','rawhic']")

    return(nomhicdf)

def multicpu_eachindex(func,rangeindex,threads): #多进程
    pool = Pool(threads)
    # 使用 map 函数并行处理任务，并获取结果
    abcdflist=pool.map(func,rangeindex)
    # 关闭进程池
    pool.close()
    pool.join()
    return(abcdflist)

def calculate_diff_and_values(hicdf):
    values_dict = {}    
    binnumber = len(hicdf)
    def process_element(i, j):
        diff = abs(i - j)
        value = hicdf[i, j]
        # 将距离差值和元素值添加到字典中
        values_dict.setdefault(diff, []).append(value)
    
    for i in range(binnumber):
        for j in range(i + 1, binnumber):
            process_element(i, j)
    
    # 创建新的 DataFrame，列出所有可能的差值和对应的元素值
    result_df = pd.DataFrame(columns=["dist_for_fit", "hic_contact"])
    for diff, values in values_dict.items():
        result_df = pd.concat([result_df,pd.DataFrame({"dist_for_fit": [diff], "hic_contact": np.mean(values)})],ignore_index=True)
    
    return result_df

def makemt(i,lenhic,res,mindistance,gamma,scale,ref_gamma,ref_scale):
    eachrow_fitdf = []
    eachrow_refdf = []
    eachrow_pseudodf = []

    for j in range(lenhic):
        distance_ij = abs(i-j) * res
        if distance_ij < mindistance:
            distance_ij = mindistance 
            
        fitcount = getpowerlaw(distance_ij,gamma,scale,mindistance)
        pseudocount = getpowerlaw(mindistance,gamma,scale,mindistance)
        refcount = getpowerlaw(distance_ij,ref_gamma,ref_scale,mindistance)

        eachrow_fitdf.append(fitcount)
        eachrow_refdf.append(refcount)
        eachrow_pseudodf.append(min(pseudocount,fitcount))

    return((eachrow_fitdf,eachrow_refdf,eachrow_pseudodf))

def normalizehic_dense(matrixfile,smoothdiagbin=1):
    print("Read input files...")
    hicrawdf = pd.read_csv(matrixfile,delimiter='\t', index_col=0)
    hicdf = np.array(hicrawdf)
    lenhic = len(hicdf)

    print("Normalize to rowsum = 1")
    rowsum =np.nansum(hicdf,axis=1)
    #summean= np.nanmean(rowsum)
    rowsum_mean = np.nansum(rowsum)/lenhic
    #if abs(rowsum_mean-1) > 0.001:
    hicdf = hicdf/rowsum_mean

    print("Adjust diagnal...")
    for i in range(lenhic):
        lowerindex = i-smoothdiagbin
        upperindex = i+smoothdiagbin
        if lowerindex<0:
            lowerindex=0
        if upperindex >= lenhic:
            upperindex=lenhic-1
        hicdf[i,i]= max(hicdf[i, lowerindex], hicdf[i, upperindex])

    normalized_df = pd.DataFrame(hicdf)
    normalized_df.columns = hicrawdf.columns
    normalized_df.index = hicrawdf.index
    return(normalized_df)

def normalizehic_old(matrixfile,res,given_gamma,given_scale,ref_gamma,ref_scale,tsslist,hicmindistance,smoothbin=1,threads=32):
    print("Read input files...")
    
    hicrawdf = pd.read_csv(matrixfile,delimiter='\t', index_col=0)
    hicdf = np.array(hicrawdf)
    lenhic = len(hicdf)

    print("Normalize to rowsum = 1")
    rowsum =np.nansum(hicdf,axis=1)
    #summean= np.nanmean(rowsum)
    rowsum_mean = np.nansum(rowsum)/lenhic
    #if abs(rowsum_mean-1) > 0.001:
    hicdf = hicdf/rowsum_mean

    
    print("Adjust diagnal...")
    for i in range(lenhic):
        lowerindex = i-smoothbin
        upperindex = i+smoothbin
        if lowerindex<0:
            lowerindex=0
        if upperindex >= lenhic:
            upperindex=lenhic-1
        hicdf[i,i]= max(hicdf[i, lowerindex], hicdf[i, upperindex])



    # 应该不需要fit
    '''
    print("Fit the powerlaw...") 
    dis_and_contact = calculate_diff_and_values(np.array(hicdf))
    pseudocount = 0.000001
    logdist = np.log(dis_and_contact["dist_for_fit"].astype(int)*res + pseudocount)
    logcontact = np.log(dis_and_contact["hic_contact"] + pseudocount)
    reg = stats.linregress(logdist,logcontact)
    gamma = -reg.slope
    scale = reg.intercept
    '''
    gamma = given_gamma
    scale = given_scale

    print("Normalize by powerlaw...")

    '''
    partial_makemt = partial(makemt,lenhic=lenhic,res=res,mindistance=hicmindistance,gamma=gamma,scale=scale,ref_gamma=ref_gamma,ref_scale=ref_scale)
    fitdf_refdf = np.array(multicpu_eachindex(partial_makemt,range(lenhic),threads))
    fitdf = fitdf_refdf[:, 0, :]
    refdf = fitdf_refdf[:, 1, :]
    pseudodf = fitdf_refdf[:, 2, :]
    #这些都不用

    # Gene promoters with insufficient hic coverage should get replaced with powerlaw

    tsslist_hicindex = np.array([int(i/res) for i in tsslist])
    promoter_threshold = 0.01
    badgenebool = np.nansum(hicdf[tsslist_hicindex],axis=1) < promoter_threshold
    badgeneindex = tsslist_hicindex[badgenebool]
    hicdf[badgeneindex,:] = fitdf[badgeneindex,:]
    hicdf[:,badgeneindex] = fitdf[:,badgeneindex]
    '''

    hicdf = np.nan_to_num(hicdf)

    '''
    for i in range(lenhic):
        for j in range(lenhic):
            distance_ij = abs(i-j) * res
            if distance_ij < res * smoothbin:
                distance_ij = res * smoothbin #1*5000
                
            fitdf[i,j] = math.exp(scale + -1 * gamma * np.log(distance_ij))
            refdf[i,j] = math.exp(ref_scale + -1 * ref_gamma * np.log(distance_ij))
    '''

    normalized_df = pd.DataFrame(hicdf)
    normalized_df.columns = hicrawdf.columns
    normalized_df.index = hicrawdf.index

    return(normalized_df)

def getpowerlaw(distance,gamma,scale,mindistance):
    distance=abs(distance)
    if not mindistance: mindistance = 5000
    if distance<mindistance:
        distance=mindistance
    powerlaw = math.exp(scale + -1 * gamma * math.log(distance))
    return(powerlaw)

def quantile_normalize_separate(dfinput,refdf,signaltype,column=3,method="rankpercent"):
    df = dfinput.copy()
    promoterindex=df.index[df[4]>0]
    nonpromoterindex=df.index[df[4]==0]

    if signaltype=="H3K27ac":
        signaltype = "H3K27ac.RPM"
    elif signaltype=="DNase":
        signaltype = "DHS.RPM"
    else:
        print("improper signaltype")
        exit(1)

    ref_raw = refdf
    ref_promoter = ref_raw[ref_raw["enh_class"]=="promoter"]
    ref_nonpromoter = ref_raw[ref_raw["enh_class"]=="nonpromoter"]

    if method == 'rank':
        interpfunc_promoter = interpolate.interp1d(ref_promoter["rank"],ref_promoter[signaltype],
                                kind="linear",fill_value="extrapolate",)
        interpfunc_nonpromoter = interpolate.interp1d(ref_nonpromoter["rank"],ref_nonpromoter[signaltype],
                                kind="linear",fill_value="extrapolate",)
        df.loc[promoterindex, "quantile"] = len(promoterindex)- df.loc[promoterindex, column].rank()
        df.loc[promoterindex, column] = interpfunc_promoter(df.loc[promoterindex,"quantile"]).clip(0)
        df.loc[nonpromoterindex, "quantile"] = len(nonpromoterindex) - df.loc[nonpromoterindex, column].rank() 
        df.loc[nonpromoterindex, column] = interpfunc_nonpromoter(df.loc[nonpromoterindex,"quantile"]).clip(0)
    else:
        interpfunc_promoter = interpolate.interp1d(ref_promoter["quantile"],ref_promoter[signaltype],
                                kind="linear",fill_value="extrapolate",)
        interpfunc_nonpromoter = interpolate.interp1d(ref_nonpromoter["quantile"],ref_nonpromoter[signaltype],
                                kind="linear",fill_value="extrapolate",)
        df.loc[promoterindex, "quantile"] = df.loc[promoterindex, column].rank() / len(promoterindex)
        df.loc[promoterindex, column] = interpfunc_promoter(df.loc[promoterindex,"quantile"]).clip(0)
        df.loc[nonpromoterindex, "quantile"] = df.loc[nonpromoterindex, column].rank() / len(nonpromoterindex)
        df.loc[nonpromoterindex, column] = interpfunc_nonpromoter(df.loc[nonpromoterindex,"quantile"]).clip(0)
    outdf = df.iloc[:,0:5]
    return(outdf)


def quantile_normalize_any(dfinput,refdf,signaltype,column=3,method="rankpercent"):
    df = dfinput.copy()

    if signaltype=="H3K27ac":
        signaltype = "H3K27ac.RPM"
    elif signaltype=="DNase":
        signaltype = "DHS.RPM"
    else:
        print("improper signaltype")
        exit(1)

    ref_raw = refdf
    ref_any= ref_raw[ref_raw["enh_class"]=="any"]

    if method=='rankpercent':
        interpfunc_any = interpolate.interp1d(ref_any["quantile"],ref_any[signaltype],
                                kind="linear",fill_value="extrapolate",) #拟合
        df.loc[:, "quantile"] = df.loc[:, column].rank() / len(df)
        df.loc[:, column] = interpfunc_any(df.loc[:,"quantile"]).clip(0) #应用
    elif method == 'rank':
        interpfunc_any = interpolate.interp1d(ref_any["rank"],ref_any[signaltype],
                                kind="linear",fill_value="extrapolate",) #拟合
        #df.loc[:, "quantile"] = (1- (df.loc[:, column].rank() / len(df)))*len(df)
        df.loc[:, "quantile"] = len(df)- df.loc[:, column].rank()
        df.loc[:, column] = interpfunc_any(df.loc[:,"quantile"]).clip(0) #应用

    outdf = df.iloc[:,0:5]
    return(outdf)
