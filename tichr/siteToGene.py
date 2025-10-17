from joblib import Parallel, delayed
from multiprocessing import Pool, cpu_count
import os
import warnings
warnings.filterwarnings('ignore')


from .PRC_ROC import *
from .highOrderStructure import *
from .context import *

def matchgold(rgxfile,golddf,outname,goldcol,percent=False,withhead=False,returnDF=False,predGeneCol=10,predScoreCol=12):
    codepath = os.path.dirname(os.path.realpath(__file__))
    print(codepath)

    if percent:
        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                  rgxfile+ " "+str(goldcol)+" 10 13 "+outname+" "+str(withhead))
    else:
        print("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                  rgxfile+ " "+str(goldcol)+" "+str(predGeneCol)+ " " +str(predScoreCol)+ " " +outname+" "+str(withhead))

        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                  rgxfile+ " "+str(goldcol)+" "+str(predGeneCol)+ " " +str(predScoreCol)+ " " +outname+" "+str(withhead))
        
    
    if returnDF:
        if not withhead:
            matched = pd.read_csv(outname,sep="\t",header=None)
        else:
            matched = pd.read_csv(outname,sep="\t")
        return(matched)


def plotMultiPrc(indir,matchcol,truecol,outname="prc.pdf",dataset="K562bengiCRISPR",
                typelist=["noadj","adjRaw","adjRank","adjRankRaw","adjStruc",
                          "adjRawStruc","adjRankStruc","adjRankStrucRaw"],matchedwithhead=False,cols=cols):
    plt.figure(figsize=(4,4))
    i=1
    for name in typelist:
        if not matchedwithhead:
            matched = pd.read_csv(indir+name+".tsv",sep="\t",header=None)
        else:
            matched = pd.read_csv(indir+name+".tsv",sep="\t")
        score = matched[matchcol]
        truelabel = matched[truecol] == True
        pltoneprc(truelabel,score,name,cols[i])
        i = i +1
        
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(dataset+" dataset")
    plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    # plt.tight_layout()
    plt.savefig(outname)

def findStrucWeight(rgxfile_raw,rgfile_raw,golddf,goldcol,outdir,
                    structureTypeList,structureFileList,structureWeightList,
                    matchcol=12+1,truecol=11+1,onlyPlot=False,goldwithhead=False,
                    typelist=["noadj","adjRaw","adjStruc","adjRawStruc",]):
    # 所有列数为真实列数，不是python里面的0-based
    if not os.path.exists(outdir): os.makedirs(outdir)
    print("Try raw adjustment")
    matchgold(rgxfile_raw,golddf,outdir+"noadj.tsv",goldcol,percent=True,withhead=goldwithhead)
    matchgold(rgxfile_raw,golddf,outdir+"adjRaw.tsv",goldcol,percent=False,withhead=goldwithhead)

    print("Try structure adjustment")
    rg_adjstruc, rgx_adjstruc = adjStructure(rgxfile_raw,
                                                structureFileList, structureTypeList, structureWeightList,
                                            rgfile=rgfile_raw,rggeneID=4,rgvalue=6,
                                                cpumode="multi",withcolname=True)

    rgx_adjstruc.to_csv(outdir+"adjStruc_beforematch.tsv",sep="\t",header=None,index=False)

    matchgold(outdir+"adjStruc_beforematch.tsv",golddf,outdir+"adjStruc.tsv",
            goldcol,percent=True,withhead=goldwithhead)
    matchgold(outdir+"adjStruc_beforematch.tsv",golddf,outdir+"adjRawStruc.tsv",
            goldcol,percent=False,withhead=goldwithhead)
    
    plotMultiPrc(outdir,matchcol=matchcol-1,truecol=truecol-1,dataset=outdir,matchedwithhead=goldwithhead,
                typelist=typelist,outname="find_structure_weight.pdf")


def allcombination(rgxfile_raw,rgfile_raw,golddf,goldcol,outdir,tpmfile,
                   structureTypeList,structureFileList,structureWeightList,
                   tmpcolrep=[2+1,3+1,],ignorehead=True,tmpgeneID=1+1,
                    matchcol=12+1,truecol=11+1,onlyPlot=False,goldwithhead=False,ranktype="diffrank",
                    typelist=["noadj","adjRaw","adjRank","adjStruc","adjRankRaw",
                          "adjRawStruc","adjStrucRank","adjStrucRankRaw"]):
    # 所有列数为真实列数，不是python里面的0-based
    if not os.path.exists(outdir): os.makedirs(outdir)

    if not onlyPlot:
        print("Try raw adjustment")
        matchgold(rgxfile_raw,golddf,outdir+"noadj.tsv",goldcol,percent=True,withhead=goldwithhead)
        matchgold(rgxfile_raw,golddf,outdir+"adjRaw.tsv",goldcol,percent=False,withhead=goldwithhead)

        print("Try rank adjustment")
        rg_adjrank, rgx_adjrank= adjtpm(rgxfile_raw,tpmfile,rgfile=rgfile_raw,rggeneID=4,rgxgeneID=9,
                                       rgxvalue=11,rgxratio=12,
                                       tmpcolrep=np.array(tmpcolrep)-1,ignorehead=ignorehead,tmpgeneID=tmpgeneID-1,
                                       ranktype=ranktype)
        rgx_adjrank.to_csv(outdir+"adjRank_beforematch.tsv",sep="\t",header=None,index=False)
        matchgold(outdir+"adjRank_beforematch.tsv",golddf,outdir+"adjRank.tsv",
                  goldcol,percent=True,withhead=goldwithhead)
        matchgold(outdir+"adjRank_beforematch.tsv",golddf,outdir+"adjRankRaw.tsv",
                  goldcol,percent=False,withhead=goldwithhead)

        print("Try structure adjustment")
        rg_adjstruc, rgx_adjstruc = adjStructure(rgxfile_raw,
                                                 structureFileList, structureTypeList, structureWeightList,
                                                rgfile=rgfile_raw,rggeneID=4,rgvalue=6,
                                                 cpumode="multi",withcolname=True)

        rgx_adjstruc.to_csv(outdir+"adjStruc_beforematch.tsv",sep="\t",header=None,index=False)

        matchgold(outdir+"adjStruc_beforematch.tsv",golddf,outdir+"adjStruc.tsv",
              goldcol,percent=True,withhead=goldwithhead)
        matchgold(outdir+"adjStruc_beforematch.tsv",golddf,outdir+"adjRawStruc.tsv",
              goldcol,percent=False,withhead=goldwithhead)

        print("Try all adjustment")
        rg_adjStrucRank, rgx_adjStrucRank = adjtpm(outdir+"adjStruc_beforematch.tsv",tpmfile,
                                    rgfile=rgfile_raw,rggeneID=4,rgxgeneID=9,ranktype=ranktype,
                                    rgxvalue=11,rgxratio=12,
                                    tmpcolrep=np.array(tmpcolrep)-1,ignorehead=ignorehead,tmpgeneID=tmpgeneID-1)
        rgx_adjStrucRank.to_csv(outdir+"adjStrucRank_beforematch.tsv",sep="\t",header=None,index=False)

        matchgold(outdir+"adjStrucRank_beforematch.tsv",golddf,outdir+"adjStrucRank.tsv",
          goldcol,percent=True,withhead=goldwithhead)

        matchgold(outdir+"adjStrucRank_beforematch.tsv",golddf,outdir+"adjStrucRankRaw.tsv",
              goldcol,percent=False,withhead=goldwithhead)

    plotMultiPrc(outdir,matchcol=matchcol-1,truecol=truecol-1,dataset=outdir,matchedwithhead=goldwithhead,
                typelist=typelist,outname="all_combination.pdf")



def cut_distance(df,rgdf,maxdis=100000): #df是RgxDf
    peakpos = (df[1] + df[2]) / 2
    tsspos = np.where(df.iloc[:, 8] == "+", df.iloc[:, 6], df.iloc[:, 7])
    peak2tss = abs(peakpos-tsspos)
    cutdf = df[peak2tss< maxdis]
    cutdf.reset_index(drop=True, inplace=True)
    
    rgdf.index = rgdf[3]
    rgdf[6] = cutdf.groupby(4)[11].sum() #这个不是基于deg文件的rgdf
    cutrgdf= rgdf.dropna(subset=[rgdf.columns[6]])
    cutrgdf.reset_index(drop=True, inplace=True)
    
    return(cutdf,cutrgdf)

def adjtpm(rgxfile,tpmfile,
           rgfile=None,rgvalue=6,rggeneID=4,
           rgxgeneID=9,rgxvalue=11,rgxratio=12,
           tmpgeneID=1,tmpcolrep=[2,],ignorehead=False,
           ranktype="diffrank"):
    tpmdf=pd.read_csv(tpmfile,sep="\t")
    tpmdf["meanTPM"] = tpmdf.iloc[:,tmpcolrep].mean(axis=1)
    geneTPMdict = tpmdf.set_index(tpmdf.columns[tmpgeneID])["meanTPM"].to_dict()
    
    if ignorehead:
        rgxraw = pd.read_csv(rgxfile,sep="\t",header=None,skiprows=1)
    else:
        rgxraw = pd.read_csv(rgxfile,sep="\t",header=None)
    rgxdf = rgxraw.copy()
    
    rgxdf["geneID"] = rgxdf[rgxgeneID].apply(lambda x: x.split('.')[0])
    rgxdf["geneTPM"] = rgxdf["geneID"].map(geneTPMdict).fillna(0)
    
    rgxscore = rgxdf[rgxvalue]
    rgxpercent = rgxdf[rgxratio]

    if ranktype=="diffrank":
        rgxadj = rgxscore * (1-abs(makediffrank(rgxscore,rgxdf["geneTPM"])))
        rgxpercentadj = rgxpercent * (1-abs(makediffrank(rgxpercent,rgxdf["geneTPM"])))
    elif ranktype=="sumrank":
        rgxadj = rgxscore * (makesumrank(rgxscore,rgxdf["geneTPM"]))
        rgxpercentadj = rgxpercent * (makesumrank(rgxpercent,rgxdf["geneTPM"]))

    rgxraw[rgxvalue]=rgxadj
    rgxraw[rgxratio]=rgxpercentadj
    
    if rgfile: #只有在rgxvalue指定为raw的时候有效
        rgdf = pd.read_csv(rgfile,sep="\t",header=None)
        rgdf.index = rgdf[rggeneID]
        rgadj = rgxraw.groupby(rgxgeneID)[rgxvalue].sum()
        rgdf[rgvalue] = rgadj

    print("finish the adjustment for Rg and RgX")
    return(rgdf,rgxraw)


def compute_structure_weight(args):
    row, structureType, structureDF, structureWeight = args
    chrnum = row["peakChr"]
    peakPos = (row["peakStart"] + row["peakEnd"]) / 2
    tssPos = row["geneStart"] if row["geneStrand"] == "+" else row["geneEnd"]
    return weightStructureFunc(structureType, structureDF, structureWeight, peakPos, tssPos, chrnum)

def adjStructure(rgxfile,structureFile, structureType, structureWeight,ignorehead=False,
                 rgfile=None,rggeneID=4,rgvalue=6,cpumode="multi",outname="structureAdjusted",withcolname=True):
    global structureDF
    
    if all(isinstance(var, list) for var in [structureFile, structureType, structureWeight]):
        print("structureDF, structureType and structureWeight are all list")
        structureDF = []
        for file in structureFile:
            df = pd.read_csv(file, header=None, sep="\t")
            df[0] = pd.Categorical(df[0]) 
            structureDF.append(df)
    elif all(isinstance(var, str) for var in [structureFile, structureType]):
        print("structureDF, structureType and structureWeight are all str")
        structureDF = pd.read_csv(structureFile,header=None,sep="\t")
        structureDF[0] = pd.Categorical(structureDF[0]) #make chromosome name to category
    else:
        print("structureDF, structureType and structureWeight are not consistent style")
        exit(1)
    
    rgxdf=pd.read_csv(rgxfile,sep="\t")
    if not withcolname:
        rgxdf.columns = ['peakChr', 'peakStart', 'peakEnd', 'epigenomeActivity',
       'geneID', 'geneChr', 'geneStart', 'geneEnd', 'geneStrand',
       'geneSymbol', 'weight', 'Rgx_rawvalue', 'Rgx_percent']
    
    # 计算 structure_weight_list
    if cpumode == "single":
        structure_weight_list = [compute_structure_weight((row, structureType, structureDF, structureWeight)) 
                                 for _, row in rgxdf.iterrows()]
    
    elif cpumode == "joblib":
        structure_weight_list = Parallel(n_jobs=48)(
            delayed(compute_structure_weight)((row, structureType, structureDF, structureWeight)) 
            for _, row in rgxdf.iterrows()
        )

    elif cpumode == "multi":
        with Pool(48) as p:
            structure_weight_list = p.map(compute_structure_weight, 
                                            [(row, structureType, structureDF, structureWeight) 
                                            for _, row in rgxdf.iterrows()])
    else:
        raise ValueError("Invalid cpumode. Use 'single', 'joblib', or 'multi'.")
    # 添加计算结果
    rgxdf["Rgx_rawvalue"] = rgxdf["Rgx_rawvalue"] * structure_weight_list
    rgxdf["Rgx_percent"] = rgxdf["Rgx_percent"] * structure_weight_list
    #rgxdf["Rgx_percent"] = rgxdf["Rgx_rawvalue"] / rgxdf.groupby('geneID')['Rgx_rawvalue'].transform('sum')
    #rgxdf["Rgx_percent"].fillna(0)
    #rgxdf.to_csv(outname+".tsv", sep="\t", index=False)

    # 处理 rgfile
    if rgfile:
        if ignorehead:
            rgdf = pd.read_csv(rgfile, sep="\t", header=None, skiprows=1)
        else:
            rgdf = pd.read_csv(rgfile, sep="\t", header=None)
        rgdf.index = rgdf[rggeneID]
        rg_adj_structure = rgxdf.groupby("geneID")["Rgx_rawvalue"].sum()
        rgdf[rgvalue] = rg_adj_structure

    return rgdf, rgxdf



def diffrentThresh(rgxfile_raw,rgfile_raw,golddf,goldcol,outdir,tpmfile,ignoreRgxhead=True,
                   threshType="tpm",threshList=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],
                   tmpcolrep=[2+1,3+1,],rggeneID=4+1,tmpgeneID=1+1,rgxgeneID=9+1,plotcols=None,rgvaluecol=6+1,
                   matchcol=12+1,truecol=11+1,onlyPlot=False,goldwithhead=False):
    
    if ignoreRgxhead:
        rgxdf=pd.read_csv(rgxfile_raw,sep="\t",header=None,skiprows=1)
    else:
        rgxdf=pd.read_csv(rgxfile_raw,sep="\t",header=None)
        
    rgdf = pd.read_csv(rgfile_raw,sep="\t",header=None)
    tpmdf = pd.read_csv(tpmfile,sep="\t")
    
    rgdfgene = rgdf[rggeneID-1].apply(lambda x: x.split('.')[0])
    rgdf.index = rgdfgene
    tpmdfgene = tpmdf.iloc[:,tmpgeneID-1]
    tpmdf.index = tpmdfgene
    commongene = np.array(set(tpmdfgene) & set(rgdfgene))
    rgdffinal = rgdf.loc[commongene]    
    rgdffinal["TPM"] = tpmdf.loc[commongene].iloc[:,np.array(tmpcolrep)-1].mean(axis=1)
    
    namelist=[]
    for th in threshList:
        if threshType=="tpm":
            definedcol = [mcolors.to_hex(plt.cm.get_cmap('inferno')(i)) for i in np.linspace(0, 1, 15)]

            name = "TPM"+"_q"+str(th)
            namelist.append(name)
            if not onlyPlot:
                rgdfbool = rgdffinal["TPM"] > rgdffinal["TPM"].quantile(th)
                selectedgene = list(rgdfbool[rgdfbool].index)
                rgxgene = rgxdf.iloc[:,rgxgeneID-1].apply(lambda x: x.split('.')[0])
                rgxdf_select = rgxdf[rgxgene.isin(selectedgene)]
        elif threshType=="sumrank":
            definedcol = [mcolors.to_hex(plt.cm.get_cmap('inferno_r')(i)) for i in np.linspace(0, 1, 15)]

            name = "sumrank_Rg-TPM"+"_q"+str(th)
            namelist.append(name)
            if not onlyPlot:
                sumrank = makesumrank(rgdffinal["TPM"],rgdffinal[rgvaluecol-1])
                rgdfbool = sumrank > sumrank.quantile(th)
                selectedgene = list(rgdfbool[rgdfbool].index)
                rgxgene = rgxdf.iloc[:,rgxgeneID-1].apply(lambda x: x.split('.')[0])
                rgxdf_select = rgxdf[rgxgene.isin(selectedgene)]
        elif threshType=="distance":
            definedcol = [mcolors.to_hex(plt.cm.get_cmap('Spectral_r')(i)) for i in np.linspace(0, 1, 20)]

            name = threshType+"_"+str(th)+"bp"
            namelist.append(name)
            if not onlyPlot:
                rgxdf_select,_ = cut_distance(rgxdf,rgdf,maxdis=th)
                
        if not onlyPlot:
            rgxdf_select.to_csv(outdir+name+".selectRgx.tsv",sep="\t",header=None,index=False)
            matchgold(outdir+name+".selectRgx.tsv",golddf,
                      outdir+name+".tsv",goldcol,percent=False,withhead=goldwithhead)
    if plotcols:
        finalcolor=plotcols
    else:
        finalcolor=definedcol

    plotMultiPrc(outdir, matchcol=matchcol-1,truecol=truecol-1,dataset=outdir,cols=finalcolor,
            typelist=namelist,matchedwithhead=goldwithhead,outname=f"diffrentThresh_{threshType}.pdf")
    


def plotdensity(rgxfile_raw,indir,matchedwithhead = False,matchcol=12,truecol=11):
    
    rgxdf=pd.read_csv(rgxfile_raw,sep="\t",header=None,skiprows=1)
    
    normrgx = np.log2(rgxdf[11] / rgxdf[11].mean() +1)
    normrgx_non0 = normrgx[normrgx>0]
    
    plt.figure(figsize=(3,3))
    kdeplot = sns.kdeplot(normrgx_non0, color="#6A5ACD", fill=False, linewidth=2, label="Density")
    # 阈值线位置
    threshold_x = np.log2(1 + 1)  # 即 log2(2) = 1
    # 添加虚线
    plt.axvline(x=threshold_x, color='#2E8B57', linestyle='--', linewidth=1.5, label="normRgx = 1")
    # 添加阴影区域（大于阈值部分）
    x_vals = np.linspace(normrgx_non0.min(), normrgx_non0.max(), 1000)
    kde_vals = kdeplot.get_lines()[0].get_data()
    x_kde, y_kde = kde_vals
    # 用 fill_between 填充大于阈值的区域
    plt.fill_between(x_kde, y_kde, where=(x_kde > threshold_x), color='orange', alpha=0.3, label="Above threshold")
    # 图形标题和标签
    plt.title(indir, fontsize=11)
    plt.xlabel("normRgx (exclude 0)", fontsize=10)
    plt.ylabel("Density", fontsize=10)
    plt.legend(fontsize=9)
    # plt.tight_layout()
    plt.show()
    
    typelist_raw=["adjRaw","adjRankRaw","adjRawStruc","adjStrucRankRaw"]
    plt.figure(figsize=(3,3))
    for name in typelist_raw:
        if not matchedwithhead:
            matched = pd.read_csv(indir+name+".tsv",sep="\t",header=None)
        else:
            matched = pd.read_csv(indir+name+".tsv",sep="\t")
    
        rgxvalue = matched[matchcol-1]
        truelabel = matched[truecol-1] == True
        normRgx = np.log2(rgxvalue/rgxvalue.mean()+1)
        
        precision, recall, thresholds = precision_recall_curve(truelabel,normRgx)
        
        recall_cutoffs = np.arange(0.1, 1.0, 0.1) 
        print(name + " 不同 recall 水平对应的阈值：")
        
        recall_threshold_list = []
        for target_recall in recall_cutoffs:
            idx = np.argmin(np.abs(recall[1:] - target_recall))  # recall[1:] 对应 thresholds
            corresponding_threshold = round(thresholds[idx],4)
            recall_threshold_list.append(corresponding_threshold)
        print(recall_cutoffs)
        print(recall_threshold_list)

        sns.kdeplot(normRgx,label=name)
    plt.xlabel("Rgx (raw)", fontsize=10)
    plt.title(indir, fontsize=11)
    plt.legend()
    plt.show()
    
    typelist_percent=["noadj","adjRank","adjStruc","adjStrucRank",]
    plt.figure(figsize=(3,3))
    for name in typelist_percent:
        if not matchedwithhead:
            matched = pd.read_csv(indir+name+".tsv",sep="\t",header=None)
        else:
            matched = pd.read_csv(indir+name+".tsv",sep="\t")
        
        rgxvalue = matched[matchcol-1]
        sns.kdeplot(rgxvalue[rgxvalue>0],label=name)
    plt.xlabel("Rgx (percent)", fontsize=10)
    plt.title(indir, fontsize=11)
    plt.legend()
    plt.show()


def adjvalue(rgxfile_raw,rgfile_raw,outdir,tpmfile,
                   structureTypeList,structureFileList,structureWeightList,
                   tmpcolrep=[2+1,3+1,],ignorehead=True,tmpgeneID=1+1,ranktype="diffrank",):
    
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Step 1: Adjust raw signal by structure
    rg_adjrawstruc, rgx_adjrawstruc = adjStructure(
        rgxfile_raw,
        structureFileList,
        structureTypeList,
        structureWeightList,
        rgfile=rgfile_raw,
        rggeneID=4,  # fixed 1-based column for gene ID
        rgvalue=6,   # fixed 1-based column for raw Rg signal
        ignorehead=ignorehead,
        cpumode="multi",
        withcolname=True
    )

    # Save adjusted raw structure signals
    rg_adjrawstruc.to_csv(outdir + "adjrawstruc.rg.tsv", sep="\t", header=None, index=False)
    rgx_adjrawstruc.to_csv(outdir + "adjrawstruc.rgx.tsv", sep="\t", header=None, index=False)

    # Step 2: Adjust by TPM signal
    rg_adjall, rgx_adjall = adjtpm(
        outdir + "adjrawstruc.rgx.tsv",
        tpmfile,
        rgfile=outdir + "adjrawstruc.rg.tsv",
        rggeneID=4,
        rgxgeneID=9,
        ranktype=ranktype,
        rgxvalue=11,
        rgxratio=12,
        tmpcolrep=np.array(tmpcolrep) - 1,
        ignorehead=ignorehead,
        tmpgeneID=tmpgeneID - 1
    )

    # Save final adjusted signals
    rg_adjall.to_csv(outdir + "adjall.rg.tsv", sep="\t", header=None, index=False)
    rgx_adjall.to_csv(outdir + "adjall.rgx.tsv", sep="\t", header=None, index=False)
    
    # Clean up intermediate files
    try:
        os.remove(outdir + "adjrawstruc.rg.tsv")
        os.remove(outdir + "adjrawstruc.rgx.tsv")
    except FileNotFoundError:
        print("Some intermediate files were not found for deletion.")