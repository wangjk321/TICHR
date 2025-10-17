import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
import warnings
warnings.filterwarnings('ignore')
from scipy.stats import wilcoxon,ttest_1samp
import os

from .PRC_ROC import *



def benchmarkvalue(RgDF_Ctrl,RgDF_Treat,degtype="all"):
    #meanRg = (np.log1p(RgDF_Treat[9]) + np.log1p(RgDF_Ctrl[9]))/2
    #meanTPM = (RgDF_Ctrl[8] + RgDF_Treat[8])/2
    deltaRg = np.log1p(RgDF_Treat[9]) -  np.log1p(RgDF_Ctrl[9])
    #deltaTPM = RgDF_Ctrl[6]
    if degtype=="all":
        degbool = (abs(RgDF_Ctrl[6])>1) & (RgDF_Ctrl[7]<0.05)
    elif degtype=="up":
        degbool = (RgDF_Ctrl[6])>1 & (RgDF_Ctrl[7]<0.05)
    elif degtype=="down":
        degbool = (RgDF_Ctrl[6])<-1 & (RgDF_Ctrl[7]<0.05)
    else:
        exit(1)
    
    return(abs(deltaRg),degbool)

def cut_distance(df,rgdf,maxdis=100000): #df是RgxDf
    peakpos = (df[1] + df[2]) / 2
    tsspos = np.where(df.iloc[:, 8] == "+", df.iloc[:, 6], df.iloc[:, 7])
    peak2tss = abs(peakpos-tsspos)
    cutdf = df[peak2tss< maxdis]
    cutdf.reset_index(drop=True, inplace=True)
    
    rgdf.index = rgdf[3]
    rgdf[9] = cutdf.groupby(4)[11].sum()
    cutrgdf= rgdf.dropna(subset=[rgdf.columns[9]])
    cutrgdf.reset_index(drop=True, inplace=True)
    
    return(cutdf,cutrgdf)




#选取变化最大的基因
def getQprc(RgDF_Ctrl,RgDF_Treat,RgxDF_Ctrl,RgxDF_Treat,selecttype="rg",thresh=0.8,returnDF=False):
    changeRgxDF = pd.DataFrame()
    changeRgxDF["geneID"] = RgxDF_Ctrl[4]
    changeRgxDF["change_epix"] = RgxDF_Treat[3] - RgxDF_Ctrl[3]
    changeRgxDF["change_3dx"] = RgxDF_Treat[10] - RgxDF_Ctrl[10]
    changeRgxDF["change_rgx"] =  RgxDF_Treat[11] - RgxDF_Ctrl[11]
    changeRgxDF["avergae_rgxpercent"] =  (RgxDF_Treat[12] + RgxDF_Ctrl[12])/2
    changeRgxDF["change_rgxpercent"] =  RgxDF_Treat[12] - RgxDF_Ctrl[12]
    changeRgxDF["change_fcepix"] = np.log1p(RgxDF_Treat[3]) - np.log1p(RgxDF_Ctrl[3])
    changeRgxDF["change_fc3dx"] = np.log1p(RgxDF_Treat[10]) - np.log1p(RgxDF_Ctrl[10])
    changeRgxDF["change_fcrgx"] =  np.log1p(RgxDF_Treat[11]) - np.log1p(RgxDF_Ctrl[11])
    
    genesum_changeEpi= changeRgxDF.groupby("geneID")["change_epix"].sum()
    genesum_change3d = changeRgxDF.groupby("geneID")["change_3dx"].sum()
    genesum_changerg = changeRgxDF.groupby("geneID")["change_rgx"].sum()
    change_rg = abs(genesum_changerg) >= abs(genesum_changerg).quantile(thresh)
    change_epi = abs(genesum_changeEpi) >= abs(genesum_changeEpi).quantile(thresh)
    change_3d = abs(genesum_change3d) >= abs(genesum_change3d).quantile(thresh)
    
    genesum_fcEpi= np.log1p(RgxDF_Treat.groupby(4)[3].sum())-np.log1p(RgxDF_Ctrl.groupby(4)[3].sum())
    genesum_fc3d= np.log1p(RgxDF_Treat.groupby(4)[10].sum())-np.log1p(RgxDF_Ctrl.groupby(4)[10].sum())
    genesum_fcRg= np.log1p(RgxDF_Treat.groupby(4)[11].sum())-np.log1p(RgxDF_Ctrl.groupby(4)[11].sum())
    change_fcEpi = abs(genesum_fcEpi) >= abs(genesum_fcEpi).quantile(thresh)
    change_fc3d = abs(genesum_fc3d) >= abs(genesum_fc3d).quantile(thresh)
    change_fcrg = abs(genesum_fcRg) >= abs(genesum_fcRg).quantile(thresh)
    
    if selecttype == 'rg':
        geneSelectBool = change_rg 
    elif selecttype == 'epi':
        geneSelectBool = change_epi
    elif selecttype == '3d':
        geneSelectBool = change_3d
    elif selecttype == 'fcrg':
        geneSelectBool = change_fcrg
    elif selecttype == 'fcepi':
        geneSelectBool = change_fcEpi
    elif selecttype == 'fc3d':
        geneSelectBool = change_fc3d
        
    geneSelect = np.array(geneSelectBool[geneSelectBool].index)
    geneSelectFinal = RgDF_Ctrl[3].isin(geneSelect)
    
    if not returnDF:
        return(geneSelectFinal)
    else:
        RgDFchange = RgDF_Ctrl.copy()
        RgDFchange.index = RgDFchange[3]
        RgDFchange['change_epi'] = genesum_changeEpi
        RgDFchange['change_3d'] = genesum_change3d
        RgDFchange['change_rg'] = genesum_changerg
        RgDFchange['change_fcepi'] = genesum_fcEpi
        RgDFchange['change_fc3d'] = genesum_fc3d
        RgDFchange['change_fcrg'] = genesum_fcRg
        
        return(geneSelectFinal,RgDFchange,changeRgxDF)
    

def randomGene(RgDF_Ctrl,RgDF_Treat,selectbool,n=100,seed=42):
    np.random.seed(seed)
    prclist=[] #长度为100
    precisionlist=[]
    recalllist=[]

    for i in range(n):
        shuffledBool = np.random.permutation(selectbool.values)

        abschange_random, degbool_random = benchmarkvalue(RgDF_Ctrl[shuffledBool],
                                                          RgDF_Treat[shuffledBool])
        prclist.append(average_precision_score(degbool_random, abschange_random))

        precision, recall, _ = precision_recall_curve(
            degbool_random, #防止abschange有相同的值
            abschange_random+np.random.normal(0, 1e-10, abschange_random.shape)
        )

        precisionlist.append(precision)
        recalllist.append(recall)

    recall_median = np.median(recalllist,axis=0) #长度为1303
    precision_median = np.median(precisionlist,axis=0)

    recall_common = np.linspace(0, 1, len(selectbool))

    recall_q5 = np.percentile(recalllist, 5, axis=0)
    recall_q95 = np.percentile(recalllist, 95, axis=0)
    precision_q5 = np.percentile(precisionlist, 5, axis=0)
    precision_q95 = np.percentile(precisionlist, 95, axis=0)

    # 按 recall_q5 排序
    sorted_indices_q5 = np.argsort(recall_q5)  # 获取 recall_q5 的排序索引
    recall_q5_sorted = recall_q5[sorted_indices_q5]  # 排序后的 recall_q5
    precision_q5_sorted = precision_q5[sorted_indices_q5]  # 对应的 precision_q5 排序

    # 对 recall_q95 和 precision_q95 按 recall_q5 排序（如果需要）
    sorted_indices_q95 = np.argsort(recall_q95)
    recall_q95_sorted = recall_q95[sorted_indices_q95]
    precision_q95_sorted = precision_q95[sorted_indices_q95]

    precision_q5_interp = np.interp(recall_common, recall_q5_sorted, precision_q5_sorted)
    precision_q95_interp = np.interp(recall_common, recall_q95_sorted, precision_q95_sorted)

    return(prclist,recall_median,precision_median,recall_common,precision_q5_interp,precision_q95_interp) 

class DiffEvent:
    def __init__(self,RgDF_Ctrl_file,RgxDF_Ctrl_file,RgDF_Treat_file,RgxDF_Treat_file,maxdistance=500000,
                 outdir=os.getcwd(),inputtype="file",seed=42, pdf=True):
        self.seed=seed
        print(inputtype)
        if inputtype=="file":
            RgDF_Ctrl_raw = pd.read_csv(RgDF_Ctrl_file,header=None,sep="\t")
            RgxDF_Ctrl_raw = pd.read_csv(RgxDF_Ctrl_file,header=None,sep="\t")
            RgDF_Treat_raw = pd.read_csv(RgDF_Treat_file,header=None,sep="\t")
            RgxDF_Treat_raw = pd.read_csv(RgxDF_Treat_file,header=None,sep="\t")
        elif inputtype=="pandas":
            RgDF_Ctrl_raw = RgDF_Ctrl_file
            RgxDF_Ctrl_raw = RgxDF_Ctrl_file
            RgDF_Treat_raw = RgDF_Treat_file
            RgxDF_Treat_raw = RgxDF_Treat_file

        self.maxdistance = maxdistance
        self.RgxDF_Ctrl, self.RgDF_Ctrl = cut_distance(RgxDF_Ctrl_raw,RgDF_Ctrl_raw,maxdis=maxdistance)
        self.RgxDF_Treat, self.RgDF_Treat = cut_distance(RgxDF_Treat_raw,RgDF_Treat_raw,maxdis=maxdistance)
        self.outdir = outdir
        self.pdf=pdf
    
    def predictDEG(self,type="PRC",label="label",custom=False,linecolor=None,degtype='all'):
        abschange, degbool = benchmarkvalue(self.RgDF_Ctrl,self.RgDF_Treat,degtype=degtype)
        if custom:
            if type=="PRC":
                pltoneprc(degbool,abschange,label,linecolor)
            elif type=="ROC":
                pltoneroc(degbool,abschange,label,linecolor)
        else:
            if type=="PRC":
                showAUPRC(degbool,abschange,label)
            elif type=="ROC":
                showROC(degbool,abschange,label)

    def quantilePRC(self,selectGeneType="rg",label="label",plotbg=True,plotpvalue='wilcox',
                    thlist=[0.99, 0.98, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0],bgalpha=0.5,
                    plotylim=[0,0.7],plotallp=False):
        plt.figure(figsize=(5.5,3.8))
        collist = [mcolors.to_hex(plt.cm.get_cmap('Spectral_r')(i)) for i in np.linspace(0, 1, len(thlist))]
        i=0
        auprclist=[]
        geneSelectFinalList = []
        for threshhold in thlist:

            geneSelectFinal = getQprc(self.RgDF_Ctrl,self.RgDF_Treat,self.RgxDF_Ctrl,self.RgxDF_Treat,
                           selecttype=selectGeneType,thresh=threshhold,)
            abschange, degbool = benchmarkvalue(self.RgDF_Ctrl[geneSelectFinal],self.RgDF_Treat[geneSelectFinal])
            pltoneprc(degbool,abschange,selectGeneType+"_q"+str(threshhold)+":",collist[i])
            auprclist.append(average_precision_score(degbool,abschange))
            geneSelectFinalList.append(geneSelectFinal)
            i=i+1
        plt.legend(loc="center left", bbox_to_anchor=(1, 0.5))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(label+"_max"+str(self.maxdistance)+"bp")
        plt.tight_layout()
        if self.outdir:
            plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_multi_prcSelect.pdf")
        
        plt.figure(figsize=(5,3.8))
        if plotbg:
            prc_random_lower = []
            prc_random_upper = []
            ttestp_list=[]
            wilcoxp_list = []
            j=0
            for threshhold in thlist:
                print(threshhold)
                prclist,recall_median,precision_median,recall_common,precision_q5_interp,precision_q95_interp = \
                randomGene(self.RgDF_Ctrl,self.RgDF_Treat,geneSelectFinalList[j],seed=self.seed)
                prc_random_lower.append(np.percentile(prclist, 5))
                prc_random_upper.append(np.percentile(prclist, 95))
                #plt.step(recall_median,precision_median,where='post',color=collist[j],label=label)
                plt.fill_between(recall_common, precision_q5_interp, precision_q95_interp, 
                    color=collist[j], alpha=bgalpha, label=selectGeneType+"_q"+str(threshhold))
                
                _,wilcoxp = wilcoxon([x - auprclist[j] for x in prclist],alternative='less')
                _,ttestp = ttest_1samp(prclist, auprclist[j], alternative='less')
                ttestp_list.append(ttestp)
                wilcoxp_list.append(wilcoxp)

                j=j+1

        plt.legend(bbox_to_anchor=(1, 0.5),loc='center left')
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(label+"_max"+str(self.maxdistance)+"bp")
        plt.tight_layout()
        if self.outdir:
            plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_multi_prcBg.pdf")
            # plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_multi_prcBg.png",dpi=300)
        
        if plotallp:
            plt.figure(figsize=(5.4,3.8))
        else:
            plt.figure(figsize=(3.5,3.8))
        plt.scatter(thlist, auprclist, c=collist, s=70, alpha=1, edgecolor='black',)
        line_main, = plt.plot(thlist, auprclist, color='black', linestyle='-', linewidth=3, alpha=0.5, 
                 label="selected genes")
        bg_patch = plt.fill_between(thlist, prc_random_lower, prc_random_upper, 
                    color='lightgrey', alpha=0.5, label="background")
        
        stat_handles = []

        if plotpvalue == "wilcox":
            medianp=np.median(np.nan_to_num(wilcoxp_list, nan=1.0))
            line_medianp,=plt.plot([],[]," ",label="median p="+str("{:.1e}".format(medianp)))
            for i,t in enumerate(thlist):
                pvalue = str("{:.1e}".format(wilcoxp_list[i]))
                dummy, = plt.plot([],[]," ", label="q"+str(t)+" p="+pvalue)
                stat_handles.append(dummy)
        elif plotpvalue == "ttest":
            medianp=np.median(np.nan_to_num(ttestp_list, nan=1.0))
            line_medianp,=plt.plot([],[]," ",label="median p="+str("{:.1e}".format(medianp)))
            for i,t in enumerate(thlist):
                pvalue = str("{:.1e}".format(ttestp_list[i]))
                dummy, = plt.plot([],[]," ", label="q"+str(t)+" p="+pvalue)
                stat_handles.append(dummy)
        else:
            pass

        main_legend = plt.legend(handles=[line_main, bg_patch,line_medianp], loc='upper left')
        plt.gca().add_artist(main_legend)
        if plotallp:
            plt.legend(handles=stat_handles,bbox_to_anchor=(1, 0.5),loc='center left',fontsize=9.5)

        plt.xlabel("select_"+selectGeneType+' quantile')
        plt.ylabel('AUPRC')
        plt.title(label+"\n(select vs. background)")
        if plotylim: plt.ylim(plotylim)
        plt.tight_layout()
        if self.outdir:
            plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_multi_auprcSelectVSbg.pdf")
    
    def selectgene(self,selectGeneType="rg",label="label",threshhold=0.9,plot=True,plotbg=True):      
        plt.figure(figsize=(3.2,3.5))
        geneSelectFinal = getQprc(self.RgDF_Ctrl,self.RgDF_Treat,self.RgxDF_Ctrl,self.RgxDF_Treat,
                           selecttype=selectGeneType,thresh=threshhold)
        if not plot:
            return(geneSelectFinal)
        
        abschange, degbool = benchmarkvalue(self.RgDF_Ctrl[geneSelectFinal],self.RgDF_Treat[geneSelectFinal])
        selectprc = average_precision_score(degbool,abschange)
        
        pltoneprc(degbool,abschange,selectGeneType+"_q"+str(threshhold),"m")

        #return(geneSelectFinal)
        if plotbg:
            prclist,recall_median,precision_median,recall_common,precision_q5_interp,precision_q95_interp = \
            randomGene(self.RgDF_Ctrl,self.RgDF_Treat,geneSelectFinal,seed=self.seed)
            prclist_q5 = np.percentile(prclist, 5)
            prclist_q95 = np.percentile(prclist, 95)
            plt.step(recall_median,precision_median,where='post',color="grey",label="background median")
            plt.fill_between(recall_common, precision_q5_interp, precision_q95_interp, 
                color="grey", alpha=0.1, label="background 95%")
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(label+"\n(select vs. background)")    
        plt.legend(loc='best')
        plt.tight_layout()
        if self.outdir:
            plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_q"+str(threshhold)+"_prcSelectVSbg.pdf")
            # plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_q"+str(threshhold)+"_prcSelectVSbg.png",dpi=300)
        
        _,wilcoxp = wilcoxon([x - selectprc for x in prclist],alternative='less')
        _,ttestp = ttest_1samp(prclist, selectprc, alternative='less')
        print(wilcoxp,ttestp)
        
        plt.figure(figsize=(3.5,3.5))
        sns.kdeplot(prclist,shade=True, color='royalblue', alpha=0.3,label="Background")
        plt.axvline(selectprc, color='indigo', linestyle='--',label=selectGeneType+"_q"+str(threshhold))
        
        ax = plt.gca()
        plt.plot([],[],' ',label="ttest p="+str("{:.2e}".format(ttestp)))
        plt.plot([],[],' ',label="wilcox p="+str("{:.2e}".format(wilcoxp)))
        plt.legend(loc="best")
        plt.title(label+"\n(select vs. background)") 
        plt.xlabel('AUPRC')
        if self.outdir:
            plt.savefig(self.outdir+"/"+label+"_select"+selectGeneType+"_q"+str(threshhold)+"_auprcSelectVSbg.pdf")
        
    
    def diffgene(self,selectGeneType="rg",quantiTh=0.9,selectQuantile=True,
                 selectDeg=True,selectFc=True,fcTh=0.5,
                 selectRank="same",sameTh=0.9,diffTh=0.9,plotxlim=None,plotxlabel=None,
                 title="title"):
        
        # 选定特征的前百分之几
        selectQuantiBool,RgDFchange,changeRgxDF=getQprc(self.RgDF_Ctrl,self.RgDF_Treat,
                                                        self.RgxDF_Ctrl,self.RgxDF_Treat,
                                selecttype=selectGeneType,thresh=quantiTh,returnDF=True)
        selectQuantiBool.index = RgDFchange[3]
        
        totalbool = RgDFchange[1] >= 0 #定义一个全部为True的pd.Series
        if selectQuantile: totalbool = totalbool & selectQuantiBool
        
        self.RgDFchange = RgDFchange
        self.changeRgxDF = changeRgxDF
        
        geneFC = RgDFchange[6]
        geneFDR = -np.log10(RgDFchange[7])
        selectTypeChange = RgDFchange["change_"+selectGeneType]
        pointsize = np.interp(geneFDR, (geneFDR.min(), geneFDR.max()), (1, 20))
        
        plt.figure(figsize=(4.5,4))
        if selectRank=="diff":
            scatter_with_rank(selectTypeChange,geneFC,"diffrank",s=pointsize,alpha=0.3,ifcolbar=False)
        else:
            scatter_with_rank(selectTypeChange,geneFC,"sumrank0",s=pointsize,alpha=0.3,ifcolbar=False)        

        plt.colorbar(shrink=0.2,aspect=10)

        if selectDeg:   #差异表达基因
            selectDegBool = (abs(RgDFchange[6])>1) & (RgDFchange[7]<0.05)
            totalbool = totalbool & selectDegBool
            
        if selectFc:  #变化倍数
            selectFcBool = abs(RgDFchange["change_fc"+selectGeneType.removeprefix("fc")])>fcTh
            totalbool = totalbool & selectFcBool
            
        if selectRank == "same": #变化趋势
            rankBool = abs(makesumrank_center0(RgDFchange["change_"+selectGeneType], RgDFchange[6]))>sameTh
        elif selectRank == "diff":
            rankBool = abs(makediffrank(RgDFchange["change_"+selectGeneType], RgDFchange[6]))>diffTh
        totalbool = totalbool & rankBool
        
        legend_sizes = [0.05,1e-5,1e-10]  # 代表性点大小
        legend_handles = [plt.scatter([], [], 
                                      s=np.interp(-np.log10(size), (geneFDR.min(), geneFDR.max()), (1, 20)), 
                                      facecolors='none', edgecolors='grey',alpha=0.5, label=f"geneFDR={size}") 
                                      for size in legend_sizes]
        #plt.legend(handles=legend_handles, loc="upper right", title="Point Size",bbox_to_anchor=(1, 0.5))
        
        plt.scatter(RgDFchange[totalbool]["change_"+selectGeneType],
                    RgDFchange[totalbool][6],alpha=0.8, edgecolors='black',facecolors='m',
                    s=pointsize[totalbool],color="m",label="Selected genes")
        if plotxlabel: plt.xlabel(plotxlabel)
        if plotxlim: plt.xlim(plotxlim)
        plt.ylabel("Gene logFC")
        plt.title(title)
        plt.legend(loc="best")
        plt.tight_layout()
        if self.outdir:
            if self.pdf:
                plt.savefig(self.outdir+"/"+title+"_select"+selectGeneType+"_"+selectRank+"_diffgene_simple.pdf")
            else:
                plt.savefig(self.outdir+"/"+title+"_select"+selectGeneType+"_"+selectRank+"_diffgene_simple.png",dpi=300)

        print(sum(totalbool)," genes selected")

        ##################画选择图
        plt.figure(figsize=(5.5,3.8))
        
        if selectDeg:
            plt.axhline(-1, color='r', linestyle='--',alpha=0.5,)  
            plt.axhline(1, color='r', linestyle='--',label="abs(geneFC)>1",alpha=0.5,) 
        
        plt.axvline(-abs(selectTypeChange).quantile(quantiTh), color='b', linestyle='--',alpha=0.5)  # 画一条 y=-1 的红色虚线
        plt.axvline(abs(selectTypeChange).quantile(quantiTh),
                    color='b', linestyle='--',label="Quantile>"+str(quantiTh),alpha=0.5) # 画一条 y=-1 的红色虚线
        if selectFc:
            plt.scatter(RgDFchange[selectFcBool]["change_"+selectGeneType],
                RgDFchange[selectFcBool][6],s=pointsize[selectFcBool],
                        alpha=0.5,label="abs(FC)>0.5",edgecolors='none')
        if selectRank in ["same","diff"]:
            plt.scatter(RgDFchange[rankBool]["change_"+selectGeneType],
                        RgDFchange[rankBool][6],s=pointsize[rankBool],alpha=0.5,
                        color="g",label="Trends",edgecolors='none')
            
        plt.scatter(selectTypeChange,geneFC,s=pointsize,
                    facecolors='none', edgecolors='grey',alpha=0.5)
        legend_sizes = [0.05,1e-5,1e-10]  # 代表性点大小
        legend_handles = [plt.scatter([], [], 
                                      s=np.interp(-np.log10(size), (geneFDR.min(), geneFDR.max()), (1, 20)), 
                                      facecolors='none', edgecolors='grey',alpha=0.5, label=f"geneFDR={size}") 
                                      for size in legend_sizes]
        #plt.legend(handles=legend_handles, loc="upper right", title="Point Size",bbox_to_anchor=(1, 0.5))
        
        plt.scatter(RgDFchange[totalbool]["change_"+selectGeneType],
                    RgDFchange[totalbool][6],alpha=0.8, edgecolors='black',facecolors='m',
                    s=pointsize[totalbool],color="m",label="Selected genes")
        if plotxlabel: plt.xlabel(plotxlabel)
        if plotxlim: plt.xlim(plotxlim)
        plt.ylabel("Gene logFC")
        plt.title(title)
        plt.legend(bbox_to_anchor=(1, 0.5),loc='center left',)
        plt.tight_layout()
        if self.outdir:
            if self.pdf:
                plt.savefig(self.outdir+"/"+title+"_select"+selectGeneType+"_"+selectRank+"_diffgene_full.pdf")
            else:
                plt.savefig(self.outdir+"/"+title+"_select"+selectGeneType+"_"+selectRank+"_diffgene_full.png",dpi=300)
        #plt.legend(bbox_to_anchor=(1, 0.5))
        self.totalbool = totalbool
        self.finalgene = RgDFchange[totalbool]
    
    def diffEP(self,selectdiffgene=True,selectRgxType="rgx",
               selectRgxQuantile=True,quantileRgxCut=0.9,
               selectRgxFold =True,RgxFcCut=0.5,
               selectRgxRank = "same",rgxRankCut=0.8,
               selectRgxRatio=True,RgxRatioCut=0.01,
               plotxlim=None,plotxlabel=None,title="title"):
        geneFCdict = self.RgDF_Ctrl.set_index(3)[6].to_dict()
        geneFDRdict = self.RgDF_Ctrl.set_index(3)[7].to_dict()
        geneTPMdict = self.RgDF_Ctrl.set_index(3)[8].to_dict()
        RgxDFchange = pd.concat([self.RgxDF_Ctrl,self.changeRgxDF],axis=1)
        RgxDFchange["geneFC"] = RgxDFchange[4].map(geneFCdict)
        RgxDFchange["geneFDR"]= RgxDFchange[4].map(geneFDRdict)
        RgxDFchange["geneTPM"] = RgxDFchange[4].map(geneTPMdict)
        rgxTotalBool = RgxDFchange[1] >= -100000000000
        
        selectRgxChange = RgxDFchange['change_'+selectRgxType]
        
        if selectRgxType=="rgxpercent": selectRgxFold=False
        
        #第一个根据基因选择
        if selectdiffgene:
            Rgx_selectgene_bool = RgxDFchange["geneID"].isin(self.totalbool[self.totalbool].index)
            print("Selected ",sum(Rgx_selectgene_bool)," site-to-gene links by selected genes")
            rgxTotalBool = rgxTotalBool & Rgx_selectgene_bool
        
        if selectRgxRatio:
            majorRgxBool = RgxDFchange["avergae_rgxpercent"] >= RgxRatioCut
            print("Selected ",sum(majorRgxBool)," site-to-gene links by major RgX")
            rgxTotalBool = rgxTotalBool & majorRgxBool

        if selectRgxQuantile:
            #第二个根据Quantile选择
            changeRgxBool = abs(selectRgxChange) >= abs(selectRgxChange).quantile(quantileRgxCut)
            print("Selected ",sum(changeRgxBool)," site-to-gene links by Rgx quantile")
            rgxTotalBool = rgxTotalBool & changeRgxBool

        if selectRgxFold:
            #第三个根据变化倍数选择
            print("selectRgxFold")
            selectRgxFc = RgxDFchange["change_fc"+selectRgxType.removeprefix("fc")]
            selectRgxFcBool = abs(selectRgxFc)>= RgxFcCut
            print("Selected ",sum(selectRgxFcBool)," site-to-gene links by Rgx log fold change")
            rgxTotalBool = rgxTotalBool & selectRgxFcBool
        
        if selectRgxRank == "same": #变化趋势
            rankBool = abs(makesumrank_center0(selectRgxChange, RgxDFchange["geneFC"]))> rgxRankCut
            print("Selected ",sum(rankBool)," site-to-gene links by sumrank")
            rgxTotalBool = rgxTotalBool & rankBool
        elif selectRgxRank == "diff":
            rankBool = abs(makediffrank(selectRgxChange, RgxDFchange["geneFC"]))> rgxRankCut
            print("Selected ",sum(rankBool)," site-to-gene links by diffrank")
            rgxTotalBool = rgxTotalBool & rankBool
        else:
            pass
        
        plt.figure(figsize=(4.5,4))
        geneFDR=-np.log10(RgxDFchange["geneFDR"])
        pointsize = np.interp(geneFDR, (geneFDR.min(), geneFDR.max()), (1, 20))
        if selectRgxRank == "diff":
            scatter_with_rank(selectRgxChange,RgxDFchange["geneFC"],"diffrank",s=pointsize,alpha=0.3,ifcolbar=False)
        else:
            scatter_with_rank(selectRgxChange,RgxDFchange["geneFC"],"sumrank0",s=pointsize,alpha=0.3,ifcolbar=False)

        plt.colorbar(shrink=0.2,aspect=10)
        
        legend_sizes = [0.05,1e-5,1e-10]  # 代表性点大小
        legend_handles = [plt.scatter([], [], 
                                      s=np.interp(-np.log10(size), (geneFDR.min(), geneFDR.max()), (1, 20)), 
                                      facecolors='none', edgecolors='grey',alpha=0.5, label=f"geneFDR={size}") 
                                      for size in legend_sizes]

        plt.scatter(selectRgxChange[rgxTotalBool],
                    RgxDFchange["geneFC"][rgxTotalBool],alpha=0.8, edgecolors='black',facecolors='m',
                    s=pointsize[rgxTotalBool],color="m",label="Selected genes")
        if plotxlim:plt.xlim(plotxlim)
        if plotxlabel: plt.xlabel(plotxlabel)
        plt.title(title)

        print("Selected ",sum(rgxTotalBool)," site-to-gene links overall")
        plt.legend(loc="best")
        plt.ylabel("Gene logFC")
        plt.tight_layout()

        if self.outdir:
            plt.savefig(self.outdir+"/"+title+"_select"+selectRgxType+"_diffs2g_simple.pdf")
            # plt.savefig(self.outdir+"/"+title+"_select"+selectRgxType+"_diffs2g_simple.png",dpi=300)

        #######
        plt.figure(figsize=(5.5,3.8))
        plt.scatter(selectRgxChange, RgxDFchange["geneFC"], s=pointsize, facecolors='none', edgecolors='grey', alpha=0.2,
                    label="All site-to-gene")

        if selectdiffgene:
            plt.scatter(selectRgxChange[Rgx_selectgene_bool], RgxDFchange["geneFC"][Rgx_selectgene_bool], s=pointsize[Rgx_selectgene_bool],
                    color="#d9fa45", alpha=0.2, label="Selected_by_gene", )

        if selectRgxQuantile:
            plt.scatter(selectRgxChange[changeRgxBool], RgxDFchange["geneFC"][changeRgxBool], s=pointsize[changeRgxBool],
                    color="#1abc9c", alpha=0.2, label="Rgx_quantile", )

        if selectRgxRatio:
            plt.scatter(selectRgxChange[majorRgxBool], RgxDFchange["geneFC"][majorRgxBool], s=pointsize[majorRgxBool],
                    color="#52c22d", alpha=0.2, label="Major_Rgx", )
        
        if selectRgxFold:
            plt.scatter(selectRgxChange[selectRgxFcBool], RgxDFchange["geneFC"][selectRgxFcBool], s=pointsize[selectRgxFcBool],
                    color="#f39c12", alpha=0.2, label="Rgx_Foldchange", )
            
        if selectRgxRank in ["same","diff"]:
            plt.scatter(selectRgxChange[rankBool], RgxDFchange["geneFC"][rankBool], s=pointsize[rankBool],
                    color="#3498db", alpha=0.2, label="Rankbool", marker="x")

        plt.scatter(selectRgxChange[rgxTotalBool], RgxDFchange["geneFC"][rgxTotalBool], s=pointsize[rgxTotalBool], color="#e74c3c", alpha=0.8, 
                    edgecolors='black', facecolors='m', label="Final site-to-gene", )
        if plotxlim:plt.xlim(plotxlim)
        plt.ylabel("Gene logFC")
        plt.title(title)
        if plotxlabel: plt.xlabel(plotxlabel)
        plt.legend(bbox_to_anchor=(1, 0.5),loc='center left',)
        plt.tight_layout()
        if self.outdir:
            #plt.savefig(self.outdir+"/"+title+"_select"+selectRgxType+"_diffs2g_full.pdf")
            plt.savefig(self.outdir+"/"+title+"_select"+selectRgxType+"_diffs2g_full.png",dpi=300)
        
        self.rgxTotalBool = rgxTotalBool
        self.finalep = RgxDFchange[rgxTotalBool]