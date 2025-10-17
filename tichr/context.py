import pandas as pd
import numpy as np
from matplotlib.gridspec import GridSpec
from scipy.stats import rankdata
from scipy.stats import spearmanr,pearsonr
from rpy2.robjects import r
from rpy2.robjects.packages import importr
from scipy.stats import mannwhitneyu
from .RPpredictDEG import *



def logneg(values):
    return np.sign(values) * np.log1p(np.abs(values))

def makediffrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1-rank2)/len(list1)

def makesumrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1+rank2)/(2*len(list1))

def makesumrank_center0(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1+rank2)/(len(list1)) -1

def mergeDF(rgCtrl_path,rgTreat_path,rgxCtrl_path,rgxTreat_path,samehic=False, minRgx=0.1,minRgxRatio=0.01,filter=True):
    rgCtrl = pd.read_csv(rgCtrl_path, sep="\t",header=None)
    rgTreat = pd.read_csv(rgTreat_path, sep="\t",header=None)
    rgxCtrl = pd.read_csv(rgxCtrl_path, sep="\t",header=None)
    rgxTreat = pd.read_csv(rgxTreat_path, sep="\t",header=None)
    
    if filter:
        goodvalue = rgxCtrl[11] > float('-inf')
        goodvalue = ((rgxCtrl[11]+rgxTreat[11])/2 > minRgx) & goodvalue
        goodvalue = ((rgxCtrl[12]+rgxTreat[12])/2 > minRgxRatio) & goodvalue

        #把真正的ABC换成两个百分比的平均值
        rgxCtrl[12] = (rgxCtrl[12]+rgxTreat[12])/2
        rgxTreat[12] = (rgxCtrl[12]+rgxTreat[12])/2
        
        if samehic:
            meanhic = (rgxCtrl[10]+rgxTreat[10])/2
            rgxCtrl[10] = meanhic
            rgxCtrl[11] = meanhic * rgxCtrl[3]
            
            rgxTreat[10] = meanhic
            rgxTreat[11] = meanhic * rgxTreat[3]
        
        rgxCtrl = rgxCtrl[goodvalue]
        rgxCtrl.reset_index(drop=True, inplace=True)
        rgCtrl.index = rgCtrl[3]
        rgCtrl[9] = rgxCtrl.groupby(4)[11].sum()
        rgCtrl= rgCtrl.dropna(subset=[rgCtrl.columns[9]])
        rgCtrl.reset_index(drop=True, inplace=True)
        
        rgxTreat = rgxTreat[goodvalue]
        rgxTreat.reset_index(drop=True, inplace=True)
        rgTreat.index = rgTreat[3]
        rgTreat[9] = rgxTreat.groupby(4)[11].sum()
        rgTreat = rgTreat.dropna(subset=[rgTreat.columns[9]])
        rgTreat.reset_index(drop=True, inplace=True)


    # 拼接最后一列（rgTreat）到 rgCtrl
    rg_merged = pd.concat([rgCtrl, rgTreat.iloc[:, -1]], axis=1, ignore_index=True)
    # 拼接第12列（rgxTreat，第11号索引）到 rgxCtrl
    rgx_merged = pd.concat([rgxCtrl, rgxTreat.iloc[:, 11]], axis=1, ignore_index=True)
    
    return(rg_merged,rgx_merged)


def prepare_select_by_rank(mergedfile, basedon = "rg",title="title",
                           label="TF",genelabel="(all genes)",filetype="file",
                           negative_cutoff=0.8,positive_cutoff=0.8,
                           rgctrl_col_num=9, rgtreat_col_num=10,tpm_col_num=8, logfc_col_num=6,
                           negshow="diffrank",negshowabs=False,posshow="sumrank",posshowabs=False,
                           outname=None,plotSelect=False,plotscatter=True,
                           negselect="diffrank",negselectabs=False,posselect="sumrank",posselectabs=False,
                           outprefix="Sample"
                           ):
    if filetype=="file":
        mergedDF = pd.read_csv(mergedfile,header=None,sep="\t")
    elif filetype=="pandas":
        mergedDF = mergedfile.copy()

    #mergedDF = mergedDF[mergedDF[rgtreat_col_num]+mergedDF[rgctrl_col_num] != 0]
    mergedDF.reset_index(drop=True, inplace=True)

    meanRg = (np.log1p(mergedDF[rgtreat_col_num]) + np.log1p(mergedDF[rgctrl_col_num])) / 2
    meanTPM = mergedDF[tpm_col_num] #已经log后了
    changeRg = np.log1p(mergedDF[rgtreat_col_num]) - np.log1p(mergedDF[rgctrl_col_num])
    changeTPM = mergedDF[logfc_col_num] #已经处理好了

    if plotscatter:
                # 创建子图
        fig, axes = plt.subplots(2, 1, figsize=(3, 6))  # 上下两图
        plt.subplots_adjust(hspace=0.4)

        # -------- 图1：meanRg vs meanTPM --------
        sns.regplot(
            x=meanRg,
            y=meanTPM,
            ax=axes[0],
            scatter_kws={'alpha': 0.5, 's': 5, 'color': "slateblue"},
            line_kws={'color': 'grey'}
        )
        corr1 = np.corrcoef(meanRg, meanTPM)[0, 1]
        axes[0].set_title(title, fontsize=12)
        axes[0].set_xlabel(label+"-Rg")
        axes[0].set_ylabel("Gene TPM")
        axes[0].text(0.05, 0.95, f'Corr(Rg vs TPM): {corr1:.3f}',
                    transform=axes[0].transAxes, fontsize=10, va='top')

        # -------- 图2：deltaRg vs deltaTPM --------
        sns.regplot(
            x=changeRg,
            y=changeTPM,
            ax=axes[1],
            scatter_kws={'alpha': 0.5, 's': 5, 'color': "olive"},
            line_kws={'color':   'grey'}
        )
        corr2 = np.corrcoef(changeRg, changeTPM)[0, 1]
        axes[1].set_title("")
        axes[1].set_xlabel(label+"-ΔRg")
        axes[1].set_ylabel("Gene logFC")
        axes[1].text(0.05, 0.95, f'Corr(\u0394Rg vs logFC): {corr2:.3f}',
                    transform=axes[1].transAxes, fontsize=10, va='top')

        # 保存
        plt.tight_layout()
        plt.savefig(f"{outprefix}_scatter.pdf")

    '''
    plot_scatter_with_fit(meanRg,meanTPM,s=5,title=title,
                     x_label="Rg(log)"+label, y_label ="TPM "+genelabel)
    #plt.savefig(label+"_meanRg_corr_meanTPM.pdf")
    
    plot_scatter_with_fit(changeRg,changeTPM,title=title,s=3,color="olive",
                     x_label="Rg changes (logFC)"+ label, y_label ="TPM changes (logFC)"+genelabel)
    #plt.savefig(label+"_changeRg_corr_changeTPM.pdf")
    '''

    if basedon == "rg":
        select_by_rank(meanRg,meanTPM,changeRg,changeTPM,
                       labellist=[label+" Rg","Gene TPM"], negative_cutoff=negative_cutoff,positive_cutoff=positive_cutoff,
                       select_label="Rg-TPM",predict_label="\u0394Rg-logFC(TPM)",
                       negshow=negshow,negshowabs=negshowabs,posshow=posshow,posshowabs=posshowabs,
                       outname=outname,plotSelect=plotSelect,
                       negselect=negselect,negselectabs=negselectabs,posselect=posselect,posselectabs=posselectabs,outprefix=outprefix)
    elif basedon == "deltarg":
        select_by_rank(changeRg,changeTPM,meanRg,meanTPM,
                       labellist=[label+" \u0394Rg","Gene logFC"], negative_cutoff=negative_cutoff,positive_cutoff=positive_cutoff,
                       predict_label="Rg-TPM",select_label="\u0394Rg-logFC(TPM)",
                       negshow=negshow,negshowabs=negshowabs,posshow=posshow,posshowabs=posshowabs,
                        outname=outname,plotSelect=plotSelect,
                       negselect=negselect,negselectabs=negselectabs,posselect=posselect,posselectabs=posselectabs,outprefix=outprefix)
        
# select_rg, 所选的Rg值，比如ctrl和treat平均的Rg值
# select_tpm, 所选的基因表达值，比如ctrl和treat平均的TPM值
# predict_rg, Rg的变化值，比如ctrl和treat的Rg的logFC
# predict_tpm, tpm的变化值,比如ctrl和treat的tpm的logFC
def select_by_rank(select_rg,select_tpm,predict_rg,predict_tpm,
                   negative_cutoff=0.8,positive_cutoff=0.8,
                   labellist=["Rg","Gene TPM"],
                   select_label="Rg-geneTPM", 
                   predict_label="\u0394Rg-geneFC",plotSelect=False,
                   negselect="diffrank",negselectabs=False,posselect="sumrank",posselectabs=False,
                   negshowabs=False,
                   negshow="diffrank",posshow="sumrank",posshowabs=False,
                   outname=None,outprefix=None):
    rank_select_rg = select_rg.rank(ascending= False)
    rank_select_tpm= select_tpm.rank(ascending= False)
    diffrank_select = (rank_select_tpm-rank_select_rg)/len(rank_select_tpm)
    sumrank_select = 1- (rank_select_tpm+rank_select_rg)/(len(rank_select_tpm)*2)
    sumrank0_select= -makesumrank_center0(rank_select_tpm,rank_select_rg)

    
    rank_predict_rg = predict_rg.rank(ascending= False)
    rank_predict_tpm= predict_tpm.rank(ascending= False)
    diffrank_predict = (rank_predict_tpm-rank_predict_rg)/len(rank_predict_tpm)
    sumrank_predict = 1 - (rank_predict_tpm+rank_predict_rg)/(len(rank_predict_tpm)*2)
    sumrank0_predict = -makesumrank_center0(rank_predict_tpm,rank_predict_rg)
    
    #1. 根据diffrank选负调控
    if negselect=="diffrank":
        negselect_score = diffrank_select
    elif negselect=="sumrank":
        negselect_score = sumrank_select
    elif negselect=="sumrank0":
        negselect_score = sumrank0_select

    if negselectabs:
        negselect_score = abs(negselect_score)

    print(negselect_score.describe())
    value_bool = negselect_score > negative_cutoff
    
    
    #fig, axes = plt.subplots(1, 3, figsize=(9, 4))  # 1 行 2 列的子图布局
    fig = plt.figure(figsize=(8, 3.8))
    gs = GridSpec(1, 3, width_ratios=[1.5, 1, 1])

    colors = negselect_score
    axes0 = fig.add_subplot(gs[0])
    sc1 = axes0.scatter(select_rg,select_tpm,s=5,c=colors,cmap="coolwarm")
    if plotSelect:
        axes0.scatter(select_rg[value_bool],
                select_tpm[value_bool],facecolors='none', edgecolors='grey',s=6,linewidth=0.5)
    axes0.set_title('Select genes by\ndiffrank ('+select_label+')')
    axes0.set_xlabel(labellist[0],fontsize=12)
    axes0.set_ylabel(labellist[1],fontsize=12)
    cbar1 = fig.colorbar(sc1, ax=axes0,shrink=0.5)
    cbar1.ax.set_xlabel('diffRank')
    cbar1.ax.xaxis.set_label_position('bottom')
    cbar1.ax.xaxis.label.set_rotation(90)
    
    # 在第一个子图上画第一个箱线图
    colors = ['snow', 'rosybrown']
    axes1 = fig.add_subplot(gs[1])
    bp=axes1.boxplot([negselect_score, negselect_score[value_bool]],widths=0.8,notch=True,patch_artist=True)
    for box, color in zip(bp['boxes'], colors):
        box.set_facecolor(color)
    axes1.set_title('diffrank\n('+select_label+')')  # 设置标题
    axes1.set_xticklabels(["All genes", "'Negatively'\nregulated\ngenes"])
    plt.setp(axes1.get_xticklabels(), rotation=45, ha="right")
    
    # 在第二个子图上画第二个箱线图
    axes2 = fig.add_subplot(gs[2])
    if negshow == "diffrank":
        plot_predict = diffrank_predict
    elif negshow == "sumrank":
        plot_predict = sumrank_predict
    elif negshow == "sumrank0":
        plot_predict = sumrank0_predict
    
    
    if negshowabs: plot_predict=abs(plot_predict)

    bp=axes2.boxplot([plot_predict, plot_predict[value_bool]],widths=0.8,notch=True,patch_artist=True)
    _, pvalue = mannwhitneyu(plot_predict,plot_predict[value_bool], alternative='less')
    axes2.text(1.5, 0.95,f"p = {pvalue:.2e}",ha='center', fontsize=10,color="forestgreen")

    colors = ['snow', 'darkviolet']
    for box, color in zip(bp['boxes'], colors):
        box.set_facecolor(color)
    axes2.set_title('diffrank\n('+predict_label+')')  # 设置标题
    axes2.set_xticklabels(["All genes", "'Negatively'\nregulated\ngenes"])
    plt.setp(axes2.get_xticklabels(), rotation=45, ha="right")

    # 调整布局，防止子图重叠
    plt.tight_layout()
    # 显示子图
    if outname:
        plt.savefig(outname+"_negative.pdf")
        plt.savefig(outname+"_negative.png",dpi=300)
    plt.show()
    plt.savefig(f"{outprefix}_diffrank_top_correlated.pdf")
    #plt.savefig(outname+"_diffrank.pdf")
    
    ##############2.根据sumrank选择positive regulation##########

    if posselect=="diffrank":
        posselect_score = diffrank_select
    elif posselect=="sumrank":
        posselect_score = sumrank_select
    elif posselect=="sumrank0":
        print(posselect)
        posselect_score = sumrank0_select
    
    if posselectabs:
        posselect_score = abs(posselect_score)

    print(posselect_score.describe())
    value_bool2 = posselect_score > positive_cutoff

    fig = plt.figure(figsize=(8, 3.8))
    gs = GridSpec(1, 3, width_ratios=[1.5, 1, 1])
    
    ax0 = fig.add_subplot(gs[0])
    ax0.set_title('Select genes by\nsumrank ('+select_label+')')
    ax0.set_xlabel(labellist[0],fontsize=12)
    ax0.set_ylabel(labellist[1],fontsize=12)
    colors = posselect_score
    if not posselectabs and posselect in ["diffrank","sumrank0"]:
        sc2 =ax0.scatter(select_rg,select_tpm,s=5,c=colors,cmap="coolwarm")
    else:
        sc2 =ax0.scatter(select_rg,select_tpm,s=5,c=colors,cmap="Purples")
    if plotSelect:
        ax0.scatter(select_rg[value_bool2],
                select_tpm[value_bool2],facecolors='none', edgecolors='grey',s=6,linewidth=0.5)
    cbar2 = fig.colorbar(sc2, ax=ax0,shrink=0.5)
    cbar2.ax.set_xlabel('sumRank')
    cbar2.ax.xaxis.set_label_position('bottom')
    cbar2.ax.xaxis.label.set_rotation(90)
    
    colors = ['snow', 'rosybrown']  # 定义填充颜色 
    ax1 = fig.add_subplot(gs[1])
    bp1 = ax1.boxplot([posselect_score, posselect_score[value_bool2]],widths=0.8,notch=True,patch_artist=True)
    for box, color in zip(bp1['boxes'], colors):
        box.set_facecolor(color)  # 设置填充颜色
        
    ax1.set_title('sumrank\n('+select_label+')')  # 设置标题
    ax1.set_xticklabels(["All genes", "'Positively'\nregulated\ngenes"])
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right")

    if posshow == "diffrank":
        plot_predict = diffrank_predict
    elif posshow == "sumrank":
        plot_predict = sumrank_predict
    elif posshow == "sumrank0":
        plot_predict = sumrank0_predict
    
    if posshowabs: plot_predict=abs(plot_predict)
    
    ax2 = fig.add_subplot(gs[2])
    bp2 = ax2.boxplot([plot_predict, plot_predict[value_bool2]],widths=0.8,notch=True,patch_artist=True)
    _, pvalue = mannwhitneyu(plot_predict,plot_predict[value_bool2], alternative='less')
    ax2.text(1.5, 0.95,f"p = {pvalue:.2e}",ha='center', fontsize=10,color="forestgreen")

    colors = ['snow', 'darkviolet']
    for box, color in zip(bp2['boxes'], colors):
        box.set_facecolor(color)  # 设置填充颜色
    ax2.set_title('sumrank\n('+predict_label+')')  # 设置标题
    ax2.set_xticklabels(["All genes", "'Positively'\nregulated\ngenes"])
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right")
    
    # 调整布局，防止子图重叠
    plt.tight_layout()
    if outname:
        plt.savefig(outname+"_positive.pdf")
        plt.savefig(outname+"_positive.png",dpi=300)
    plt.show()
    plt.savefig(f"{outprefix}_sumrank_top_correlated.pdf")
    #plt.savefig(outname+"_sumrank.pdf")






def logneg(values):
    return np.sign(values) * np.log1p(np.abs(values))

def makediffrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1-rank2)/len(list1)

def makesumrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1+rank2)/(2*len(list1))


def extractNeg(mergedRgFile, mergedRgxFile, rgctrl_col_num=9, rgtreat_col_num=10,
               tpm_col_num=8, logfc_col_num=6, geneID_col_num=3,geneFDR_col_num=7, 
               rgx_geneID_col=4,rgx_ctrl_col=11,rgx_treat_col=13,negboolRgX=None,iteration=True,outname="outname",
               showInteration=False,iteration_count=0,outdir="identify_context",filetype="file",corrtype="pearson",
               extractType="negative",minRgxRatio=0,epifdr=False,
               geneFDR_cutoff=0.1,geneFC_cutoff=0.5,rgxFC_cutoff=0.2,outprefix="Sample"):
    if filetype=="file":
        mergedRg = pd.read_csv(mergedRgFile,header=None,sep="\t")
        mergedRgx = pd.read_csv(mergedRgxFile,header=None,sep="\t")
    elif filetype=="pandas":
        mergedRg = mergedRgFile.copy()
        mergedRgx = mergedRgxFile.copy()

    mergedRg = mergedRg[mergedRg[rgtreat_col_num]+mergedRg[rgctrl_col_num] != 0]
    mergedRg.reset_index(drop=True, inplace=True)
    mergedRgx = mergedRgx[mergedRgx[rgx_ctrl_col]+mergedRgx[rgx_treat_col] !=0 ]
    mergedRgx.reset_index(drop=True, inplace=True)

    df_rg = pd.DataFrame()
    df_rg["Rg_ctrl"] = np.log1p(np.array(mergedRg.iloc[:,rgctrl_col_num]))
    df_rg["Rg_treat"] = np.log1p(np.array(mergedRg.iloc[:,rgtreat_col_num]))
    df_rg["meanRg"] = (df_rg["Rg_treat"] + df_rg["Rg_ctrl"]) / 2
    df_rg['logFC-Rg'] = df_rg["Rg_treat"] - df_rg["Rg_ctrl"]
    df_rg["TPM_ctrl"] = np.array(mergedRg.iloc[:,tpm_col_num])
    df_rg["TPM_treat"] =np.array(mergedRg.iloc[:,tpm_col_num])
    df_rg["meanTPM"] = (df_rg["TPM_ctrl"] +  df_rg["TPM_treat"]) /2
    df_rg['logFC-CPM'] = mergedRg.iloc[:,logfc_col_num]
    df_rg.index = np.array(mergedRg.iloc[:,geneID_col_num])
    
    # map gene information to Rgx
    geneFCdict = mergedRg.set_index(geneID_col_num)[logfc_col_num].to_dict()
    geneFDRdict = mergedRg.set_index(geneID_col_num)[geneFDR_col_num].to_dict()
    geneTPMdict =  mergedRg.set_index(geneID_col_num)[tpm_col_num].to_dict()
    Rgx_order_FC = mergedRgx[rgx_geneID_col].map(geneFCdict)
    Rgx_order_FDR = mergedRgx[rgx_geneID_col].map(geneFDRdict)
    Rgx_order_TPM = mergedRgx[rgx_geneID_col].map(geneTPMdict)

    df_rgx = pd.DataFrame()
    df_rgx["Rgx_ctrl"] = mergedRgx[rgx_ctrl_col]
    df_rgx["Rgx_treat"] =  mergedRgx[rgx_treat_col]
    df_rgx["meanRgx"] = np.log1p((df_rgx["Rgx_ctrl"]+df_rgx["Rgx_treat"])/2)
    df_rgx["geneTPM"] = Rgx_order_TPM
    df_rgx["geneFC"] = Rgx_order_FC
    df_rgx["logFC_Rgx"] = np.log1p(df_rgx["Rgx_treat"]) - np.log1p(df_rgx["Rgx_ctrl"])
    df_rgx["geneID"] =  mergedRgx[rgx_geneID_col]
    if epifdr: df_rgx["epifdr"] = mergedRgx[14]
    
    if iteration_count == 0 and showInteration:
        plt.figure(figsize=(3, 3))
        plt.scatter(df_rgx["logFC_Rgx"],df_rgx["geneFC"],s=5,alpha=0.5,color="grey")
        plt.title("Changes of links",fontsize=14)
        plt.xlabel("\u0394RgX")
        plt.ylabel("gene-logFC")
        plt.tight_layout()
        plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRgx-geneFC.pdf")
        plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRgx-geneFC.png",dpi=300)

        plt.figure(figsize=(8.5, 3))
        plt.subplot(1, 3, 3)
        plot_scatter_with_fit(df_rg["meanRg"],df_rg["meanTPM"],insideplot=True,
                                title="Genes (Rg vs TPM)",x_label="Initial Rg", y_label ="Gene-TPM",s=5)
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.subplot(1, 3, 2)
        plot_scatter_with_fit(df_rg['logFC-Rg'],df_rg['logFC-CPM'],insideplot=True,
                      title="Changes of genes",s=5,color="olive",
                      x_label="\u0394Rg", y_label ="Gene-logFC")
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.subplot(1, 3, 1)
        plt.scatter(df_rgx["meanRgx"],df_rgx["geneTPM"],s=5,alpha=0.5,color="grey")
        plt.title("Site-to-gene links",fontsize=14)
        plt.xlabel("\u0394RgX")
        plt.ylabel("gene-logFC")
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.tight_layout()
        plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRg-geneFC.pdf")
        # plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRg-geneFC.png",dpi=300)
        
    #sumrank = makesumrank(df_rgx["logFC_Rgx"],df_rgx["geneFC"])
    sumrank0 = makesumrank_center0(df_rgx["logFC_Rgx"],df_rgx["geneFC"])
    diffrank = makediffrank(df_rgx["logFC_Rgx"],df_rgx["geneFC"])
    
    # define negtive E-P links
    tmpvalue = df_rgx["logFC_Rgx"] * df_rgx["geneFC"]
    degbool = (Rgx_order_FDR < geneFDR_cutoff) & (abs(Rgx_order_FC)>geneFC_cutoff)
    if negboolRgX is None:
        if extractType == "negative":
            negboolRgX = degbool & (abs(df_rgx["logFC_Rgx"])>rgxFC_cutoff) & (tmpvalue < 0)
            if epifdr: negboolRgX = negboolRgX & (df_rgx["epifdr"]<0.05)
        elif extractType == "positive":
            negboolRgX = degbool & (abs(df_rgx["logFC_Rgx"])>rgxFC_cutoff) & (tmpvalue > 0)
            if epifdr: negboolRgX = negboolRgX & (df_rgx["epifdr"]<0.05)
    
    print(f"----------Iteration {iteration_count} --------------")
    iteration_count += 1 

    #make adjust Rgx
    if extractType == "negative":
        neg_weight = np.where(negboolRgX, -1*abs(diffrank) , 1)
    elif extractType == "positive":
        neg_weight = np.where(negboolRgX, 1*abs(sumrank0), 1) 
    EPregulate = neg_weight
    df_rgx["adj_Rgx_ctrl"] = df_rgx["Rgx_ctrl"] * EPregulate
    df_rgx["adj_Rgx_treat"] = df_rgx["Rgx_treat"] * EPregulate
    df_rgx["adj_meanRgx"] = logneg((df_rgx["adj_Rgx_ctrl"] + df_rgx["adj_Rgx_treat"])/2)
    df_rgx["adj_logFC_Rgx"] = logneg(df_rgx["adj_Rgx_treat"]) - logneg(df_rgx["adj_Rgx_ctrl"])

    #make adjust Rg
    genesWithNeg = df_rgx["geneID"][negboolRgX].unique()
    
    # Re-calculate Rg
    adj_Rg_ctrl = df_rgx.groupby("geneID")["adj_Rgx_ctrl"].sum()
    adj_Rg_treat = df_rgx.groupby("geneID")["adj_Rgx_treat"].sum()
    
    adj_logFC_Rg = logneg(adj_Rg_treat) - logneg(adj_Rg_ctrl)
    adj_logFC_gene = df_rgx.groupby("geneID")["geneFC"].mean()

    adj_meanRg = (logneg(adj_Rg_treat) + logneg(adj_Rg_ctrl))/2
    adj_geneTPM = df_rgx.groupby("geneID")["geneTPM"].mean()
    
    adj_RgDF = pd.DataFrame()
    adj_RgDF["adj_Rg_ctrl"] = adj_Rg_ctrl
    adj_RgDF["adj_Rg_treat"] = adj_Rg_treat
    adj_RgDF["adj_logFC_Rg"] = adj_logFC_Rg
    adj_RgDF["adj_logFC_gene"] = adj_logFC_gene
    adj_RgDF["adj_meanRg"] = adj_meanRg
    adj_RgDF["adj_geneTPM"] = adj_geneTPM
    
    adj_RgDF_withNeg = adj_RgDF.loc[genesWithNeg]
    raw_RgDF_withNeg = df_rg.loc[genesWithNeg]
    
    #try a loop for negtive
    rawdiffrank_neggene = makediffrank(df_rg["meanRg"],df_rg["meanTPM"])[genesWithNeg]
    adjdiffrank_neggene = makediffrank(adj_meanRg,adj_geneTPM)[genesWithNeg]
    improvedgene = genesWithNeg[(abs(adjdiffrank_neggene) < abs(rawdiffrank_neggene))]
    improvedgeneratio = np.mean((abs(adjdiffrank_neggene) < abs(rawdiffrank_neggene)))
    improvedgenebool = [i in improvedgene for i in df_rgx["geneID"] ]
    further_negtive = negboolRgX & improvedgenebool
    
    if corrtype == "pearson":
        adj_corr = pearsonr(adj_meanRg,adj_geneTPM)[0]
        raw_corr = pearsonr(df_rg["meanRg"],df_rg["meanTPM"])[0]
    elif corrtype == "spearman":
        adj_corr = spearmanr(adj_meanRg,adj_geneTPM)[0]
        raw_corr = spearmanr(df_rg["meanRg"],df_rg["meanTPM"])[0]
    
    if showInteration or (np.array_equal(further_negtive, negboolRgX) and outdir) or not iteration:
        plt.figure(figsize=(3.5, 3))
        color_values  = np.where(negboolRgX, abs(diffrank),abs(diffrank)[negboolRgX].min())
        plt.scatter(df_rgx["logFC_Rgx"][negboolRgX],df_rgx["geneFC"][negboolRgX],s=3,alpha=1,
                    c=color_values[negboolRgX],cmap="jet")
        plt.colorbar(shrink=0.4)
        plt.scatter(df_rgx["logFC_Rgx"][~negboolRgX],df_rgx["geneFC"][~negboolRgX],s=5,alpha=0.5,
                    color="grey")
        plt.title("Changes of links",fontsize=14)
        plt.xlabel("\u0394RgX")
        plt.ylabel("gene-logFC")
        plt.tight_layout()
        if outdir:
            plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRgx-geneFC.pdf")
            # plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRgx-geneFC.png",dpi=300)

        plt.figure(figsize=(8.5, 3))
        plt.subplot(1, 3, 1)
        plt.title("Site-to-gene links",fontsize=14)
        plt.scatter(df_rgx["adj_meanRgx"][negboolRgX],df_rgx["geneTPM"][negboolRgX],s=5,alpha=1,color="m",label="adjusted links")
        plt.scatter(df_rgx["meanRgx"],df_rgx["geneTPM"],s=2,alpha=0.2,color="#689673",label="original links")
        plt.legend(fontsize=8)
        plt.xlabel("Adjusted RgX")
        plt.ylabel("Gene-TPM")
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.subplot(1, 3, 2)
        plt.title("Changes of genes",fontsize=14)
        plt.scatter(adj_logFC_Rg,adj_logFC_gene,color="indigo",s=4,alpha=0.5,label="adjusted genes")
        plt.scatter(df_rg["logFC-Rg"],df_rg["logFC-CPM"],s=4,alpha=0.5,color="#939668",label="original genes")
        plt.legend(fontsize=8)
        plt.xlabel("\u0394Rg")
        plt.ylabel("Gene-logFC")
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        plt.subplot(1, 3, 3)
        plot_scatter_with_fit(adj_meanRg,adj_geneTPM,color="darkorange",insideplot=True,
                                title="Genes (Rg vs TPM)",x_label="Adjusted Rg", y_label ="Gene-TPM",s=5)
        ax = plt.gca()
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        plt.tight_layout()
        if outdir:
            plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRg-geneFC.pdf")
            plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-deltaRg-geneFC.png",dpi=300)

        # 创建主图
        fig, axs = plt.subplots(1, 3, figsize=(6, 3), gridspec_kw={'width_ratios': [1, 1, 1]})

        # Boxplot 1: Gene TPM
        axs[0].boxplot([raw_RgDF_withNeg["TPM_ctrl"], df_rg["TPM_ctrl"]],
                        widths=0.8, notch=True, patch_artist=True, showfliers=False,boxprops=dict(facecolor='lightblue'))
        axs[0].set_title("Gene TPM")
        axs[0].set_ylabel("TPM")
        axs[0].set_xticklabels(['Genes with\nnegative', 'All Genes'],rotation=45,ha="right")
        # Boxplot 2: Gene raw Rg
        axs[1].boxplot([raw_RgDF_withNeg["Rg_ctrl"], df_rg["Rg_ctrl"]],
                        widths=0.8, notch=True, patch_artist=True, showfliers=False,boxprops=dict(facecolor='lightgreen'))
        axs[1].set_title("Gene raw Rg")
        axs[1].set_ylabel("Raw Rg")
        axs[1].set_xticklabels(['Genes with\nnegative', 'All Genes'],rotation=45,ha="right")
        # Boxplot 3: Gene adjusted Rg
        axs[2].boxplot([adj_RgDF_withNeg["adj_meanRg"], adj_RgDF["adj_meanRg"]],
                        widths=0.8, notch=True, patch_artist=True, showfliers=False, boxprops=dict(facecolor='lightcoral'))
        axs[2].set_title("Gene adjusted Rg")
        axs[2].set_ylabel("Adjusted Rg")
        axs[2].set_xticklabels(['Genes with\nnegative', 'All Genes'],rotation=45,ha="right")

        for spine in ['right', 'top']:
            axs[0].spines[spine].set_visible(False)
            axs[1].spines[spine].set_visible(False)
            axs[2].spines[spine].set_visible(False)
            
        plt.tight_layout()
        if outdir:
            plt.savefig(outdir+"/Iteration"+str(iteration_count)+"-adjRg-geneTPM.pdf")
        
    if not np.array_equal(further_negtive, negboolRgX) and adj_corr>raw_corr and iteration:

        if showInteration:
            print(f"{corrtype} correlation - Before: {raw_corr:.4f}, After: {adj_corr:.4f}, Difference: {adj_corr - raw_corr:.4f}")
            print("Ratio of improved genes : ", improvedgeneratio)
            finalRgx=mergedRgx[negboolRgX]
            numNegGene = len(finalRgx[4].unique())
            numAllGene = len(mergedRgx[degbool][4].unique())
            numNegSite = len(finalRgx[[0,1,2]].drop_duplicates())
            numAllSite = len(mergedRgx[degbool][[0,1,2]].drop_duplicates())
            print(outname)
            print("Ratio of negative site-to-genes:",np.mean(negboolRgX)/np.mean(degbool))
            print("Number of genes with negative sites:",numNegGene)
            print("Percentage of genes  with negative sites:",numNegGene/numAllGene)
            print("Number of negative sites:",numNegSite)
            print("Percentage of negative sites for degs:",numNegSite/numAllSite)


        return extractNeg(mergedRgFile, mergedRgxFile, rgctrl_col_num, rgtreat_col_num,tpm_col_num, 
               logfc_col_num, geneID_col_num,geneFDR_col_num, rgx_geneID_col,rgx_ctrl_col,rgx_treat_col,outname=outname,outdir=outdir,
               negboolRgX=further_negtive,showInteration=showInteration,iteration_count=iteration_count,filetype=filetype,
               minRgxRatio=minRgxRatio,epifdr=epifdr,
               geneFDR_cutoff=geneFDR_cutoff,geneFC_cutoff=geneFC_cutoff,rgxFC_cutoff=rgxFC_cutoff)
        
    elif adj_corr<raw_corr and iteration:
        print(outname)
        print("0 genes with 0 negative sites found")
        print(f"{corrtype} correlation - Before: {raw_corr:.4f}, After: {adj_corr:.4f}, Difference: {adj_corr - raw_corr:.4f}")
        return 0,adj_corr - raw_corr
    else:   #final output
        print(f"{corrtype} correlation - Before: {raw_corr:.4f}, After: {adj_corr:.4f}, Difference: {adj_corr - raw_corr:.4f}")

        finalRgx=mergedRgx[negboolRgX]
        if outdir:
            finalRgx.to_csv(outdir+"/"+outname+"_negSite.tsv",sep="\t",header=None, index=None)
            #mergedRg[posboolRgX].to_csv(outdir+"/pos_effect_site.tsv",sep="\t",header=None, index=None)
            finalRgx[finalRgx[13]>=minRgxRatio].to_csv(outdir+"/"+outname+"_negSite_filterRgxRatio.tsv",sep="\t",header=None, index=None)

        numNegGene = len(finalRgx[4].unique())
        numAllGene = len(mergedRgx[degbool][4].unique())
        numNegSite = len(finalRgx[[0,1,2]].drop_duplicates())
        numAllSite = len(mergedRgx[degbool][[0,1,2]].drop_duplicates())
        print(outname)
        print("Ratio of improved genes : ", improvedgeneratio)
        print("Ratio of negative site-to-genes:",np.mean(negboolRgX)/np.mean(degbool))
        print("Number of genes with negative sites:",numNegGene)
        print("Percentage of genes  with negative sites:",numNegGene/numAllGene)
        print("Number of negative sites:",numNegSite)
        print("Percentage of negative sites for degs:",numNegSite/numAllSite)

        print("------------END------------")
        print("---------------------------")
        print("---------------------------")
        print("---------------------------")
        print("---------------------------")
  
        return numNegSite/numAllSite, adj_corr-raw_corr 



def mergeDFmany(rgCtrl_file,rgTreat_file,rgxCtrl_file,rgxTreat_file,
                samehic=False, minRgx=0.1,minRgxRatio=0.01):
    
    rgvalue_ctrl = []
    for i in rgCtrl_file:
        rgi = pd.read_csv(i, sep="\t",header=None)
        rgvalue_ctrl.append(rgi[9])
    rgCtrl = rgi.copy()
    rgCtrl[9] = pd.concat(rgvalue_ctrl, axis=1).mean(axis=1)

    rgvalue_treat = []
    for i in rgTreat_file:
        rgi = pd.read_csv(i, sep="\t", header=None)
        rgvalue_treat.append(rgi[9])
    rgTreat = rgi.copy()
    rgTreat[9] = pd.concat(rgvalue_treat, axis=1).mean(axis=1)

    rgxvalue_ctrl = []
    epi_ctrl = []
    rgxratio_ctrl = []
    hic_ctrl = []
    for i in rgxCtrl_file:
        rgxi = pd.read_csv(i, sep="\t",header=None)
        epi_ctrl.append(rgxi.iloc[:, 3])
        hic_ctrl.append(rgxi.iloc[:, 10])
        rgxvalue_ctrl.append(rgxi.iloc[:, 11])
        rgxratio_ctrl.append(rgxi.iloc[:, 12])
    rgxvalue_ctrl_df = pd.concat(rgxvalue_ctrl, axis=1)
    rgxCtrl =rgxi.copy()
    rgxCtrl[3]= pd.concat(epi_ctrl, axis=1).mean(axis=1)
    rgxCtrl[10]= pd.concat(hic_ctrl, axis=1).mean(axis=1)
    rgxCtrl[11]= pd.concat(rgxvalue_ctrl, axis=1).mean(axis=1)
    rgxCtrl[12]= pd.concat(rgxratio_ctrl, axis=1).mean(axis=1)

    epi_treat = []
    hic_treat = []
    rgxvalue_treat = []
    rgxratio_treat = []
    for i in rgxTreat_file:
        rgxi = pd.read_csv(i, sep="\t",header=None)
        epi_treat.append(rgxi.iloc[:, 3])
        hic_treat.append(rgxi.iloc[:, 10])
        rgxvalue_treat.append(rgxi.iloc[:, 11])
        rgxratio_treat.append(rgxi.iloc[:, 12])
    rgxvalue_treat_df = pd.concat(rgxvalue_treat, axis=1)
    rgxTreat = rgxi.copy() 
    rgxTreat[3] = pd.concat(epi_treat, axis=1).mean(axis=1)
    rgxTreat[10] = pd.concat(hic_treat, axis=1).mean(axis=1)
    rgxTreat[11] = pd.concat(rgxvalue_treat, axis=1).mean(axis=1)
    rgxTreat[12] = pd.concat(rgxratio_treat, axis=1).mean(axis=1)

    goodvalue = rgxCtrl[11] > float('-inf')
    goodvalue = ((rgxCtrl[11]+rgxTreat[11])/2 > minRgx) & goodvalue
    goodvalue = ((rgxCtrl[12]+rgxTreat[12])/2 > minRgxRatio) & goodvalue
    print(sum(goodvalue))
    
    rgxvalue_ctrl_df.columns = [f"ctrl_{i+1}" for i in range(rgxvalue_ctrl_df.shape[1])]
    rgxvalue_treat_df.columns = [f"treat_{i+1}" for i in range(rgxvalue_treat_df.shape[1])]
    counts_df = pd.concat([rgxvalue_ctrl_df, rgxvalue_treat_df], axis=1)[goodvalue]
    counts_df.index = [f"peak_{i}" for i in range(1, len(counts_df) + 1)]
    sample_names = list(rgxvalue_ctrl_df.columns) + list(rgxvalue_treat_df.columns)
    conditions = ["ctrl"] * rgxvalue_ctrl_df.shape[1] + ["treat"] * rgxvalue_treat_df.shape[1]
    meta_df = pd.DataFrame({
        "sample": sample_names,
        "condition": conditions
    })
    # 导出 counts 表（第一列是 peak ID）
    counts_df.insert(0, "peak", counts_df.index)
    counts_df.to_csv("counts.tsv", sep="\t", index=False)
    # 导出 metadata 表
    meta_df.to_csv("metadata.tsv", sep="\t", index=False)
    run_deseq2()

    #把真正的ABC换成两个百分比的平均值
    rgxCtrl[12] = (rgxCtrl[12]+rgxTreat[12])/2
    rgxTreat[12] = (rgxCtrl[12]+rgxTreat[12])/2
    
    if samehic:
        meanhic = (rgxCtrl[10]+rgxTreat[10])/2
        rgxCtrl[10] = meanhic
        rgxCtrl[11] = meanhic * rgxCtrl[3]
        
        rgxTreat[10] = meanhic
        rgxTreat[11] = meanhic * rgxTreat[3]
    
    rgxCtrl = rgxCtrl[goodvalue]
    rgxCtrl.reset_index(drop=True, inplace=True)
    rgCtrl.index = rgCtrl[3]
    rgCtrl[9] = rgxCtrl.groupby(4)[11].sum()
    rgCtrl= rgCtrl.dropna(subset=[rgCtrl.columns[9]])
    rgCtrl.reset_index(drop=True, inplace=True)
    
    rgxTreat = rgxTreat[goodvalue]
    rgxTreat.reset_index(drop=True, inplace=True)
    rgTreat.index = rgTreat[3]
    rgTreat[9] = rgxTreat.groupby(4)[11].sum()
    rgTreat = rgTreat.dropna(subset=[rgTreat.columns[9]])
    rgTreat.reset_index(drop=True, inplace=True)

    # 拼接最后一列（rgTreat）到 rgCtrl
    rg_merged = pd.concat([rgCtrl, rgTreat.iloc[:, -1]], axis=1, ignore_index=True)
    # 拼接第12列（rgxTreat，第11号索引）到 rgxCtrl
    rgx_merged = pd.concat([rgxCtrl, rgxTreat.iloc[:, 11]], axis=1, ignore_index=True)
    
    padj_df = pd.read_csv("deseq2_padj.tsv",sep="\t")
    rgx_merged[14] = padj_df['padj']
    
    return(rg_merged,rgx_merged)


def run_deseq2(counts_file="counts.tsv", metadata_file="metadata.tsv", output_prefix="deseq2"):
    r(f"""
    library(DESeq2)
    library(apeglm)
    library(readr)
    library(BiocParallel)

    # 读取 counts 表
    counts <- read_tsv("{counts_file}")
    counts_mat <- as.matrix(counts[,-1])
    rownames(counts_mat) <- counts[[1]]
    counts_mat <- round(counts_mat)
    storage.mode(counts_mat) <- "integer"

    # 读取 metadata
    meta <- read_tsv("{metadata_file}")
    meta$condition <- factor(meta$condition, levels = c("ctrl", "treat"))

    # 创建 DESeq2 数据集
    dds <- DESeqDataSetFromMatrix(countData = counts_mat,
                                  colData = meta,
                                  design = ~ condition)

    # 差异分析
    dds <- DESeq(dds)

    register(MulticoreParam(32)) 

    # shrink log2 fold changes
    resLFC <- lfcShrink(dds, coef="condition_treat_vs_ctrl", type="apeglm", parallel=TRUE)

    padj_df <- data.frame(peak = rownames(resLFC), padj = resLFC$padj)
    padj_df$padj[is.na(padj_df$padj)] <- 1

    # 保存为 tsv 文件
    write.table(padj_df, file = "{output_prefix}_padj.tsv", row.names = FALSE, quote = FALSE, sep="\\t")
    write.table(resLFC, file = "{output_prefix}_resLFC.tsv", row.names = FALSE, quote = FALSE, sep="\\t")
    """)
