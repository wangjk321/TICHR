import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr


def corrfunc(resultdf,pair_list,corrtype="pearson",suffix="_RgDf.tsv"):
    corrrg_list = []
    corrdelta_list = []
    for pair in pair_list:
        print(pair)
        #chip = pair[0].split("_")[0]+"_"+pair[0].split("_")[3]
        RgDF_Ctrl = pd.read_csv(resultdf+pair[1]+suffix,header=None,sep="\t")
        RgDF_Treat = pd.read_csv(resultdf+pair[0]+suffix,header=None,sep="\t")

        meanRg = (np.log1p(RgDF_Treat[9]) + np.log1p(RgDF_Ctrl[9]))/2
        meanTPM = (RgDF_Ctrl[8] + RgDF_Treat[8])/2
        deltaRg = np.log1p(RgDF_Treat[9]) -  np.log1p(RgDF_Ctrl[9])
        deltaTPM = RgDF_Ctrl[6]

        if corrtype == "pearson":
            corrrg,_ = pearsonr(meanRg, meanTPM)
            corrdelta, _ = pearsonr(deltaTPM, deltaRg)
        elif corrtype == "spearman":
            corrrg,_ = spearmanr(meanRg, meanTPM)
            corrdelta, _ = spearmanr(deltaTPM, deltaRg)
            
        corrrg_list.append(corrrg)
        corrdelta_list.append(corrdelta)
        
    return(corrrg_list,corrdelta_list)


from adjustText import adjust_text
import matplotlib.pyplot as plt

def plot_correlation_scatter(corrrg_list, corrdelta_list,factorlist,corrtype="pearson",title="title",outname="correlation_scatter.pdf"):
    plt.figure(figsize=(4, 4))
    
    # 散点图
    plt.scatter(corrrg_list, corrdelta_list, color="#AB47BC", s=100, edgecolor="black", linewidth=0.8)

    # 添加文本标签
    texts = []
    for i, wtype in enumerate(factorlist):
        x, y = corrrg_list[i], corrdelta_list[i]
        texts.append(plt.text(x, y, wtype, fontsize=9))

    # 自动调整标签位置
    adjust_text(
        texts,
        only_move={'points': 'y', 'text': 'xy'},
        arrowprops=dict(arrowstyle='-', color='gray', lw=0.6, shrinkA=5, shrinkB=5),
        expand_points=(1.2, 1.2),
        expand_text=(1.2, 1.2),
        force_text=0.5,
        force_points=0.3,
    )

    plt.xlabel("Corr(Rg vs TPM)")
    plt.ylabel("Corr(\u0394Rg vs gene-logFC)")
    plt.title(title)
    plt.grid(True, linestyle="--", alpha=0.5)
    plt.axhline(0, color="gray", linestyle="--", linewidth=1)
    plt.axvline(0, color="gray", linestyle="--", linewidth=1)
    plt.tight_layout()
    plt.savefig(outname)