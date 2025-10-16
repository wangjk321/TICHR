import pandas as pd
import numpy as np
import pyranges as pr
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from scipy.stats import wilcoxon,ttest_1samp
from statsmodels.stats.multitest import multipletests


def merge_and_average(df1, df2, column_index):
    if df1.shape != df2.shape:
        raise ValueError("两个 DataFrame 的形状不匹配，请检查。")
    new_df = df1.copy()
    new_df.iloc[:, column_index] = (df1.iloc[:, column_index] + df2.iloc[:, column_index]) / 2
    return new_df

def bed_intersect(df1: pd.DataFrame, df2: pd.DataFrame,mode: str = 'intersect') -> pd.DataFrame:
    """
    对两个 DataFrame 做 bedtools 样式的 intersect 操作。
    前三列必须是 BED 格式 (chr/start/end)，会自动重命名为 Chromosome/Start/End。
    
    返回一个包含交集的新 DataFrame。
    """
    # 重命名前3列为标准BED格式
    df1 = df1.rename(columns={
        df1.columns[0]: "Chromosome",
        df1.columns[1]: "Start",
        df1.columns[2]: "End"
    })
    df2 = df2.rename(columns={
        df2.columns[0]: "Chromosome",
        df2.columns[1]: "Start",
        df2.columns[2]: "End"
    })

    # 转换 Start/End 为整数
    df1["Start"] = df1["Start"].astype(int)
    df1["End"] = df1["End"].astype(int)
    df2["Start"] = df2["Start"].astype(int)
    df2["End"] = df2["End"].astype(int)

    # 转换为 PyRanges
    gr1 = pr.PyRanges(df1)
    gr2 = pr.PyRanges(df2)

    # 执行操作
    if mode == 'subtract':
        result = gr1.subtract(gr2)
    elif mode == 'intersect':
        result = gr1.intersect(gr2)
    elif mode == "intersectV":
        result = gr1.overlap(gr2, invert=True)
    else:
        raise ValueError("mode 只能是 'subtract' 或 'intersect'")

    # 返回结果 DataFrame
    return result.df



def plotFDRbox(selectp,otherp):
    plt.figure(figsize=(3, 3))

    box = plt.boxplot(
        [selectp, otherp],
        patch_artist=True,
        widths=0.6,
        boxprops=dict(facecolor="#26A69A", color="#333333",alpha=0.5),  # 青绿色
        medianprops=dict(color="#FF6F61", linewidth=2),       # 珊瑚粉 median
        whiskerprops=dict(color="#999999", linestyle="--"),
        capprops=dict(color="#999999"),
        flierprops=dict(marker='o', color='gray', alpha=0.3)
    )

    plt.xticks([1, 2], ['DEGs', 'Other genes'], fontsize=11)
    plt.ylabel('-log10(qvalue)', fontsize=11)
    ymax = max(max(selectp), max(otherp))
    plt.ylim(0, ymax * 1.1)
    # 去掉右边和顶边框线
    for spine in ['right', 'top']:
        plt.gca().spines[spine].set_visible(False)
    # 统计检验
    stat, pval = mannwhitneyu(selectp, otherp, alternative='two-sided')
    # 写入 p 值
    plt.text(1.5, max(max(selectp), max(otherp)) * 1.05,f'p = {pval:.2e}', ha='center', fontsize=10)
    plt.title("in silico deletion FDR", fontsize=12)
    plt.tight_layout()
    plt.show()


def silico(RgxDF_Ctrl,deletepeakDF,RgDF_Ctrl=None,degtype="all",degfdr=0.01,nrandom=20,fixrandom=False):
    beforeDeletionRg = RgxDF_Ctrl.groupby(4)[11].sum()
    afterDeletion = bed_intersect(RgxDF_Ctrl,deletepeakDF,mode= 'intersectV')
    afterDeletionRg = afterDeletion.groupby(4)[11].sum()
    afterDeletionRgReindex = afterDeletionRg.reindex(beforeDeletionRg.index, fill_value=0) 
    diffrg = abs(afterDeletionRgReindex-beforeDeletionRg)
    
    rng = np.random.default_rng(777)
    diffrg_random_list=[]
    for i in range(nrandom):
        seed_i = rng.integers(0, 1e9) 
        print(i)
        randomRgx = RgxDF_Ctrl.sample(n=len(afterDeletion), random_state=seed_i)
        randomRg = randomRgx.groupby(4)[11].sum()
        randomRgReindex = randomRg.reindex(beforeDeletionRg.index, fill_value=0)
        diffrg_random = abs(randomRgReindex-beforeDeletionRg)
        diffrg_random_list.append(diffrg_random)
    
    diffrg_random_DF = pd.DataFrame(diffrg_random_list).T
    wilcoxplist=[]
    for index, eachrow in diffrg_random_DF.iterrows():
        randomdiffi = abs(eachrow)
        selectdiffi = abs(diffrg[index])
        if np.all(randomdiffi == 0):
            print("差值全为 0，Wilcoxon 检验不适用。")
            _, wilcoxp = None, 1.0
        else:
            _,wilcoxp = wilcoxon([selectdiffi - x for x in randomdiffi],alternative="greater") #选中的差异比随机值大
        wilcoxplist.append(wilcoxp)

    _, fdr_pvals, _, _ = multipletests(wilcoxplist, method='fdr_bh')
    qvaluelist= pd.Series(fdr_pvals).fillna(1)
    qvaluelist.index= diffrg.index
    
    if RgDF_Ctrl is not None:
        RgDF_Ctrl.index=RgDF_Ctrl[3]
        RgDF_Ctrl['qvalue'] = qvaluelist
        RgDF_Ctrl['qvalue'] = RgDF_Ctrl['qvalue'].fillna(1)
        
        if degtype== 'all':
            degbool = (abs(RgDF_Ctrl[6])>1) & (RgDF_Ctrl[7]<degfdr)
        elif degtype == 'down':
            degbool = (RgDF_Ctrl[6]<-1) & (RgDF_Ctrl[7]<degfdr)
        elif degtype == 'up':
            degbool = (RgDF_Ctrl[6]>1) & (RgDF_Ctrl[7]<degfdr)
        selectp = RgDF_Ctrl[degbool]['qvalue']
        otherp = RgDF_Ctrl[~degbool]['qvalue']
        print(otherp)

        print("真实DEGs中受干扰后判定为靶基因的比例:",(selectp<0.05).mean())
        print("非DEGs中受干扰后判定为靶基因的比例:",(otherp<0.05).mean())
        plotFDRbox(-np.log10(selectp),-np.log10(otherp))
    
    return(qvaluelist)
    '''pltoneprc(degbool,-np.log10(RgDF_Ctrl['qvalue'].fillna(1)),"auprc",cols[0])
    plt.legend()
    plt.show()
    
    logfc = RgDF_Ctrl[6]
    fdr = -np.log10(RgDF_Ctrl[7])
    deletep = -np.log10(RgDF_Ctrl['qvalue'])
    plt.scatter(logfc,fdr,alpha=0.5,c=deletep,s=deletep,cmap="Greens")
    plt.xlim(-5,5)
    plt.colorbar()'''