import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from matplotlib import cm
import matplotlib.gridspec as gridspec
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score, classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

from .PRC_ROC import *



def cut_distance(df,rgdf,maxdis=500000,natype="zero"):   #df是RgxDf
    peakpos = (df[1] + df[2]) / 2
    tsspos = np.where(df.iloc[:, 8] == "+", df.iloc[:, 6], df.iloc[:, 7])
    peak2tss = abs(peakpos-tsspos)
    cutdf = df[peak2tss< maxdis]
    cutdf.reset_index(drop=True, inplace=True)
    
    rgdf.index = rgdf[3]
    rgdf[9] = cutdf.groupby(4)[11].sum()
    if natype == "drop":
        cutrgdf= rgdf.dropna(subset=[rgdf.columns[9]])
    elif natype == "zero":
        rgdf[9] = rgdf[9].fillna(0)
        cutrgdf = rgdf
    cutrgdf.reset_index(drop=True, inplace=True)
    return(cutdf,cutrgdf)


def benchmarkvalue(RgDF_Ctrl,RgDF_Treat):
    meanRg = (np.log1p(RgDF_Treat[9]) + np.log1p(RgDF_Ctrl[9]))/2
    meanTPM = (RgDF_Ctrl[8] + RgDF_Treat[8])/2
    
    deltaRg = np.log1p(RgDF_Treat[9]) -  np.log1p(RgDF_Ctrl[9])
    deltaTPM = RgDF_Ctrl[6]
    
    degbool = (abs(RgDF_Ctrl[6])>1) & (RgDF_Ctrl[7]<0.05)
    
    return(abs(deltaRg),degbool)


class manyTFpredict:
    def __init__(self,Rgdir,pair_list,namelist,maxdistance=None,foldchange=False):
        i=0
        for sample in pair_list:
            if maxdistance:
                print(sample)
                RgDF_Ctrl_raw = pd.read_csv(Rgdir+sample[1]+"_RgDf.tsv",header=None,sep="\t")
                RgDF_Treat_raw = pd.read_csv(Rgdir+sample[0]+"_RgDf.tsv",header=None,sep="\t")
                RgxDF_Ctrl_raw = pd.read_csv(Rgdir+sample[1]+"_RgxDf.tsv",header=None,sep="\t")
                RgxDF_Treat_raw = pd.read_csv(Rgdir+sample[0]+"_RgxDf.tsv",header=None,sep="\t")
                
                RgxDF_Ctrl, RgDF_Ctrl = cut_distance(RgxDF_Ctrl_raw,RgDF_Ctrl_raw,maxdis=maxdistance)
                RgxDF_Treat, RgDF_Treat = cut_distance(RgxDF_Treat_raw,RgDF_Treat_raw,maxdis=maxdistance)
            else:
                RgDF_Ctrl = pd.read_csv(Rgdir+sample[1]+"_RgDf.tsv",header=None,sep="\t")
                RgDF_Treat = pd.read_csv(Rgdir+sample[0]+"_RgDf.tsv",header=None,sep="\t")
                
            if i == 0: rgcombine =  RgDF_Ctrl.iloc[:,0:9].copy()
            if foldchange:
                rgcombine[namelist[i]] = np.log1p(RgDF_Treat[9]) -  np.log1p(RgDF_Ctrl[9])
            else:
                rgcombine[namelist[i]] = RgDF_Treat[9] -  RgDF_Ctrl[9]
            i=i+1
        self.rgcombine = rgcombine
        self.namelist = namelist
    
    def makemodel(self,predicttype="all"):
        self.predicttype = predicttype
        #只考虑绝对值
        if predicttype == 'all':
            X = abs(self.rgcombine[self.namelist])
            y = (abs(self.rgcombine[6])>1) & (self.rgcombine[7]<0.05)
        elif predicttype == 'up':
            X = self.rgcombine[self.namelist]
            y = (self.rgcombine[6]>1) & (self.rgcombine[7]<0.05)
        elif predicttype == 'down':
            X = self.rgcombine[self.namelist]
            y = (self.rgcombine[6]<-1) & (self.rgcombine[7]<0.05)
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, 
                                                                                test_size=0.2, random_state=777)
        model = LogisticRegression(penalty='l2')
        model.fit(self.X_train, self.y_train)
        self.model = model
        self.X, self.y = X, y
    
    def plotresult(self,showTF="H3K27ac_rep0",showName="H3K27ac Only",suptitle="suptitle"):
        fig = plt.figure(figsize=(8, 6))
        gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1])
        # 设置左图为当前绘图区域
        ax0 = fig.add_subplot(gs[0, 0])
        plt.sca(ax0)
        pltoneprc(self.y_train, self.model.predict_proba(self.X_train)[:, 1], "train data", "#273F4F")
        pltoneprc(self.y_test, self.model.predict_proba(self.X_test)[:, 1], "test data", "#A76545")
        pltoneprc(self.y, self.X[showTF], showName, "#8C8C8C")
        plt.title("PRC")
        plt.xlabel("Recall")
        plt.ylabel("Precision")
        plt.legend()
        # 设置右图为当前绘图区域
        ax1 = fig.add_subplot(gs[1, 0])
        plt.sca(ax1)
        pltoneroc(self.y_train, self.model.predict_proba(self.X_train)[:, 1], "train data", "#273F4F")
        pltoneroc(self.y_test, self.model.predict_proba(self.X_test)[:, 1], "test data", "#A76545")
        pltoneroc(self.y, self.X[showTF], showName, "#8C8C8C")
        plt.title("ROC")
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.legend()

        '''importance_df = pd.DataFrame({
            "feature": namelist,
            "coefficient":  model.coef_[0],
            "abs_importance": np.abs(model.coef_[0])
        }).sort_values(by="abs_importance", ascending=False)

        plt.figure(figsize=(6, 8))
        plt.barh(importance_df["feature"], importance_df["coefficient"])
        plt.xlabel("Coefficient")
        plt.title(f"Features importance (Logistic Regression)")
        plt.gca().invert_yaxis()
        plt.tight_layout()
        plt.show()'''
    
        # 获取逻辑回归模型的系数
        if self.predicttype == 'all':
            feature_importance = abs(self.model.coef_[0])
            map_vir = cm.get_cmap(name='Greens')  # 使用PiYG颜色映射
        else:
            feature_importance = self.model.coef_[0]
            map_vir = cm.get_cmap(name='PiYG')  # 使用PiYG颜色映射
        feature_importance = feature_importance / abs(feature_importance).max()
        # 按照重要性排序
        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + 0.5
        # 生成颜色映射
        ncolor = np.arange(0, len(self.X_test.columns), step=1)
        norm = plt.Normalize(ncolor.min(), ncolor.max())
        norm_color = norm(ncolor)
        colors = map_vir(norm_color)
        
        # 创建绘图
        ax2 = fig.add_subplot(gs[:, 1])
        plt.sca(ax2)
        ax2.barh(pos, feature_importance[sorted_idx], align='center', color=colors)
        if self.predicttype == 'all':
            ax2.set_xlim(0, 1)
        else:
            ax2.set_xlim(-1, 1)
        ax2.set_yticks(pos)
        ax2.set_yticklabels(np.array(self.X_test.columns)[sorted_idx], fontsize=8)
        ax2.set_xlabel("Relative Feature Importance")
        
        ax2.set_title("Feature Importance")
        plt.suptitle("Predict "+self.predicttype+" DEGs " +suptitle, fontsize=14, y=1.02)
        plt.tight_layout()
        plt.show()