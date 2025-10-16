import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, auc
import seaborn as sns
import warnings
import matplotlib.colors as mcolors
cols = list(mcolors.TABLEAU_COLORS.keys())


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

def scatter_with_rank(x_data, y_data,ranktype="diffrank",ifcolbar=True,**kwargs,):
    if ranktype == "diffrank":
        cor = makediffrank(x_data, y_data)
        plt.scatter(x_data,y_data,c=cor,cmap="coolwarm",**kwargs)
    elif ranktype == "sumrank":
        cor = makesumrank(x_data, y_data)
        plt.scatter(x_data,y_data,c=cor,cmap="Purples",**kwargs)
    elif ranktype == "sumrank0":
        cor = makesumrank_center0(x_data, y_data)
        plt.scatter(x_data,y_data,c=cor,cmap="RdBu_r",**kwargs)
        
    if ifcolbar:
        plt.colorbar()


def pltoneprc(stand,score,stand_label,color):
    precision, recall, _ = precision_recall_curve(stand,score)
    average_precision = average_precision_score(stand, score)
    label = stand_label+" "+str(round(average_precision,3))
    plt.step(recall, precision, where='post',color=color,label=label)
    
def pltoneroc(stand,score,stand_label,color):
    fpr, tpr, thresholds = roc_curve(stand, score)
    roc_auc = auc(fpr, tpr)
    if roc_auc>0.5:
        label = stand_label+" "+str(round(roc_auc,3))
        plt.plot(fpr, tpr, lw=2, label=label,color=color,alpha=0.8)
    else:
        label = stand_label+" "+str(round(1-roc_auc,3))
        plt.plot(1-fpr, 1-tpr, lw=2, label=label,color=color,alpha=0.8)


def showAUPRC(stand1,score1,label1='curve1',
            stand2=None,score2=None,label2=None,
            stand3=None,score3=None,label3=None,
            stand4=None,score4=None,label4=None,
            stand5=None,score5=None,label5=None,
            stand6=None,score6=None,label6=None,
            stand7=None,score7=None,label7=None,
            stand8=None,score8=None,label8=None,
            stand9=None,score9=None,label9=None,
            stand10=None,score10=None,label10=None,
            figsize=(4,4),legendloc="best",title='Precision-Recall curve',
            cols = list(mcolors.TABLEAU_COLORS.keys()),
            outname=None
            ):
    warnings.filterwarnings("ignore")

    plt.figure(figsize=figsize)

    pltoneprc(stand1,score1,label1,cols[1-1])
    if stand2 is not None and score2 is not None: pltoneprc(stand2,score2,label2,cols[2-1])
    if stand3 is not None and score3 is not None: pltoneprc(stand3,score3,label3,cols[3-1])
    if stand4 is not None and score4 is not None: pltoneprc(stand4,score4,label4,cols[4-1])
    if stand5 is not None and score5 is not None: pltoneprc(stand5,score5,label5,cols[5-1])
    if stand6 is not None and score6 is not None: pltoneprc(stand6,score6,label6,cols[6-1])
    if stand7 is not None and score7 is not None: pltoneprc(stand7,score7,label7,cols[7-1])
    if stand8 is not None and score8 is not None: pltoneprc(stand8,score8,label8,cols[8-1])
    if stand9 is not None and score9 is not None: pltoneprc(stand9,score9,label9,cols[9-1])
    if stand10 is not None and score10 is not None: pltoneprc(stand10,score10,label10,cols[10-1])
    
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(title)
    plt.legend(loc=legendloc)
    plt.tight_layout()
    if outname != None:
        plt.savefig(outname)
    #print(average_precision1,average_precision2)
    
def showROC(stand1,score1,label1='curve1',
            stand2=None,score2=None,label2=None,
            stand3=None,score3=None,label3=None,
            stand4=None,score4=None,label4=None,
            stand5=None,score5=None,label5=None,
            stand6=None,score6=None,label6=None,
            stand7=None,score7=None,label7=None,
            stand8=None,score8=None,label8=None,
            stand9=None,score9=None,label9=None,
            stand10=None,score10=None,label10=None,
            figsize=(4,4),legendloc="best",title='ROC Curve',outname=None):
    plt.figure(figsize=figsize)
    warnings.filterwarnings("ignore")

    pltoneroc(stand1,score1,label1,cols[1-1])
    if stand2 is not None and score2 is not None: pltoneroc(stand2,score2,label2,cols[2-1])
    if stand3 is not None and score3 is not None: pltoneroc(stand3,score3,label3,cols[3-1])
    if stand4 is not None and score4 is not None: pltoneroc(stand4,score4,label4,cols[4-1])
    if stand5 is not None and score5 is not None: pltoneroc(stand5,score5,label5,cols[5-1])
    if stand6 is not None and score6 is not None: pltoneroc(stand6,score6,label6,cols[6-1])
    if stand7 is not None and score7 is not None: pltoneroc(stand7,score7,label7,cols[7-1])
    if stand8 is not None and score8 is not None: pltoneroc(stand8,score8,label8,cols[8-1])
    if stand9 is not None and score9 is not None: pltoneroc(stand9,score9,label9,cols[9-1])
    if stand10 is not None and score10 is not None: pltoneroc(stand10,score10,label10,cols[10-1])

    '''
    fpr, tpr, thresholds = roc_curve(stand, abcscore)
    roc_auc = auc(fpr, tpr)
    label1 = stand1_label+" "+str(round(roc_auc,3))
    plt.plot(fpr, tpr, lw=2, label=label1 )
    

    if stand2 is not None and abcscore2 is not None:
        fpr2, tpr2, thresholds2 = roc_curve(stand2, abcscore2)
        roc_auc2 = auc(fpr2, tpr2)
        label2 = stand2_label+" "+str(round(roc_auc2,3))
        plt.plot(fpr2, tpr2, lw=2, label='NewMethod' )
    '''
        
    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(title)
    plt.legend(loc=legendloc)
    plt.tight_layout()
    #plt.show()
    if outname != None:
        plt.savefig(outname)


def corrABC(myabc_common,realabc_common,iflog=False,xlim=1,ylim=1):

    plt.figure(figsize=(4,4))
    # 绘制散点图
    if iflog:
        plt.scatter(np.log1p(myabc_common*100), np.log1p(realabc_common*100), 
                s=10, alpha=0.1, color='navy')
    else:
        plt.scatter(myabc_common, realabc_common, 
                s=10, alpha=0.1, color='navy')
    plt.xlim([0, xlim])
    plt.ylim([0, ylim])

    plt.xlabel("Built-in ABC (SCALE)")
    plt.ylabel("Original method (KR)")

    # 添加标题
    plt.title("Comparison of built-in ABC and original ABC")

    # 显示网格线
    plt.grid(True)

    # 显示图形
    plt.show()
    