
import sys
import seaborn as sns 
import warnings
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score

from .RPmodel import *
from .PRC_ROC import *


#def 
def makediffrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1-rank2)/len(list1)

def makesumrank(list1,list2):
    rank1 = list1.rank()
    rank2 = list2.rank()
    return (rank1+rank2)/(2*len(list1))


def plot_scatter_with_fit(x_data, y_data, x_label='X', y_label='Y',height=4, aspect=1,
                          insideplot=False,title='Scatter Plot with Linear Fit',**kwargs,):
    import warnings
    warnings.filterwarnings('ignore')
    # 创建数据框
    data = pd.DataFrame({x_label: x_data, y_label: y_data})
    
    # 使用 seaborn 绘制带有拟合直线的散点图
    if not insideplot:
        sns.lmplot(x=x_label, y=y_label, data=data, line_kws={'color': 'grey'},
                height=height, aspect=aspect,scatter_kws={**kwargs})
    else:
        sns.regplot(x=x_label, y=y_label, data=data, line_kws={'color': 'grey'},
                scatter_kws={**kwargs})
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title,fontsize=14)
    # 计算并显示相关性系数
    corr = np.corrcoef(x_data, y_data)[0, 1]
    plt.text(0.05, 0.95, f'Correlation: {corr:.3f}', transform=plt.gca().transAxes, 
             fontsize=14, verticalalignment='top')
    plt.tight_layout()

def lgfun(X,Y,printlm=False):
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=42)
    model = LogisticRegression()
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)

    if printlm:
        print(accuracy_score(y_test, y_pred))
        print("Coefficients:", model.coef_)
        print("Intercept:", model.intercept_)

    return(model.decision_function(X),model.coef_)

def importRumballResult(filepath):
    rnaseq = pd.read_csv(filepath,sep="\t")
    rnaseq = rnaseq[rnaseq["type"]=="protein_coding"]
    rnaseq_bed = rnaseq[['chromosome',"start",'end',"Ensembl ID","genename","strand","logFC","FDR","logCPM"]]
    rnaseq_bed.to_csv("allgene_withfdr.bed",sep="\t",header=None, index=None)

def obtainRP(genebedfile,gtfile,bam_ctrl,bam_treat,species):
    rpobj_ctrl = myRP(genebedfile,bam_ctrl,gtfile,padding = int(1e5),
             decay_distance=10000,rpweight="classic",species=species)
    rpdf_ctrl = rpobj_ctrl.calculateRP_wang_allgene(16)

    rpobj_treat = myRP(genebedfile,bam_treat,gtfile,padding = int(1e5),
             decay_distance=10000,rpweight="classic",species=species)
    rpdf_treat = rpobj_treat.calculateRP_wang_allgene(16)

    rpdf_treat.to_csv("rpdf_treat.tsv",sep="\t",
                  header=None,index=None)
    rpdf_ctrl.to_csv("rpdf_ctrl.tsv",sep="\t",
                    header=None,index=None)

def predictDEG(ctrlcsv,treatcsv,logRg=True,fdrcutoff=0.01,
               fccol=6,fdrcol=7,cpmcol= 8,rpcol = 9,titleend=None,
               plotScatter=False,plotAllDeg=False,plotROC=False,returnlist=False):
    
    rpdf_ctrl = pd.read_csv(ctrlcsv,sep="\t",header=None)
    rpdf_treat = pd.read_csv(treatcsv,sep="\t",header=None)

    if logRg:
        print("log the Rg Score")
        rpdifflist = np.log1p(rpdf_treat[rpcol])- np.log1p(rpdf_ctrl[rpcol])
        rpctrllist = np.log1p(rpdf_ctrl[rpcol])
        rptreatlist = np.log1p(rpdf_treat[rpcol])
    else:
        rpdifflist = rpdf_treat[rpcol]- rpdf_ctrl[rpcol]
        rpctrllist = rpdf_ctrl[rpcol]
        rptreatlist = rpdf_treat[rpcol]

    alldeglist = (abs(rpdf_ctrl[fccol])>1) & (rpdf_ctrl[fdrcol]<fdrcutoff)
    updeglist = (rpdf_ctrl[fccol]>1) & (rpdf_ctrl[fdrcol]<fdrcutoff)
    downdeglist = (rpdf_ctrl[fccol]<-1) & (rpdf_ctrl[fdrcol]<fdrcutoff)

    logfclist = rpdf_ctrl[fccol]
    fdrlist = rpdf_ctrl[fdrcol]
    cpmlist = rpdf_ctrl[cpmcol]

    if plotScatter:
        print("correlation between Gene-logFC and Rg-diff")
        plot_scatter_with_fit(logfclist, rpdifflist, 
                        x_label='Gene TPM changes(logFC)', 
                        y_label='Gene Rg changes(logFC)', title=None)
        
        print("correlation between Gene-TPM and Rg")
        plot_scatter_with_fit(cpmlist, (rpctrllist+rptreatlist)/2, 
                        x_label='Gene CPM (sample averaged)', 
                        y_label='Gene Rg (sample averaged)', title=None)
    
    upgenebool = (pd.Series(fdrlist)<0.01) & (pd.Series(logfclist) >1)
    downgenebool = (pd.Series(fdrlist)<0.01) & (pd.Series(logfclist) <-1)
    nondegbool = (pd.Series(fdrlist) >0.01) & (abs(pd.Series(logfclist)) <1)
    alldegbool = (pd.Series(fdrlist) < 0.01) & (abs(pd.Series(logfclist)) >1)

    X1 = np.column_stack((rpctrllist, rptreatlist, rpdifflist))  # 合并特征为二维数组
    X2 = np.column_stack((rpdifflist,))  # 合并特征为二维数组
    X3 = np.column_stack((rpctrllist,))  # 合并特征为二维数组
    X4 = np.column_stack((rptreatlist,))  # 合并特征为二维数组

    if plotAllDeg:
        plotdeglist = [alldeglist,updeglist,downdeglist]
    else:
        plotdeglist = [updeglist,downdeglist]

    for Y in plotdeglist:
        lg1,coef1 = lgfun(X1,Y)
        lg2,_ = lgfun(X2,Y)
        lg3,_ = lgfun(X3,Y)
        lg4,_ = lgfun(X4,Y)

        text_str = (f"ctrl_coef: {coef1[0][0]:.4f}\n"
            f"treat_coef: {coef1[0][1]:.4f}\n"
            f"diff_coef: {coef1[0][2]:.4f}")
        if np.array_equal(Y,alldeglist):
            title="Predict All DEGs"+" " +titleend
        elif np.array_equal(Y,updeglist):
            title="Predict Up DEGs"+" " +titleend
        elif np.array_equal(Y,downdeglist):
            title="Predict Down DEGs"+" " +titleend
        
        if plotROC:
            showROC(Y,lg3,'OnlyCtrl',
                Y,lg4,'OnlyTreat',
                Y,lg2,'OnlyDiff',
                Y,lg1,'Combined',title=title)
            plt.text(0.95, 0.5, text_str, fontsize=9, ha='right', va='center')
        
        showAUPRC(Y,lg3,'OnlyCtrl',
            Y,lg4,'OnlyTreat',
            Y,lg2,'OnlyDiff',
            Y,lg1,'Combined',title=title)
        plt.text(0.95, 0.5, text_str, fontsize=9, ha='right', va='center')

    if returnlist:
        return rpdifflist,updeglist,downdeglist
    
    