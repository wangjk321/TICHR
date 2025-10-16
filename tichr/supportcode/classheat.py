import os
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
import matplotlib.cm
import matplotlib.colors as mcolors
import warnings
from scipy import stats
warnings.filterwarnings('ignore')

def plotheat(df,type,width=4,height=7,ncluster=3,clustermethod='kmeans',
             iflog=True,normType="zscore",vmax=None,vmin=None,outmatrixname="out_matrix",
             labellist=None,col_cluster=True,ifpca=False,num_pca=10):

    df = normalizedf(df,iflog=iflog,ntype=normType)

    if type == "raw":
        sns.clustermap(df,cmap="RdBu_r",center=0,yticklabels=False,xticklabels=True,vmax=vmax,vmin=vmin,
                       figsize=(width,height),row_cluster=True)
    elif type == "sort":
        sortDF = df.sort_values(by=df.columns[0])
        sns.clustermap(sortDF, cmap="RdBu_r", center=0,vmax=vmax,vmin=vmin,
               row_cluster=False, yticklabels=False, figsize=(width,height), xticklabels=True)
    
    elif type == "cluster":
        clusterDF = df.copy(deep=True)
        clusterDF["clusters"] = DimReduction(clusterDF, ncluster,ifpca=ifpca,num_pca=num_pca,seed=0,clustermethod=clustermethod)
        clusterDF = clusterDF.sort_values(by=['clusters'])

        clusterDF.to_csv(outmatrixname+".tsv",sep="\t")

        # prepare for plot
        colist = list(mcolors.TABLEAU_COLORS.keys())
        lut = dict(zip(clusterDF["clusters"].unique(), colist))
        row_colors = clusterDF["clusters"].map(lut)

        clusterDF_nolabel = clusterDF.iloc[:,0:clusterDF.shape[1]-1]
        sns.clustermap(clusterDF_nolabel, row_colors=row_colors, row_cluster=False,
                         cmap="RdBu_r", yticklabels=False, center=0,vmax=vmax,vmin=vmin,
                         figsize=(width,height), xticklabels=True)
    elif type == "raw+label":
        colist = list(mcolors.TABLEAU_COLORS.keys())
        lut = dict(zip(labellist.unique(), colist))
        row_colors = labellist.map(lut)

        sns.clustermap(df, row_colors=row_colors, row_cluster=False,
                         cmap="RdBu_r", yticklabels=False, center=0,vmax=vmax,vmin=vmin,
                         figsize=(width,height), xticklabels=True)


def normalizedf(df,iflog=True,ntype="zscore"):
    if iflog:
        df = np.log1p(df)
    if ntype == "zscore":
        normdf = (df - df.mean())/df.std()
    elif ntype == "scale0to1":
        normdf = df/(df.max().max())
    elif ntype == "MAD":
        normdf = (df - df.median())/stats.median_abs_deviation(df)
    elif ntype == "no":
        normdf = df
    else:
        print("please choose the correct type")
        exit(1)
    return(normdf)

def DimReduction(data, ncluster,ifpca=False,num_pca=10,seed=0,clustermethod='kmeans'):
    if ifpca:
        pca = PCA(n_components=num_pca)
        pca.fit(data)
        matrix = pca.transform(data)
    else:
        matrix = data

    if clustermethod=='minikmeans':
        model = MiniBatchKMeans(random_state=seed, n_clusters=ncluster, max_iter=10000, batch_size=100)
    elif clustermethod=='kmeans':
        from sklearn.cluster import KMeans
        model = KMeans(random_state=seed, n_clusters=ncluster)
    elif clustermethod=='spectral':
        from sklearn.cluster import SpectralClustering
        model = SpectralClustering(random_state=seed, n_clusters=ncluster)
    elif clustermethod=='meanshift':
        from sklearn.cluster import MeanShift
        model = MeanShift()
    elif clustermethod=='dbscan':
        from sklearn.cluster import DBSCAN
        model = DBSCAN(eps=3, min_samples=2)
    elif clustermethod=='affinity':
        from sklearn.cluster import AffinityPropagation
        AffinityPropagation(random_state=seed)
    else:
        raise("please use the correst cluster method")
        exit(1)
    outlabels = model.fit_predict(matrix)
    return outlabels