from operator import gt
import pandas as pd
import seaborn as sns
import numpy as np
from scipy.stats import pearsonr,spearmanr
from sklearn.linear_model import LinearRegression
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm  
import xgboost as xgb
from functools import partial
from .tichr import *
from .preprocess_hic import *


class largescale:
    def __init__(self,epiFile,genefile,hicFile,gt,hicRes=25000):
        self.dhs_df = pd.read_csv(epiFile,sep="\t")
        self.gene_expr_df = pd.read_csv(genefile,sep="\t")
        self.hicRes = hicRes
        self.nomhicdf = gethicfile(hicFile,hicRes,'rawhic_sparse',list(self.gene_expr_df["chr"].unique()),
                                hicnorm="VC_SQRT",gt=gt,threads=12,
                                further_normalize_type="no_further")
    
    def calculate_importance(self,gene_row,method="Expon",maxdistance=200000,threads=8,halfDistance=10000):
        dhs_df=self.dhs_df
        gene_chr = gene_row["chr"]
        gene_tss = int(gene_row["tss"])
        gene_symbol = gene_row["geneSymbol"]
        gene_values = gene_row[3:].to_numpy(dtype=float)
        
        region_start = max(0, gene_tss - maxdistance)
        region_end = gene_tss + maxdistance
        
        # 查找上下游 200kb 范围内的 DHS 位点
        dhs_in_range = dhs_df[(dhs_df["chr"] == gene_chr) & 
                            (dhs_df["start"] <= region_end) & 
                            (dhs_df["end"] >= region_start)]
        numpeak = len(dhs_in_range)
        if numpeak <1:
            return [], []
            
        if method in ["Pearson","Spearman"]:
            best_correlation = None  # 初始化最佳相关性
            best_dhs_values = None  # 初始化最佳 DHS 行

            result = []
            for _, dhs_row in dhs_in_range.iterrows():
                dhs_values = dhs_row[3:].to_numpy(dtype=float)
                # 检查数组是否为常数
                if len(set(gene_values)) == 1 or len(set(dhs_values)) == 1:
                    correlation = 0 
                else:
                    if method == "Pearson":
                        correlation = pearsonr(dhs_values, gene_values)[0]
                    elif method == "Spearman":
                        correlation = spearmanr(dhs_values, gene_values)[0]
                
                # 更新最佳相关性和最佳 DHS 行
                if best_correlation is None or abs(correlation) > abs(best_correlation):
                    best_correlation = correlation
                    best_dhs_values = dhs_values
                
                if not best_correlation:
                    correlation = 0
                    best_dhs_values = 0 * dhs_values
                
                result.append({
                    "dhs_chr": dhs_row["chr"],
                    "dhs_start": dhs_row["start"],
                    "dhs_end": dhs_row["end"],
                    "gene_chr": gene_chr,
                    "gene_tss": gene_tss,
                    "gene_symbol": gene_symbol,
                    "importance": correlation
                })
            prediction = best_dhs_values
                
        elif method in ["OLS","XGB","Expon","Expon+OLS","Expon+XGB","HiC","HiC+Linear","HiC+XGB"]:
            if method in ["XGB","Expon+XGB","HiC+XGB"]:
                X = dhs_in_range.iloc[:, 3:].T  # DHS信号
                y = gene_row[3:].values  # 基因表达值
            elif method in ["Expon","OLS","Expon+OLS","HiC","HiC+OLS"]:
                
                X = dhs_in_range.iloc[:, 3:].values.T  # DHS信号
                y = gene_row[3:].values  # 基因表达值
            
            if method == "OLS":
                # 初始化线性回归模型
                model = LinearRegression()
                # 拟合模型
                model.fit(X, y)
                # 提取特征系数，表示DHS位点的重要性
                importance = model.coef_
                prediction = model.predict(X)
                
            elif method == "XGB":
                # 将数据转换为 DMatrix 格式
                dtrain = xgb.DMatrix(X, label=y)
                # 设置 XGBoost 回归模型参数
                params = {
                    'objective': 'reg:squarederror',  # 使用均方误差作为目标函数
                    'max_depth': 6,
                    'eta': 0.3,
                    'eval_metric': 'rmse',  # 根均方误差作为评估指标
                    'nthread': threads
                }
                # 训练模型
                num_boost_round = 100  # 您可以根据需要调整迭代轮数
                bst = xgb.train(params, dtrain, num_boost_round)
                # 提取特征重要性
                importance_dict = bst.get_score(importance_type='gain')
                ordered_importance_dict = {}
                for col in X.columns:
                    ordered_importance_dict[str(col)] = importance_dict.get(str(col), 0)  # 如果不在字典中，默认值为0
                importance = list(ordered_importance_dict.values())
                
                prediction = bst.predict(dtrain)
                
            elif method == "Expon":
                z = abs(gene_tss - (dhs_in_range["start"]+dhs_in_range["end"])/2 )
                weight = np.array(2 ** -(z/halfDistance))
                rgx = X * weight
                importance = rgx.sum(axis=0)
                prediction = rgx.sum(axis=1)
            
            elif method == "HiC":
                weight = []
                for peaki in range(numpeak):
                    peakPos = np.array((dhs_in_range["start"]+dhs_in_range["end"])/2)
                    tssPos = gene_tss
                    weight_peaki = makeWeightFunction("hic",peakPos[peaki],tssPos,hicProcessedData=self.nomhicdf,
                                    hicRes=self.hicRes,geneChr=gene_chr,hicProcessedDataType="rawhic_sparse",ifUseHiCRef=False)
                    weight.append(weight_peaki)
                    
                rgx = X * np.array(weight)
                importance = rgx.sum(axis=0)
                prediction = rgx.sum(axis=1)
                
            elif method == "HiC+OLS":
                weight = []
                for peaki in range(numpeak):
                    peakPos = np.array((dhs_in_range["start"]+dhs_in_range["end"])/2)
                    tssPos = gene_tss
                    weight_peaki = makeWeightFunction("hic",peakPos[peaki],tssPos,hicProcessedData=self.nomhicdf,
                                   hicRes=self.hicRes,geneChr=gene_chr,hicProcessedDataType="rawhic_sparse",ifUseHiCRef=False)
                    weight.append(weight_peaki)
                    
                X = X * np.array(weight)
                model = LinearRegression()
                # 拟合模型
                model.fit(X, y)
                # 提取特征系数，表示DHS位点的重要性
                importance = model.coef_
                prediction = model.predict(X)
            
            elif method == "HiC+XGB":
                weight = []
                for peaki in range(numpeak):
                    peakPos = np.array((dhs_in_range["start"]+dhs_in_range["end"])/2)
                    tssPos = gene_tss
                    weight_peaki = makeWeightFunction("hic",peakPos[peaki],tssPos,hicProcessedData=self.nomhicdf,
                    hicRes=self.hicRes,geneChr=gene_chr,hicProcessedDataType="rawhic_sparse",ifUseHiCRef=False)
                    weight.append(weight_peaki)
                    
                X = X * np.array(weight)
                
                # 将数据转换为 DMatrix 格式
                dtrain = xgb.DMatrix(X, label=y)
                # 设置 XGBoost 回归模型参数
                params = {
                    'objective': 'reg:squarederror',  # 使用均方误差作为目标函数
                    'max_depth': 6,
                    'eta': 0.3,
                    'eval_metric': 'rmse',  # 根均方误差作为评估指标
                    'nthread': threads
                }
                # 训练模型
                num_boost_round = 100  # 您可以根据需要调整迭代轮数
                bst = xgb.train(params, dtrain, num_boost_round)
                
                # 提取特征重要性
                importance_dict = bst.get_score(importance_type='gain')
                ordered_importance_dict = {}
                for col in X.columns:
                    ordered_importance_dict[str(col)] = importance_dict.get(str(col), 0)  # 如果不在字典中，默认值为0
                importance = list(ordered_importance_dict.values())
                
                prediction = bst.predict(dtrain)
                
            elif method == "Expon+OLS":
                z = abs(gene_tss  - (dhs_in_range["start"]+dhs_in_range["end"])/2 )
                weight = np.array(2 ** -(z/halfDistance))
                X = X * weight
                model = LinearRegression()
                # 拟合模型
                model.fit(X, y)
                # 提取特征系数，表示DHS位点的重要性
                importance = model.coef_
                prediction = model.predict(X)
                
            elif method == "Expon+XGB":
                z = abs(gene_tss  - (dhs_in_range["start"]+dhs_in_range["end"])/2 )
                weight = np.array(2 ** -(z/halfDistance))
                X = X * weight
                # 将数据转换为 DMatrix 格式
                dtrain = xgb.DMatrix(X, label=y)
                # 设置 XGBoost 回归模型参数
                params = {
                    'objective': 'reg:squarederror',  # 使用均方误差作为目标函数
                    'max_depth': 6,
                    'eta': 0.3,
                    'eval_metric': 'rmse',  # 根均方误差作为评估指标
                    'nthread': threads
                }
                # 训练模型
                num_boost_round = 100  # 您可以根据需要调整迭代轮数
                bst = xgb.train(params, dtrain, num_boost_round)
                
                # 提取特征重要性
                importance_dict = bst.get_score(importance_type='gain')
                ordered_importance_dict = {}
                for col in X.columns:
                    ordered_importance_dict[str(col)] = importance_dict.get(str(col), 0)  # 如果不在字典中，默认值为0
                importance = list(ordered_importance_dict.values())
                
                prediction = bst.predict(dtrain)
            
            result = []
            for i, (_, dhs_row) in enumerate(dhs_in_range.iterrows()):
                result.append({
                    "dhs_chr": dhs_row["chr"],
                    "dhs_start": dhs_row["start"],
                    "dhs_end": dhs_row["end"],
                    "gene_chr": gene_chr,
                    "gene_tss": gene_tss,
                    "gene_symbol": gene_symbol,
                    "importance": importance[i]  # DHS 位点的重要性
                })

        prediction_list = [gene_chr, str(gene_tss),gene_symbol]
        prediction_list.extend(list(prediction))
        
        return(result,[prediction_list])
    
    def testfun(self,):
        self.calculate_importance(self.gene_expr_df.loc[0],method="Expon")

    def process(self, method="Expon",maxdistance=200000,threads=8,halfDistance=10000,
                outname="outname"):
        if method in ["Expon+Linear","HiC","HiC+Linear","Pearson","Spearman","Linear","Expon"]:
            results = []
            prediction_results = [] 
            with ProcessPoolExecutor(threads) as executor:
                calculate_with_method = partial(self.calculate_importance, method=method,maxdistance=maxdistance,
                                                threads=threads,halfDistance=halfDistance)
                futures = {executor.submit(calculate_with_method, row): row for _, row in self.gene_expr_df.iterrows()}

                # 使用 tqdm 包装 as_completed 以显示进度条
                for future in tqdm(as_completed(futures), total=len(futures), desc="Processing genes"):
                    results.extend(future.result()[0])
                    prediction_results.extend(future.result()[1])
                    
                    # 移除已完成的 Future
                    del futures[future]

            # 转换结果为 DataFrame
            importance_df = pd.DataFrame(results)
            prediction_df = pd.DataFrame(prediction_results)
            prediction_df.columns = self.gene_expr_df.columns

        elif method in ["HiC+XGB","Expon+XGB","XGB"]:
            results = []
            prediction_results = []
            for _, gene_row in tqdm(self.gene_expr_df.iterrows(), total=len(self.gene_expr_df), desc="Processing genes"):
                output = self.calculate_importance(gene_row, method=method,maxdistance=maxdistance,
                                                threads=threads,halfDistance=halfDistance)
                results.extend(output[0])
                prediction_results.extend(output[1])
            
            importance_df = pd.DataFrame(results)
            prediction_df = pd.DataFrame(prediction_results)
            prediction_df.columns = self.gene_expr_df.columns
        
        importance_df.to_csv(outname+"_"+method+"_importance.tsv", index=False,sep="\t")
        prediction_df.to_csv(outname+"_"+method+"_predicted.tsv", index=False,sep="\t")
        
        #return(importance_df,prediction_df)



import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import umap
from sklearn.preprocessing import StandardScaler
import warnings
warnings.filterwarnings("ignore")
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score

def dimensionality_reduction_plot(df, method='UMAP',label="Gene expression",
                                  seed=40,usePCAkmean=True,legend=False):
    # 提取基因表达值部分 (从第5列开始)
    data = df.iloc[:, 5:]

    # 数据标准化
    scaler = StandardScaler()
    data_scaled = scaler.fit_transform(data.T)
    
    # 根据选择的降维方法应用不同的模型
    if method == 'PCA':
        model = PCA(n_components=2,random_state=42)
    elif method == 'UMAP':
        model = umap.UMAP(n_components=2, random_state=42)
    elif method == 'TSNE':
        model = TSNE(n_components=2, random_state=42)
    else:
        raise ValueError("Invalid method. Choose from 'PCA', 'UMAP', or 'TSNE'.")
    
    # 执行降维
    embedding = model.fit_transform(data_scaled)    
    # 将 cellDF 的 CANONICAL 列与数据列名匹配，提取对应的 TISSUE 列
    tissue_colors = cellDF.loc[data.columns, 'TISSUE'] 
    
    if not usePCAkmean:
        kmeans = KMeans(n_clusters=15, random_state=seed)
        predicted_labels = kmeans.fit_predict(embedding)
    else:
        model2 = PCA(n_components=30,random_state=seed)
        embedding2 = model2.fit_transform(data_scaled)
        
        kmeans = KMeans(n_clusters=15, random_state=seed)
        predicted_labels = kmeans.fit_predict(embedding2)
    
    # 计算 ARI
    if cellDF is not None:
        ari = adjusted_rand_score(tissue_colors, predicted_labels)
        print("ARI Score:", ari)
    
    
    # 使用 Seaborn 的调色板生成颜色映射
    unique_tissues = tissue_colors.unique()
    palette = sns.color_palette("tab20", len(unique_tissues))
    color_map = dict(zip(unique_tissues, palette))
    colors = tissue_colors.map(color_map)

    # 绘图
    plt.figure(figsize=(3, 3.2))
    plt.scatter(embedding[:, 0], embedding[:, 1], c=colors, alpha=0.7,)
    #plt.text(0.5,0.9,f"ARI Score: {ari:.3f}",transform=plt.gca().transAxes,fontsize=12)
    
    # 图例
    handles = [plt.Line2D([0], [0], marker='o', color=color, linestyle='', markersize=8) for color in color_map.values()]
    
    if legend:
        plt.legend(handles, color_map.keys(), title="Tissue Type", loc='best', fontsize=9,
               bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    
    # 添加基因标签（例如 geneSymbol列中的名称）
    #for i, label in enumerate(data.columns):
    #    plt.annotate(label, (embedding[i, 0], embedding[i, 1]), fontsize=10, alpha=0.7)

    plt.xlabel(f'{method} 1',fontsize=10)
    plt.ylabel(f'{method} 2',fontsize=10)
    plt.title(f'{label}',fontsize=14)
    #plt.legend()
    plt.tight_layout()
    plt.savefig("/home/wang/Tichr/2025May/pdf/fig6/RED-umap/"+label+".pdf")
    return(ari)