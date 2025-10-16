import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def plot_one_gene(file_path, timepoint_counts, title, time_points=None, outdir=None):

    # 固定时间点顺序
    if time_points==None:
        time_points=[]
        for i in range(timepoint_counts):
            time_points.append(f"T{str(i+1)}")

    names = ['peakChr', 'peakStart', 'peakEnd', 'epigenomeActivity',
            'geneID', 'geneChr', 'geneStart', 'geneEnd', 'geneStrand',
            'geneSymbol'] + time_points + ['mark']
    df = pd.read_csv(file_path, sep='\t', names=names)

    columns = time_points

    # 创建 figure + 双 y 轴
    fig, ax1 = plt.subplots(figsize=(3, 3))
    #ax2 = ax1.twinx()
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    all_logs = []
    total_vals = np.zeros(len(time_points))

    # 第一步：画个体曲线到 ax1
    for _, row in df.iterrows():
        vals = [row[c] + 1 for c in columns]
        log_vals = np.log10(vals)
        all_logs.extend(log_vals)

        if row['mark'] == 'Target':
            ax1.plot(time_points, log_vals,
                     color='#7E9BB7', alpha=0.95,marker='^',markersize=3,
                     linewidth=5, zorder=3,label="selected")
        elif row['mark'] == 'N':
            ax1.plot(time_points, log_vals,
                     color='#979998', alpha=0.2,
                     linewidth=1, zorder=1)
        else:
            pass

        # 同时累加原始值用于 sum
        total_vals += [row[c] for c in columns]
    

    # 关闭网格
    ax1.grid(False)

    # x 轴紧贴边界
    ax1.set_xlim(time_points[0], time_points[-1])
    y_min, y_max = min(all_logs), max(all_logs)
    ax1.set_ylim(y_min, y_max)

    # 左侧 y 轴——个体曲线
    ax1.set_ylabel(f'Rgx (log1p)', fontsize=12)
    # 匹配 Target 曲线颜色
    ax1.yaxis.label.set_color('#01579B')
    # ax1.tick_params(axis='y', colors='#C69287')

    ax1.set_xlabel('Time', fontsize=12)
    plt.legend()
    plt.title(title + "", fontsize=12)

    fig.tight_layout()
    if outdir: fig.savefig(f"{outdir}/{title.replace(' ', '_')}.pdf",)
    plt.show()

import pandas as pd
import numpy as np

def lfc_two_timepoints(rgx1, rgx2, extra_df=None, extra_name=None):
    rgx1df = pd.read_csv(rgx1, sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])
    rgx2df = pd.read_csv(rgx2, sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])

    # 1) basedf 去重 + 重置索引
    basedf = (
        rgx1df[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol']]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    rgx1_sum = (
        rgx1df.groupby('geneID', as_index=False)
        .agg({'epigenomeActivity': 'sum', 'weight': 'sum', 'Rgx_rawvalue': 'sum'})
        .rename(columns={'epigenomeActivity': 'epi1', 'weight': 'hic1', 'Rgx_rawvalue': 'rgx1'})
    )
    rgx2_sum = (
        rgx2df.groupby('geneID', as_index=False)
        .agg({'epigenomeActivity': 'sum', 'weight': 'sum', 'Rgx_rawvalue': 'sum'})
        .rename(columns={'epigenomeActivity': 'epi2', 'weight': 'hic2', 'Rgx_rawvalue': 'rgx2'})
    )

    basedf = basedf.merge(rgx1_sum, on='geneID', how='left')
    basedf = basedf.merge(rgx2_sum, on='geneID', how='left')

    with np.errstate(divide='ignore', invalid='ignore'):
        basedf['EnhancerLFC'] = np.where(
            basedf['epi1'].eq(0), np.nan, np.log2(basedf['epi2'] / basedf['epi1'])
        )
        basedf['ContactLFC'] = np.where(
            basedf['hic1'].eq(0), np.nan, np.log2(basedf['hic2'] / basedf['hic1'])
        )
        basedf['RgxLFC'] = np.where(
            basedf['rgx1'].eq(0), np.nan, np.log2(basedf['rgx2'] / basedf['rgx1'])
        )

    if extra_df: 
        extra_df = pd.read_csv(extra_df, sep="\t", header=None, names=[
            'geneChr', 'geneStart', 'geneEnd', 'geneSymbol', 'geneID', 'geneStrand', f'{extra_name}1', f'{extra_name}2'
        ])
        extra_df[f'{extra_name}LFC'] = np.where(
            extra_df[f'{extra_name}1'].eq(0), np.nan, np.log2(extra_df[f'{extra_name}2'] / extra_df[f'{extra_name}1'])
        )
        extra_df = extra_df.drop_duplicates(subset=['geneSymbol'])
        basedf = basedf.merge(extra_df[['geneSymbol', f'{extra_name}LFC']], on='geneSymbol', how='left')
        basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC", f"{extra_name}LFC"]]
    else:
        basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC"]]
    return basedf

def change_or_not(lfcdf, names, change_marks, thresholds=None, plot_it=False, title="LFC Boxplot", ylim=([-2,2]), outdir=None):
    if thresholds is None:
        thresholds = [0.1] * len(names)

    for name, mark, th in zip(names, change_marks, thresholds):
        if mark:
            lfcdf = lfcdf[lfcdf[name + 'LFC'].abs() > th]
        else:
            lfcdf = lfcdf[lfcdf[name + 'LFC'].abs() <= th]

    if plot_it:
        features = names
        data = [lfcdf[f + 'LFC'].dropna() for f in features]

        fig, ax = plt.subplots(figsize=(4, 3.75))

        box = ax.boxplot(
            data,
            patch_artist=True,
            showfliers=False,
            medianprops={'linestyle': '--', 'color': 'black'},
            widths=0.7
        )

        colors = ['#E07F86', '#8FA2CD', '#F8BC7E', '#6BC179']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)

        ax.set_xticks(range(1, len(features) + 1))
        ax.set_xticklabels(features, rotation=0, ha='center', fontsize=11)
        ax.set_title("", fontsize=12)
        ax.set_ylim(ylim)
        ax.set_ylabel('Log₂ Fold Change', fontsize=13)

        plt.tight_layout()
        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(f"{outdir}/ChangeOrNot.png", dpi=300)
        plt.show()
    return lfcdf

def two_conflict(lfcdf, names, change_marks, thresholds=None, plot_it=False, title="LFC Boxplot", ylim=([-2,2]), outdir=None):
    if thresholds is None:
        thresholds = [0.1] * len(names)

    for name, mark, th in zip(names, change_marks, thresholds):
        if mark:
            lfcdf = lfcdf[lfcdf[name + 'LFC'].abs() > th]
        else:
            lfcdf = lfcdf[lfcdf[name + 'LFC'].abs() <= th]

    if plot_it:
        features = names
        data = [lfcdf[f + 'LFC'].dropna() for f in features]

        fig, ax = plt.subplots(figsize=(4, 3.75))

        box = ax.boxplot(
            data,
            patch_artist=True,
            showfliers=False,
            medianprops={'linestyle': '--', 'color': 'black'},
            widths=0.7
        )

        colors = ['#E07F86', '#8FA2CD', '#F8BC7E', '#6BC179']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)

        ax.set_xticks(range(1, len(features) + 1))
        ax.set_xticklabels(features, rotation=0, ha='center', fontsize=11)
        ax.set_title("", fontsize=12)
        ax.set_ylim(ylim)
        ax.set_ylabel('Log₂ Fold Change', fontsize=13)

        plt.tight_layout()
        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(f"{outdir}/ChangeOrNot.png", dpi=300)
        plt.show()
    return lfcdf

def two_conflict(lfcdf, names, change_marks, thresholds=None, plot_it=False, title="LFC Boxplot", ylim=([-2,2]), outdir=None):
    if thresholds is None:
        thresholds = [0.1] * len(names)

    for name, mark, th in zip(names, change_marks, thresholds):
        if mark == "Up":
            lfcdf = lfcdf[lfcdf[name + 'LFC'] > th]
        elif mark == "Down":
            lfcdf = lfcdf[lfcdf[name + 'LFC'] < -th]
        elif mark == "NoChange":
            lfcdf = lfcdf[(lfcdf[name + 'LFC'] >= -th) & (lfcdf[name + 'LFC'] <= th)]
        else:
            print("Error: change_marks should be 'Up', 'Down' or 'NoChange'")
            exit(1)

    if plot_it:
        features = names
        data = [lfcdf[f + 'LFC'].dropna() for f in features]

        fig, ax = plt.subplots(figsize=(4, 3.75))

        box = ax.boxplot(
            data,
            patch_artist=True,
            showfliers=False,
            medianprops={'linestyle': '--', 'color': 'black'},
            widths=0.7
        )

        colors = ['#E07F86', '#8FA2CD', '#F8BC7E', '#6BC179']
        for patch, color in zip(box['boxes'], colors):
            patch.set_facecolor(color)

        ax.set_xticks(range(1, len(features) + 1))
        ax.set_xticklabels(features, rotation=0, ha='center', fontsize=11)
        ax.set_title("", fontsize=12)
        ax.set_ylim(ylim)
        ax.set_ylabel('Log₂ Fold Change', fontsize=13)

        plt.tight_layout()
        if outdir:
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            plt.savefig(f"{outdir}/Conflict.png", dpi=300)
        plt.show()
    return lfcdf

import pandas as pd

def calculate_lfc(df, features, timepoint_counts):
    """给lfc_multi_timepoints用的"""
    lfc_cols = []  
    lfc_cols_list = []

    for feat in features:
        tmp_lfc_cols = []
        t1_col = f"{feat}T1"
        if t1_col not in df.columns:
            continue

        lfc_t1 = f"{feat}_LFC_T1"
        df[lfc_t1] = 0.0
        lfc_cols.append(lfc_t1)
        tmp_lfc_cols.append(lfc_t1)

        for t in range(2, timepoint_counts + 1):
            t_col = f"{feat}T{t}"
            if t_col not in df.columns:
                continue

            new_col = f"{feat}_LFC_T{t}"
            df[new_col] = np.where(
                df[t1_col] != 0,
                np.log2(df[t_col] / df[t1_col]),
                np.nan
            )
            lfc_cols.append(new_col)
            tmp_lfc_cols.append(new_col)
        lfc_cols_list.append(tmp_lfc_cols)

    final_cols = list(df.columns[:6]) + lfc_cols
    df_final = df[final_cols]

    return df_final, lfc_cols_list

def lfc_multi_timepoints(RgxDfs, extra_file=None, extra_name=None, time_points=None):
    #建一个基础的df，后面对它左连接
    rgx1df = pd.read_csv(RgxDfs[0], sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])
    basedf = (
        rgx1df[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol']]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    if time_points==None:
        time_points=[]
        for i in range(len(RgxDfs)):
            time_points.append(f"T{str(i+1)}")

    for rgxfile, time in zip(RgxDfs, time_points):
        rgxdf = pd.read_csv(rgxfile, sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])
        rgx_sum = (
            rgxdf.groupby('geneID', as_index=False)
            .agg({'epigenomeActivity': 'sum', 'weight': 'sum', 'Rgx_rawvalue': 'sum'})
            .rename(columns={'epigenomeActivity': f'Enhancer{time}', 'weight': f'Contact{time}', 'Rgx_rawvalue': f'Rgx{time}'})
        )
        basedf = basedf.merge(rgx_sum, on='geneID', how='left')

    if extra_file: 
        extra_time_points=[]
        for i in range(len(RgxDfs)):
            extra_time_points.append(f"{extra_name}T{str(i+1)}")
        names = ["geneChr", "geneStart", "geneEnd", "geneSymbol", "geneID", "geneStrand"] + extra_time_points
        extra_df = pd.read_csv(extra_file, sep='\t', header=None, names=names)
        extra_df = extra_df.drop_duplicates(subset=['geneSymbol'])
        basedf = basedf.merge(extra_df[['geneSymbol'] + extra_time_points], on='geneSymbol', how='left')
        # basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC", f"{extra_name}LFC"]]
    # else:
    # basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC"]]

    resdf, lfc_cols_list = calculate_lfc(basedf, ["Enhancer", "Contact", "Rgx"] + ([extra_name] if extra_file else []), timepoint_counts=len(RgxDfs))
    return resdf, lfc_cols_list


def calculate_sts_series(df, timepoint_cols):
    """给compute_sts用的"""
    sts = []
    baseline_vals = df[timepoint_cols[0]].values.astype(float)
    y_vals = df[timepoint_cols].values.astype(float)

    for baseline, y in zip(baseline_vals, y_vals):
        ymax = max(y)
        ymin = min(y)
        targetvalue = max(abs(ymax - baseline), abs(ymin - baseline))

        diffs = y - baseline
        diffs = diffs / targetvalue
        
        area = np.trapz(diffs)
        if area == 0:
            sts_val = float('nan')
        else:
            sts_val = area

        differences = np.diff(y)
        sum_diff = np.sum(differences)
        sum_abs_diff = np.sum(np.abs(differences))
        change_factor = abs(sum_diff / sum_abs_diff)

        if change_factor == 0:
            sts.append(float('nan'))
        else:
            sts.append(sts_val * change_factor)
    
    return np.array(sts)

def compute_sts(df, features, lfc_cols_list):
    for feature, lfc_cols in zip(features, lfc_cols_list):
        sts_values = calculate_sts_series(df, lfc_cols)
        df[f"{feature}_STS"] = sts_values
    return df


import matplotlib.pyplot as plt

def plot_timecourses_grid(df, features, timepoints, alpha, ylims,outlabel="wangNoRNAsig"):
    metrics = [
        f"{features[0]}_LFC",
        f"{features[1]}_LFC",
        f"{features[2]}_LFC",
    ]
    title_map = {
        f"{features[0]}_LFC": features[0],
        f"{features[1]}_LFC": features[1],
        f"{features[2]}_LFC": features[2],
    }
    color_map = {
        f"{features[0]}_LFC": "#1565C0",
        f"{features[1]}_LFC": "#A4C400",
        f"{features[2]}_LFC": "#00BCD4",
    }
    # df = pd.read_csv(summary_tsv, sep="\t", header=0)

    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(7, 2.5))
    axes = axes.flatten()

    metric_id = 0
    
    for idx, metric in enumerate(metrics):
        #print(metric)
        ax = axes[idx]
        cols = [f"{metric}_{tp}" for tp in timepoints]
        if set(cols) - set(df.columns):
            raise KeyError(f"无法在文件中找到以下列：{set(cols) - set(df.columns)}")
        sub = df[cols].astype(float)
        # data = np.log10(sub + 1)
        data=sub
        

        x = np.arange(len(timepoints))
        base_color = color_map[metric]
        for _, row in data.iterrows():
            if np.any(np.isnan(row.values)):
                pass
            else:
                ax.plot(x, row.values.tolist(), color=base_color, alpha=alpha , lw=0.8)

        arr = data.values
        q1 = np.nanpercentile(arr, 25, axis=0)
        med = np.nanpercentile(arr, 50, axis=0)
        q3 = np.nanpercentile(arr, 75, axis=0)
        mean = np.nanmean(arr, axis=0)

        ax.set_xticks(x)
        ax.set_xticklabels(timepoints, rotation=45, ha='right',fontsize=8)
        ax.tick_params(axis='y', labelsize=8)
        ax.set_xlim(0, len(timepoints) - 1)
        ymin, ymax = np.nanmin(arr), np.nanmax(arr)
        yrange = ymax - ymin if ymax > ymin else 1
        
        ax.set_ylim(ylims[metric_id])
            
        ax.set_title(metric, fontsize=10)
        #ax.set_xlabel('Time', fontsize=8)
        ax.grid(axis='y', linestyle=':', linewidth=0.7, alpha=0.6)
        metric_id += 1

def multi_points_change(df, features, change_marks, thresholds, times, plotit=False, alpha=0.1, ylim=[(-3,3), (-3,3), (-3,3)]):
    for feature, change_mark, threshold in zip(features, change_marks, thresholds):
        if change_mark:
            df = df[df[f"{feature}_STS"].abs() > threshold]
        else:
            df = df[df[f"{feature}_STS"].abs() <= threshold]
    
    if plotit:
        plot_timecourses_grid(df, features, times, alpha, ylim)

def calculate_tsr_series(df, lfc_cols, accumulation_coefficient=0.7):
    time_columns = lfc_cols
    
    tsr_list = []
    for index, row in df.iterrows():
        time_series = row[time_columns].values
        t0h_value = time_series[0]

        topdelta = abs(np.max(time_series) - t0h_value)
        bottomdelta = abs(np.min(time_series) - t0h_value)
        if topdelta >= bottomdelta:
            min_value = np.min(time_series)
            max_value = np.max(time_series)
        else:
            max_value = np.min(time_series)
            min_value = np.max(time_series)
        if (max_value - min_value) != 0:
            normalized_series = (time_series - min_value) / (max_value - min_value)
        else:
            normalized_series = np.zeros_like(time_series)  
        
        area = 0
        for i in range(1, len(normalized_series)):
            area += (normalized_series[i] + normalized_series[i-1]) * (accumulation_coefficient ** (i-1)) / 2
        tsr_list.append(area)
    return np.array(tsr_list)

def compute_tsr(df, features, lfc_cols_list, plotit=False, times=None):
    for feature, lfc_cols in zip(features, lfc_cols_list):
        sts_values = calculate_tsr_series(df, lfc_cols)
        df[f"{feature}_TSR"] = sts_values

        if plotit:
            tmp_df = df.copy()
            tmp_df = tmp_df.filter(like=feature, axis=1)
            tmp_df = tmp_df[tmp_df[f'{feature}_STS'] > 1]
            tmp_df = tmp_df.sort_values(by=f'{feature}_TSR', ascending=False)
            value_cols = []
            for time in times:
                value_cols.append(f"{feature}_LFC_{time}")


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def combine_multi_timepoints(RgxDfs, extra_file=None, extra_name=None, time_points=None):
    #建一个基础的df，后面对它左连接
    rgx1df = pd.read_csv(RgxDfs[0], sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])
    basedf = (
        rgx1df[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol']]
        .drop_duplicates()
        .reset_index(drop=True)
    )

    if time_points==None:
        time_points=[]
        for i in range(len(RgxDfs)):
            time_points.append(f"T{str(i+1)}")

    for rgxfile, time in zip(RgxDfs, time_points):
        rgxdf = pd.read_csv(rgxfile, sep="\t", header=None, names=[
        "peakChr", "peakStart", "peakEnd", "epigenomeActivity",
        "geneSymbol", "geneChr", "geneStart", "geneEnd", "geneStrand",
        "geneID", "weight", "Rgx_rawvalue", "Rgx_percent"
    ])
        rgx_sum = (
            rgxdf.groupby('geneID', as_index=False)
            .agg({'epigenomeActivity': 'sum', 'weight': 'sum', 'Rgx_rawvalue': 'sum'})
            .rename(columns={'epigenomeActivity': f'Enhancer{time}', 'weight': f'Contact{time}', 'Rgx_rawvalue': f'Rgx{time}'})
        )
        basedf = basedf.merge(rgx_sum, on='geneID', how='left')

    if extra_file: 
        extra_time_points=[]
        for i in range(len(RgxDfs)):
            extra_time_points.append(f"{extra_name}T{str(i+1)}")
        names = ["geneChr", "geneStart", "geneEnd", "geneSymbol", "geneID", "geneStrand"] + extra_time_points
        extra_df = pd.read_csv(extra_file, sep='\t', header=None, names=names)
        extra_df = extra_df.drop_duplicates(subset=['geneSymbol'])
        basedf = basedf.merge(extra_df[['geneSymbol'] + extra_time_points], on='geneSymbol', how='left')
        # basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC", f"{extra_name}LFC"]]
    # else:
    # basedf = basedf[['geneChr', 'geneStart', 'geneEnd', 'geneStrand', 'geneID', 'geneSymbol', "EnhancerLFC", "ContactLFC", "RgxLFC"]]

    # resdf, lfc_cols_list = calculate_lfc(basedf, ["Enhancer", "Contact", "Rgx"] + ([extra_name] if extra_file else []), timepoint_counts=len(RgxDfs))
    return basedf

def plot_TSR(df, feature, timepoint_counts, times):
    for i in range(2):
        # df = pd.read_csv("/home/sunpx/my_project/timecourseMultiomics/Tichr/TichrApplication/GSE201376/Advalce_UseMean/TargetGene_PercentAreaMuch.tsv", sep="\t")
        if i == 1:
            tmpdf = df[df[f'{feature}T{timepoint_counts}'] > df[f'{feature}T1']]
            colors = cm.get_cmap('Oranges', 8)(np.linspace(0.2, 0.9, 8))
            outname=f"{feature}_TSR_up.pdf"
        else:
            # print(i)
            # tmpdf = df[df[f'{feature}_STS'] < 0]
            tmpdf = df[df[f'{feature}T{timepoint_counts}'] < df[f'{feature}T1']]
            colors = cm.get_cmap('Greens', 8)(np.linspace(0.2, 0.9, 8))
            outname=f"{feature}_TSR_down.pdf"

        timepoints = []
        for i in range(timepoint_counts):
            timepoints.append(f'{feature}T{i + 1}')
        timepoints_labels = times

        tmpdf[f'{feature}_TSR'] = tmpdf[f'{feature}_TSR'].astype(float)
        df_sorted = tmpdf.sort_values(by=f'{feature}_TSR', ascending=False)

        group_size = len(df_sorted) // 8
        groups = [df_sorted.iloc[i * group_size:(i + 1) * group_size] for i in range(8)]

        plt.figure(figsize=(4, 4))

        for idx, group in enumerate(groups):
            # 计算median_tpm（归一化后）
            median_tpm = group[timepoints].apply(
                lambda x: (x - x.min()) / (x.max() - x.min()) * 100 if x.max() != x.min() else [50]*len(x), axis=1
            ).median(axis=0)
            
            lower_pct = 100- idx * 12.5
            upper_pct = 100- (idx + 1) * 12.5
            label=f'TSR {lower_pct:.1f}-{upper_pct:.1f}%'
            plt.plot(timepoints_labels, median_tpm, label=label, linewidth=2,color=colors[idx])

        #plt.ylim(0, 100)
        #plt.xlabel("Time Points")
        plt.ylabel(f"{feature} relative to (max-min)")
        plt.title(f"Different order of {feature} changes")
        plt.legend(fontsize=9)
        plt.grid(False)
        plt.tight_layout()
        plt.savefig(outname)



def calculate_tsr_series(df, lfc_cols, accumulation_coefficient=0.8):
    time_columns = lfc_cols
    
    tsr_list = []
    for index, row in df.iterrows():
        time_series = row[time_columns].values
        t0h_value = time_series[0]

        topdelta = abs(np.max(time_series) - t0h_value)
        bottomdelta = abs(np.min(time_series) - t0h_value)
        if topdelta >= bottomdelta:
            min_value = np.min(time_series)
            max_value = np.max(time_series)
        else:
            max_value = np.min(time_series)
            min_value = np.max(time_series)
        if (max_value - min_value) != 0:
            normalized_series = (time_series - min_value) / (max_value - min_value)
        else:
            normalized_series = np.zeros_like(time_series)  
        
        area = 0
        for i in range(1, len(normalized_series)):
            area += (normalized_series[i] + normalized_series[i-1]) * (accumulation_coefficient ** (i-1)) / 2
        tsr_list.append(area)
    return np.array(tsr_list)

import re
def compute_tsr(RgxDfs, extra_file=None, extra_name=None, plotit=False, times=None):
    df = combine_multi_timepoints(RgxDfs, extra_file=extra_file, extra_name=extra_name)
    if extra_file == None:
        features = ["Enhancer", "Contact", "Rgx"]
    else:
        features = ["Enhancer", "Contact", "Rgx", extra_name]
    for feature in features:
        columns = [col for col in df.columns.tolist() if re.match(rf'^{feature}T\d+', col)]
        columns.sort(key=lambda x: int(re.search(r'\d+', x).group()))

        sts_values = calculate_tsr_series(df, columns)
        df[f"{feature}_TSR"] = sts_values
        # display(df)

        if plotit:
            plot_TSR(df, feature, len(columns), times)

    return df