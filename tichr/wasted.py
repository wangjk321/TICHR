def convert_bedgraph(input_file, output_file, resolution=100):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            chr, start, end, value = line.strip().split()
            start, end = int(start), int(end)
            for i in range(start, end, resolution):
                new_end = min(i + resolution, end)
                outfile.write(f"{chr}\t{i}\t{new_end}\t{value}\n")

input_file = "input.bedgraph"   # 输入的BEDGraph文件
output_file = "output.bedgraph" # 输出的100bp分辨率BEDGraph文件

convert_bedgraph(input_file, output_file)


#原myabc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import average_precision_score
from sklearn.metrics import roc_curve, auc
import seaborn as sns
import os
from preprocess_hic import *
from sys import exit
import time
import subprocess



def processcoverage(codepath,candidatesite_file,read_bamfile_list,refgene_file,tssrange,quantileref,signaltype,separatepromoter,quantile_method):
    peakdf_values = []
    for read_bamfile in read_bamfile_list:
        print("Processing "+str(read_bamfile))
        os.system("bash "+codepath+"/bashcode/processcoverage.sh "+ candidatesite_file+ " " +read_bamfile+ " " + 
                refgene_file+ " " +str(tssrange)+" readcount_oversite.temp")
        peakdf=pd.read_csv("readcount_oversite.temp",sep="\t",header=None)
        peakdf_values.append(list(peakdf[3]))
    peakdf[3] = np.array(peakdf_values).mean(axis=0)

    if quantileref and not signaltype:
        print("Please specify singaltype as DNase or H3K27ac")
        exit(1)
    elif quantileref and signaltype:
        print("Using reference file to make quantiles of ChIP-seq signals")
        refdf = pd.read_csv(quantileref, sep='\t')
        if separatepromoter:
            print("separate promoter")
            peakdf = quantile_normalize_separate(peakdf,refdf,signaltype,method=quantile_method)
        else:
            print("not separate promoter")
            peakdf = quantile_normalize_any(peakdf,refdf,signaltype,method=quantile_method)
    else:
        peakdf = peakdf

    return(peakdf)
        

class calculate_abc:
    def __init__(self,candidatesite_file,read_bamfile_list,candidategene_file,refgene_file,tssrange=500,
                 hicfile=None,hicres=0,hictype="rawhic_sparse",gt=None,juicertool=None,
                 hicnorm="SCALE",hicmindistance=5000,macs2gs=None,blackregion=None,
                quantileref=None,signaltype=None,separatepromoter=False,threads=8,
                read_bamfile2_list=None,signaltype2=None,quantile_method="rankpercent",
                ):
        
        codepath = os.path.dirname(os.path.realpath(__file__))

        '''
        genedf=pd.read_csv(candidategene_file,sep="\t",header=None) #chr,start,end,symbol,value,strands,ID
        genechrlist=genedf[0].unique()
        if hicfile:
            print("Extracting Hi-C information...")
            self.nomhicdf = gethicfile(hicfile,hicres,hictype,genechrlist,hicnorm=hicnorm,gt=gt,juicertool=juicertool,threads=threads)
        '''

        if not candidatesite_file:
            if not blackregion:
                with open("blackregion.blank.tmp", 'w') as file:
                    pass  # 不写入任何内容
                blackregion="blackregion.blank.tmp"
            
            print("Call candidate enhancer sites from bam files")
            subprocess.run(["bash", codepath + "/bashcode/ABC_make_candidate.sh",macs2gs,'MACS2tmp', " ".join(read_bamfile_list), 
                            gt, blackregion,refgene_file], 
                           stdout=open("info.txt", "w"), stderr=subprocess.STDOUT)
            candidatesite_file="MACS2tmp/candidatebed.enhancer.bed"

        print("Extract read count at the given loci")
        peakdf = processcoverage(codepath,candidatesite_file,read_bamfile_list,refgene_file,tssrange,quantileref,signaltype,separatepromoter,quantile_method)
        
        if read_bamfile2_list:
            peakdf2 = processcoverage(codepath,candidatesite_file,read_bamfile2_list,refgene_file,tssrange,quantileref,signaltype2,separatepromoter,quantile_method)
            peakdf[3] = (peakdf[3] * peakdf2[3]) ** 0.5
        
        self.peakdf = peakdf

        print("reading candidate gene files...")
    
        genedf=pd.read_csv(candidategene_file,sep="\t",header=None) #chr,start,end,symbol,value,strands,ID
        genechrlist=genedf[0].unique()
        self.genedf = genedf
        self.contact_data = hicfile
        self.tssrange = tssrange
        self.threads = threads
        self.hicres = int(hicres)
        self.hicmindistance = hicmindistance
        self.hictype = hictype

        os.system("rm readcount_oversite.temp")

        self.given_gamma = 1.024238616787792
        self.given_scale = 5.9594510043736655

        self.ref_gamma = 0.87
        self.ref_scale = -4.80 + 11.63 * 0.87

        if hicfile:
            print("Extracting Hi-C information...")
            self.nomhicdf = gethicfile(hicfile,hicres,hictype,genechrlist,hicnorm=hicnorm,gt=gt,juicertool=juicertool,threads=threads)
        
        self.peakdf[0] = pd.Categorical(self.peakdf[0])

    def getABC_genei(self,genei,maxdistance=5000000,abctype = 'abcep',setpromoter1=True): # type in ['abcep','abcgene']
        if genei % 500 == 0:
            print(str(genei)+" genes processed.")

        strand= self.genedf.iloc[genei,5]
        geneichr = self.genedf.iloc[genei,0]
        if strand == "+":
            tssi = self.genedf.iloc[genei,1]
        elif strand == "-":
            tssi = self.genedf.iloc[genei,2]
        else:
            print("Please confirm the strand of a gene is + or -")
            exit(1)
        
        geneipeaklist = self.peakdf[(self.peakdf[0]==geneichr) & (self.peakdf[2]> tssi-maxdistance) & (self.peakdf[1]< tssi+maxdistance)]
        geneipeakcenterlist = (geneipeaklist[1]+geneipeaklist[2])/2
        ifpromoter = (geneipeaklist[2]> tssi-self.tssrange) & (geneipeaklist[1]< tssi+self.tssrange)
        
        if not self.contact_data:
            geneicontactlist = [getpowerlaw(i-tssi,self.given_gamma,self.given_scale,self.hicmindistance) for i in geneipeakcenterlist]
            geneiabc = np.array(geneicontactlist)*np.array(geneipeaklist[3])
        else:
            tssi_hicindex = int(tssi/self.hicres)
            plfit_contact=[]
            plref_contact=[]
            pseudo_contact=[]
            geneipeaklist_hicindex=[]
            sparseindex = [] #only used in sparse mode
            '''
            pseudo_contact = np.array([getpowerlaw(i-tssi,self.given_gamma,self.given_scale,self.hicres) 
                                       for i in geneipeakcenterlist])
            plref_contact = np.array([getpowerlaw(i-tssi,self.ref_gamma,self.ref_scale,self.hicmindistance) 
                                      for i in geneipeakcenterlist])
            plfit_contact = np.array([getpowerlaw(i-tssi,self.given_gamma,self.given_scale,self.hicmindistance) 
                                      for i in geneipeakcenterlist])
            sparseindex = [str(tssi_hicindex*self.hicres)+"to"+str(i*self.hicres) if tssi_hicindex<= i 
                               else str(i*self.hicres)+"to"+str(tssi_hicindex*self.hicres) 
                               for i in geneipeaklist_hicindex]
            '''

            for x in geneipeakcenterlist:
                pseudo_contact.append(getpowerlaw(x-tssi,self.given_gamma,self.given_scale,self.hicres))
                plref_contact.append(getpowerlaw(x-tssi,self.ref_gamma,self.ref_scale,self.hicmindistance))
                plfit_contact.append(getpowerlaw(x-tssi,self.given_gamma,self.given_scale,self.hicmindistance))
                
                geneipeaklist_hicindex_x = int(x/self.hicres)
                geneipeaklist_hicindex.append(geneipeaklist_hicindex_x)
                if tssi_hicindex <= geneipeaklist_hicindex_x:
                    sparseindex.append(str(tssi_hicindex*self.hicres)+"to"+str(geneipeaklist_hicindex_x * self.hicres))
                else:
                    sparseindex.append(str(geneipeaklist_hicindex_x * self.hicres)+"to"+str(tssi_hicindex*self.hicres))
            

            plfit_contact=np.array(plfit_contact)
            plref_contact=np.array(plref_contact)
            pseudo_contact=np.array(pseudo_contact)

            if self.hictype == "rawhic_dense" or self.hictype == "matrix_dense":
                genei_hicrow = self.nomhicdf[geneichr][tssi_hicindex,:]
                hic_contact = genei_hicrow[geneipeaklist_hicindex]
            elif self.hictype == 'rawhic_sparse':
                hic_contact = self.nomhicdf[geneichr].reindex(sparseindex)["counts"]

            #qc gene
            badgene_threshold = 0.01 #排除没有链接的基因
            if np.nansum(hic_contact) < badgene_threshold:
                hic_contact = plfit_contact
            
            #fillna
            hic_contact = np.nan_to_num(hic_contact)

            geneicontactlist = hic_contact * (plref_contact/plfit_contact) + pseudo_contact
            #geneicontactlist = hic_contact
            geneiabc = np.array(geneicontactlist)*np.array(geneipeaklist[3])


        geneiabcsum = geneiabc.sum()
        if geneiabcsum >0: #判断是否全为0
            geneiabc_percent = geneiabc/geneiabcsum
        else:
            geneiabc_percent = geneiabc * 0 

        if setpromoter1:
            geneiabc_percent[ifpromoter] = 1 #启动子设置为1

        geneiabc_df=geneipeaklist.iloc[:,0:4].copy()
        geneiabc_df.columns = ["peakchr",'peakstart','peakend','activity']
        geneiabc_df["genesymbol"] = self.genedf.iloc[genei,3]
        geneiabc_df["genechr"] = self.genedf.iloc[genei,0]
        geneiabc_df["genestart"] = self.genedf.iloc[genei,1]
        geneiabc_df["geneend"] = self.genedf.iloc[genei,2]
        geneiabc_df["genestrand"] = self.genedf.iloc[genei,5]
        geneiabc_df["geneid"] = self.genedf.iloc[genei,4]
        geneiabc_df["contact"] = list(geneicontactlist)
        geneiabc_df["abc_rawvalue"] = list(geneiabc)
        geneiabc_df["abc_final"] = list(geneiabc_percent)

        if abctype == 'abcep':
            return(geneiabc_df)
        elif abctype == 'abcgene':
            return(geneiabcsum)
    
    def getABC_allgene(self,abctype = 'abcep',setpromoter1=True):
        if abctype == 'abcep':
            abcdf = pd.DataFrame()
            for genei in range(self.genedf.shape[0]):
                abcdf_genei = self.getABC_genei(genei,abctype = 'abcep',setpromoter1=setpromoter1)
                abcdf = pd.concat([abcdf,abcdf_genei], ignore_index=True)
            return(abcdf)
        elif abctype == 'abcgene':
            outdf = self.genedf.copy()
            abcgenelist=[]
            for genei in range(self.genedf.shape[0]):
                abcsum_genei = self.getABC_genei(genei,abctype = 'abcgene',setpromoter1=setpromoter1)
                abcgenelist.append(abcsum_genei)
            outdf['ABCsum_gene']=abcgenelist
            return(outdf)
            

    def getABC_allgene_multiprocess(self,nthread): #多进程--> 不共享内存，适合计算密集型，如果变量很大的话，占用很多内存
        function = self.getABC_genei
        rangeindex = range(self.genedf.shape[0])
        abcdf = pd.concat(multicpu_eachindex(function,rangeindex,threads=nthread), ignore_index=True)
        return(abcdf)
    
    def getABC_allgene_multithread(self,nthread):  #多线程--> 共享内存，适合大量小的运算
        abcdf = pd.DataFrame()
        function = self.getABC_genei
        rangeindex = range(self.genedf.shape[0])
        with concurrent.futures.ThreadPoolExecutor(max_workers=nthread) as executor:
            futures = {executor.submit(function, genei): genei for genei in rangeindex}
            for future in concurrent.futures.as_completed(futures):
                abcdf_genei = future.result()
                abcdf = pd.concat([abcdf,abcdf_genei], ignore_index=True)
        
        return(abcdf)
    


#adjustRgx_old.py
import os
import pandas as pd
import numpy as np
from PRC_ROC import *


def adjrank(matched,tpmfile,tpmcol=[2,3],matchcol="matchedScore",
            matchgene="gene",adjtype="sumrank"):
    tpmdf=pd.read_csv(tpmfile,sep="\t")
    tpmdf["meanTPM"] = tpmdf.iloc[:,tpmcol].mean(axis=1)
    geneTPMdict = tpmdf.set_index("Unnamed: 1")["meanTPM"].to_dict()
    if matchgene != "gene":
        matched["gene"] = matched[matchgene].apply(lambda x: x.split('.')[0])
        matched["geneTPM"] = matched["gene"].map(geneTPMdict).fillna(0)
    else:
        matched["geneTPM"] = matched[matchgene].map(geneTPMdict).fillna(0)
    diffrank = makediffrank(matched[matchcol],matched["geneTPM"])
    sumrank = makesumrank(matched[matchcol],matched["geneTPM"])
    if adjtype == "sumrank":
        return(matched[matchcol]*sumrank)
    elif adjtype == "diffrank":
        return(matched[matchcol]*(1-abs(diffrank)))

def compareAdjust(rgxdir,golddf,outdir,tpmfile,tpmcol,goldcol,
                  withhead=True,title="title",if1stmatch=True,
                  truecol="Significant",matchcol="matchedScore",matchgene="gene",
                 adjtype="sumrank"):
    
    codepath = os.path.dirname(os.path.realpath(__file__))
    print(codepath)

    if not os.path.exists(outdir): os.makedirs(outdir)

    if if1stmatch:
        predictdfname= "combine5"
        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                      rgxdir+predictdfname+"_Rgx_beforematch.tsv"+ 
                      " "+str(goldcol)+" 10 13 "+outdir+predictdfname+"_percent.tsv"+" "+str(withhead))
        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                      rgxdir+predictdfname+"_Rgx_beforematch.tsv"+ 
                      " "+str(goldcol)+" 10 12 "+outdir+predictdfname+"_rawRgx.tsv"+" "+str(withhead))

        predictdfname= "noStructure"
        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                      rgxdir+predictdfname+"_Rgx_beforematch.tsv"+ 
                      " "+str(goldcol)+" 10 13 "+outdir+predictdfname+"_percent.tsv"+" "+str(withhead))
        os.system("bash "+codepath+"/supportcode/overlap_predict_true.sh "+ golddf +" "+
                      rgxdir+predictdfname+"_Rgx_beforematch.tsv"+ 
                      " "+str(goldcol)+" 10 12 "+outdir+predictdfname+"_rawRgx.tsv"+" "+str(withhead))
    
    if withhead:
        nostruc_percent = pd.read_csv(outdir+"noStructure_percent.tsv",sep="\t")
        nostruc_raw = pd.read_csv(outdir+"noStructure_rawRgx.tsv",sep="\t")
        combine5_percent = pd.read_csv(outdir+"combine5_percent.tsv",sep="\t")
        combine5_raw = pd.read_csv(outdir+"combine5_rawRgx.tsv",sep="\t")
    else:
        nostruc_percent = pd.read_csv(outdir+"noStructure_percent.tsv",sep="\t",header=None)
        nostruc_raw = pd.read_csv(outdir+"noStructure_rawRgx.tsv",sep="\t",header=None)
        combine5_percent = pd.read_csv(outdir+"combine5_percent.tsv",sep="\t",header=None)
        combine5_raw = pd.read_csv(outdir+"combine5_rawRgx.tsv",sep="\t",header=None)
    truelabel = nostruc_percent[truecol] == True

    showAUPRC(
        truelabel,nostruc_percent[matchcol],"no-adj",
        truelabel,adjrank(nostruc_percent,tpmfile,tpmcol=tpmcol,matchcol=matchcol,
                          matchgene=matchgene,adjtype=adjtype),"adjRank",
        truelabel,nostruc_raw[matchcol],"adjRaw",
        truelabel,combine5_percent[matchcol],"adjStruc",
        truelabel,adjrank(nostruc_raw,tpmfile,tpmcol=tpmcol,matchcol=matchcol,
                          matchgene=matchgene,adjtype=adjtype),"adjRank+adjRaw",
        truelabel,adjrank(combine5_percent,tpmfile,tpmcol=tpmcol,matchcol=matchcol,
                          matchgene=matchgene,adjtype=adjtype),"adjStruc+adjRank",
        truelabel,combine5_raw[matchcol],"adjStruc+adjRaw",
        truelabel,adjrank(combine5_raw,tpmfile,tpmcol=tpmcol,matchcol=matchcol,
                          matchgene=matchgene,adjtype=adjtype),"adjStruc+adjRank+adjRaw",
        title=title,
             )
    plt.legend(loc="center left",bbox_to_anchor=(1, 0.5))


#makeSurrondingBin.py

def extract_values_from_bedgraph_pandas(df, chromosome, start, end):
    # df 是给定bedgraph文件，
    start=start-1
    #end=end+1

    # 过滤特定染色体的数据
    #！！！！！！！！！！！！！！这一步特别特别慢
    df['chrom'] = pd.Categorical(df['chrom'])
    df_chrom = df[df['chrom'] == chromosome]
    
    # 过滤特定区间的数据
    df_filtered = df_chrom[(df_chrom['end'] >= start) & (df_chrom['start'] <= end)]
    
    # 提取每个bp的值

    result_length = end - start
    values = np.zeros(result_length)

    for _, row in df_filtered.iterrows():
        overlap_start = max(row['start'], start)
        overlap_end = min(row['end'], end)

        # Calculate the relative positions in the result array
        rel_start = overlap_start - start
        rel_end = overlap_end - start

        values[rel_start:rel_end] = row['value']
    
    return values

