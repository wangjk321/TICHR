import argparse
import os

from .tichr import *
from .context import *
from .siteToGene import *

def main():
    parser = argparse.ArgumentParser(
        description=(
    "TICHR is software to analyse transcriptional regulation "
    "by integrating Epigenome (ChIP-seq etc.), 3D genome (Hi-C) and Transcriptome (RNA-seq).\n\n"
    "See the project page at https://github.com/wangjk321/tichr.\n\n"
    "This command-line provides basic usages; more functions are available within the Python API."
)
                    )
    subparsers = parser.add_subparsers(help="Choose the mode to use sub-commands")

#------------------------------------------------------------------
    #Function1 One line command to calculate Rg and RgX
    def func_calcu(args):
        print("Creating Tichr object...")
        args.readFileList = args.readFileList.split(",")
        if args.readFileList2 != None:
            args.readFileList2 = args.readFileList2.split(",")
        
        #if not os.path.exists(args.outdir): os.makedirs(args.outdir)
        
        tichobj = Tichr(args.candidateSite,args.readFileList,args.gtfile,args.candidateGeneFile,refGeneFile=args.refGeneFile,
                        ifTSSrange=args.TSSrange,S2Gmax=args.S2Gmax,
                        hicfilepath=args.hicfilepath,readFileList2=args.readFileList2)
        print("Start makeSiteBed...")


        tichobj.makeSiteBed(macs2species=args.macs2species,binResolution=args.binResolution,
                    blackregion=args.blackregion,tmpdir=args.tmpdir,fixPeakWidth=args.fixPeakWidth)
        print("Finish makeSiteBed")

        

        print("Start makeSiteBdg...")
        tichobj.makeSiteBdg(args.coverageMethod,spmr = args.spmr,multiBAMmerge=args.multiBAMmerge,file_type=args.file_type,)
        print("Finish makeSiteBdg")
        
        if args.hicfilepath:
            print("Start process HiC...")
            tichobj.proceessHiC(args.hicRes,args.hicDataType,args.hicNormType,juicertool=args.juicertool,
                                threads=args.threads,contactNorm=args.contactNorm,)
            print("Finish process HiC")


        print("Start Computing...")
        tichobj.computeAllGene(args.weightType,fixedFunctionType=args.fixedFunctionType,halfDistance=args.halfDistance,
                               setpromoter1=args.setpromoter1,threads=1,ifUseHiCRef=args.ifUseHiCRef,
                               noise_ratio=args.noise_ratio, noise_quantile=args.noise_quantile)
    

        print("Saving files...")
        tichobj.save(outname=args.outname, header=args.outhead)

        print("Finish Computing...")
        tichobj.clean()




    


        
        

    #input file
    parser_calcu = subparsers.add_parser("calcu", help="Calculate Rg and RgX based on multiomics data")
    input_group = parser_calcu.add_argument_group("Input argument")
    input_group.add_argument("readFileList",help='''A list of input files for epigenomic data such as ChIP-seq, 
                                                   ATAC-seq, or CUT&Tag.\ Supported formats include BAM, BigWig, 
                                                and BedGraph. Multiple files should be provided as a list, e.g., 
                                                ["testdata/DNase_rep1.bam","testdata/DNase_rep2.bam"]''',type=str)
    input_group.add_argument("candidateSite",help='''Candidate regulatory sites, this could be a BED3 file (e.g, "testdata/candidatePeak.bed"),\
                            or a predefined type specified by a string: "denovo_peak","surronding_bin","onlypromoter".''',type=str,)
    input_group.add_argument("candidateGeneFile",help="A tab-separated file listing candidate genes. \
                             It must contain at least six columns in the following order: \
                             [chromosome, start, end, gene symbol, gene ID, strand (+/-)]",type=str)
    input_group.add_argument("gtfile",help="A tab-separated genome table file, with column 1 specifying chromosome names and column 2 indicating chromosome lengths.",type=str)
    input_group.add_argument("--refGeneFile",help="Reference gene file in the same format as candidateGeneFile. This file is used to define gene promoter regions. Typically, it can be the same as candidateGeneFile",type=str)
    input_group.add_argument("--TSSrange",help="Defines the promoter region as transcription start site (TSS) ± this range.",type=int,default=500)
    input_group.add_argument("--S2Gmax",help="Maximum distance (in base pairs) allowed between a peak and a gene for linking",type=int,default=100000)
    input_group.add_argument("--hicfilepath",help="Path to hic files in Juicer .hic format.",type=str,default=None)
    input_group.add_argument("--readFileList2",default=None,help="A second set of epigenome data files, in the same format as readFileList. \
                             For example, DNase signals can be provided in readFileList and H3K27ac signals in readFileList2. \
                             These two signals are combined using the geometric mean",type=str)
    
    

    #process epigenome command
    processEpi_group = parser_calcu.add_argument_group("Process epigenome data arguments")
    processEpi_group.add_argument("--macs2species",help='''Used only when candidateSite is set to "denovo_peak". Specifies the effective genome size for MACS2 peak calling. It can be a numeric value (e.g., 1000000000) or a shortcut string (‘hs’ for human, ‘mm’ for mouse, ‘ce’ for C. elegans, ‘dm’ for Drosophila)''',type=str,default="hs")
    processEpi_group.add_argument("--binResolution",help='''Used only when candidateSite is set to “surronding_bin”. Defines the bin size (in base pairs) for creating windowed candidate sites.''',type=int,default=100)
    processEpi_group.add_argument("--blackregion",help="Regions to exclude, such as ENCODE blacklist sites, provided in Bed3 format",type=str,default=None)
    processEpi_group.add_argument("--fixPeakWidth",dest='fixPeakWidth', action='store_true',help=": Applicable only for given BED3 candidate sites. If set to True, each peak’s width is fixed to 500 bp by centering and extending ±250 b")
    processEpi_group.add_argument("--tmpdir",help="Temporary directory name for intermediate files. Default is a randomly generated name like tichr_tmp_rsDuchihKJ",type=str,default=None)
    processEpi_group.add_argument("--coverageMethod",help='''Method used to compute coverage. For most users, “coverageBed” is recommended.''',type=str,default="coverageBed")
    processEpi_group.add_argument('--spmr', dest='spmr', action='store_true', help="Whether to normalize signal by total mapped reads (Signal Per Million Reads). Set to True if you plan to compare Rg or RgX across samples")
    processEpi_group.add_argument("--multiBAMmerge",type=str,default='mean',help="Strategy to merge multiple replicates. Options: 'mean' (default) or 'sum'.")
    processEpi_group.add_argument("--file_type",help="Format of epigenomic input files. Supported types: “bam”, “bigwig”, or “bedGraph”",default="bam",type=str)

    #process hic command
    processHiC_group = parser_calcu.add_argument_group("Process Hi-C data arguments")
    processHiC_group.add_argument("--hicRes",help="Resolution for Hi-C contact, for example, 10000 for 10kb resolution",type=int,default=25000)
    processHiC_group.add_argument("--hicDataType",type=str,default="rawhic_sparse",help="could be rawhic_sparse (recommended), matrix_dense (dense matrix for each chromosome), \
                                  or rawhic_dense (used for ‘strange’ hic files such as that generated by juicertools >2.0. \
                                  This is the last choice if there are any bugs for the rawhic_sparse mode)",)
    processHiC_group.add_argument("--hicNormType",type=str,default="VC_SQRT",help="Normalization type for Hi-C data. Options: 'KR' (Knight-Ruiz), 'VC' (vanilla coverage), 'VC+S' (vanilla coverage + sparse), 'none' (no normalization).")
    processHiC_group.add_argument("--juicertool",type=str,default=None,help="Path to juicer_tools.jar file. Only for hicDataType=rawhic_dense. Give a user-difined juicertools jar file to process the hic files.")
    processHiC_group.add_argument("--threads",type=int,default=1,help="Number of threads to use for processing Hi-C data. Default is 1.")
    processHiC_group.add_argument("--contactNorm",type=str,default=None,help="default: default normalize; abc: similar normalization to the ABC model; \
                                  oe: observed/expected normalize; 0to1: divide by 95 quantile values; total: divide by the sum of all values, then muliply 1e7. ")
    processHiC_group.add_argument("--ifUseHiCRef",action='store_true',default=False,help="If set, uses the Hi-C reference file to calculate RgX.")

    #Calculation arguments
    calculate_group = parser_calcu.add_argument_group("Calculation arguments")
    calculate_group.add_argument("--weightType",help="Determines how the site-to-gene weight is calculated. Options: \
                                 'hic' (based on Hi-C contact frequency) or 'fixed_function' (based on genomic distance).",default="fixed_function")
    calculate_group.add_argument("--fixedFunctionType",help='''Specifies the function used if weightType=”fixed_function”. Options include: “Sigmoid”, “Exponential”, “Powerlaw”, “NormPL”, “Linear”, “Constant”, “Closest”, or “OnlyPromoter”.''',
                                 default="Exponential",type=str)
    calculate_group.add_argument("--halfDistance",default=100000,type=int,help="Distance (in bp) at which the weight decays to 0.5 for supported functions [sigmoid,exponential,powerlaw,linear-half]")
    calculate_group.add_argument("--setpromoter1",action='store_true',default=False,help="If set, sets the RgX ratio of promoter regions to 1")
    calculate_group.add_argument("--threadscalcu",type=int,default=1,help="Not recommended. Number of threads for calculation")

    

    # outgroup
    out_group = parser_calcu.add_argument_group("Output argument")
    out_group.add_argument("--outname",help="Output name prefix",default="outdir",type=str)
    out_group.add_argument("--outhead",dest='outhead', action='store_true', default=False, help="Define if the output file with column name")
    out_group.add_argument("--noise_ratio",type=float,default=0,help=" If the proportion of the Rgx value to the Rg value of the gene is less than this value, the Rgx value will be set to 0.")
    out_group.add_argument("--noise_quantile",type=float,default=0,help="If the percentile of the Rgx value among all Rgx values of the gene is less than this value, the Rgx value will be set to 0.")


    parser_calcu.set_defaults(func=func_calcu)
# #------------------------------------------------------------------
#     #Function2 DEG analysis for Rg and RgX
#     parser_deg = subparsers.add_parser("deg", help="Differential analysis based on Rg and RgX")

# #------------------------------------------------------------------
#     #Function3 Predict candidate enhancers or target genes
#     parser_ep = subparsers.add_parser("ep", help="Predict candidate enhancers or target genes")


# #------------------------------------------------------------------
#     #Function4 identification of context-specific functions
#     parser_diff = subparsers.add_parser("context", help="Identification of context-specific functions")
#     context_input = parser_diff.add_argument_group("Input file argument for context-specific analysis")
#     context_input.add_argument("type",help="choose a mode from [test,extract]. Use 'test' to exam if a factor has context-specific function. \
#                                Use 'extract' to extract CRM pairs with negative function",default="test",type=str)
#     context_input.add_argument("mergedRgFile",help="Files providing geneRg_ctrl, geneRg_treat, geneTPM, geneLogfc, geneID and geneFDR")
#     context_input.add_argument("--rg_ctrl_col",type=int,)
#     context_input.add_argument("--rg_treat_col",type=int)
#     context_input.add_argument("--tpm_col",type=int)
#     context_input.add_argument("--logfc_col",type=int)
    
#     context_test = parser_diff.add_argument_group("Argument for 'test' analysis")
#     context_test.add_argument("--basedon",default="rg")
    
#     context_extract = parser_diff.add_argument_group("Argument for 'extract' analysis")
#     context_extract.add_argument("--mergedRgxFile",help="Files providing sites(chr,start,end),rgx_geneID_col,rgx_ctrl_col,rgx_treat_col")
#     context_extract.add_argument("--geneid_col",type=int)
#     context_extract.add_argument("--geneFDR_col",type=int)
#     context_extract.add_argument("--rgx_geneID_col",type=int)
#     context_extract.add_argument("--rgx_ctrl_col",type=int)
#     context_extract.add_argument("--rgx_treat_col",type=int)

#     context_output = parser_diff.add_argument_group("Output argument for context-specific analysis")
#     context_output.add_argument("--outname",default="TF")

    def func_adjust(args):

        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)

        #if args.tpmFile and args.strucTypeList and args.strucFileList and args.strucWeightList:
        if args.tpmFile:
            print("Start adjust RgX and Rg...")

            args.tpmFile = os.path.abspath(args.tpmFile)

            args.tpmCols= [int(i) for i in args.tpmCols.split(",")]

            if args.strucTypeList and args.strucFileList and args.strucWeightList:
                args.strucTypeList = args.strucTypeList.split(",")
                print("strucTypeList:",args.strucTypeList)

                args.strucFileList = args.strucFileList.split(",")
                print("strucFileList:",args.strucFileList)

                args.strucWeightList = [float(i) for i in args.strucWeightList.split(",")]
                print("strucWeightList:",args.strucWeightList)

                if len(args.strucTypeList) != len(args.strucFileList) or len(args.strucTypeList) != len(args.strucWeightList):
                    raise ValueError("The length of strucTypeList, strucFileList, and strucWeightList must be the same.")

            adjS2G(args.inputRgx,args.inputRg,args.tpmFile,
                    args.strucTypeList,args.strucFileList,args.strucWeightList,
                    tpmCols=args.tpmCols,tpmHead=args.tpmHead,tpmGeneID=args.tpmGeneID,
                    rankType=args.rankType,outdir=args.outdir,RgXhead=False)
            print("Finish adjust RgX and Rg...")


    # adjust grounp
    parser_adjust = subparsers.add_parser("adjust", help="Adjust S2G regulation")
    adjust_group = parser_adjust.add_argument_group("Adjustment Arguments")

    adjust_group.add_argument("inputRgx",type=str, help="The input RgX DF file, should be standard output from tichr calcu")
    adjust_group.add_argument("inputRg",type=str, help="The input Rg DF file, should be standard output from tichr calcu ")
    adjust_group.add_argument("outdir",type=str, help="The output directory where the adjusted results will be saved.")
    adjust_group.add_argument("--RgXhead",action='store_true',default=False,help="Set it if the RgX files contain header lines.")


    adjust_group.add_argument("--tpmFile",type=str, default=None, help="The TPM file contains gene expression values in a single-column format, \
                              where each row corresponds to the TPM value of a gene.")
    adjust_group.add_argument("--tpmCols",type=str,default="3,4",help="The column number in the TPM file that contains the TPM, could be multiple columns. like 1,2")
    adjust_group.add_argument("--tpmHead",action='store_true',default=False,help="If set, the first row of the TPM file is ignored. \
                              This is useful when the first row contains column headers.")
    adjust_group.add_argument("--tpmGeneID",type=int,default=2,help="The column name in the TPM file that contains gene IDs." )
    adjust_group.add_argument("--rankType",type=str,default=0,help="ranktype could be sumrank or diffrak")

    adjust_group.add_argument("--strucTypeList",type=str,default=None,help="A comma-separated list of structure types to be used for adjustment. \
                              must be supplied in this way 'boundary','tad','loop','stripe','compartmentSame'")
    adjust_group.add_argument("--strucFileList",type=str,default=None,help="A comma-separated list of files containing structure information. \
                              must be supplied in this way 'boundary.bed','tad.bed','loop.bed','stripe.bed','compartmentSame.bed'")
    adjust_group.add_argument("--strucWeightList",type=str,default=None,help="A comma-separated list of weights corresponding to each structure type. \
                              must be supplied in this way 0.5,1.2,5,2,2 ")
    

    parser_adjust.set_defaults(func=func_adjust)

#---------------------------------------

    def func_negative(args):
        print("Merging data frames...")
        rg_merged,rgx_merged = mergeDF(args.rgCtrl,args.rgTreat,args.rgxCtrl,args.rgxTreat,
                                       minRgx=args.minRgx,minRgxRatio=args.minRgxRatio,minRgxQuantile=args.minRgxQuantile)

        print("Computing negative regulation...")

        extractNeg(rg_merged, rgx_merged,showInteration=False,
                corrtype=args.corrtype,filetype="pandas",
                outdir=args.outdir,outname=args.outname,
                geneFC_cutoff=args.geneFC_cutoff,geneFDR_cutoff=args.geneFDR_cutoff,
                rgxFC_cutoff=args.rgxFC_cutoff)

    parser_neg = subparsers.add_parser("neg", help="Identify context-specific repressive functions")
    neg_nes = parser_neg.add_argument_group("Necessary Arguments")

    # neg_nes = neg_input.add_argument_group("Necessary arguments")
    
    neg_nes.add_argument("rgCtrl", type=str,help="Direction to control group Rg file.")
    neg_nes.add_argument("rgTreat", type=str,help="Direction to treated group Rg file.")
    neg_nes.add_argument("rgxCtrl", type=str,help="Direction to control group Rgx file.")
    neg_nes.add_argument("rgxTreat", type=str,help="Direction to treated group Rgx file.")
    neg_nes.add_argument("outdir", type=str,help="Direction for output files.")
    neg_nes.add_argument("outname", type=str,help="Output file prefix.")

    neg_opt = parser_neg.add_argument_group("Options")
    neg_opt.add_argument("--corrtype", type=str, default="pearson", 
                               help="Could be 'spearman' or 'pearson' (Default: pearson)")
    neg_opt.add_argument("--minRgx", type=float, default=0, help="Filter the site-to-gene links by RgX value > min Rgx. (Default: 0)")
    neg_opt.add_argument("--minRgxQuantile", type=float, default=0, 
                         help="If the percentile of the Rgx value among all Rgx values of the gene is less than this value, the Rgx value will be set to 0.")
    neg_opt.add_argument("--minRgxRatio", type=float, default=0, 
                         help="Filter the site-to-gene links by RgX Ratio > min RgxRatio (Default: 0)")
    
    neg_opt.add_argument("--geneFC_cutoff", type=float,default=0.5,
                         help="Cutoff for the absolute log fold-change (|logFC|) of gene expression. Default=0.5")
    neg_opt.add_argument("--geneFDR_cutoff", type=float,default=0.05,
                         help="Cutoff for the FDR of gene expression. Default=0.05")
    neg_opt.add_argument("--rgxFC_cutoff", type=float,default=0.5,
                         help="Cutoff for the absolute log fold-change (|logFC|) of gene-level regulation.")

    parser_neg.set_defaults(func=func_negative)
    

#     def contextfunc(args):
#         if args.type ==  "test":
#             prepare_select_by_rank(args.mergedRgFile, args.rg_ctrl_col, args.rg_treat_col, args.tpm_col, args.logfc_col, 
#                                basedon = args.basedon,label=args.outname)
#         elif args.type ==  "extract":
#             extractNeg(args.mergedRgFile, args.mergedRgxFile, args.rg_ctrl_col, args.rg_treat_col, args.tpm_col, args.logfc_col, 
#                        args.geneid_col,args.geneFDR_col, args.rgx_geneID_col,args.rgx_ctrl_col,args.rgx_treat_col,
#                        negboolRgX=None,iteration=True,showInteration=True,iteration_count=0,outdir="identify_context")

#     parser_diff.set_defaults(func=contextfunc)


# #------------------------------------------------------------------
#     #Function5 large-scale analysis of Rg and RgX
#     parser_large = subparsers.add_parser("large", help="Large-scale analysis of Rg and RgX")

# #------------------------------------------------------------------
#     #Function6 time series analysis of Rg and RgX
#     parser_time = subparsers.add_parser("time", help="Time series analysis of Rg and RgX")

#------------------------------------------------------------------
    parser.add_argument("-V","--version",help="Show tichr version",action='store_true',default=False)
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        print('\nerror: No command specified')
        sys.exit(0)
        
    if args.version:
        print("tichr version 0.1.4")
        exit(0)
    try:
        func = args.func
    except AttributeError:
        parser.error("Too few arguments, please specify more parameters")
    func(args)

    

if __name__ == '__main__':
    main()  
