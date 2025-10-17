import argparse
import os
from .tichr import *
from .context import *
from .siteToGene import *

def main():
    parser = argparse.ArgumentParser(description="TICHR is software \n \
                                                to analyse transcriptional regulation \n \
                                                by integrating Epigenome (ChIP-seq etc.), \n \
                                                3D genome (Hi-C) and Transcriptome (RNA-seq) \n \
                                                (https://github.com/wangjk321/tichr) ")
    subparsers = parser.add_subparsers(help="Choose the mode to use sub-commands")

#------------------------------------------------------------------
    #Function1 One line command to calculate Rg and RgX
    def func_calcu(args):
        print("Creating Tichr object...")
        args.readFileList = args.readFileList.split(",")
        args.readFileList2 = args.readFileList2.split(",")
        
        if not os.path.exists(args.outdir): os.makedirs(args.outdir)
        
        tichobj = Tichr(args.candidateSite,args.readFileList,args.gtfile,args.candidateGeneFile,refgene_file=args.refgene_file,
                        ifTSSrange=args.TSSrange,peakToGeneMaxDistance=args.peakToGeneMaxDistance,
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
                                threads=args.threads,further_normalize_type=args.further_normalize_type,)
            print("Finish process HiC")


        print("Start Computing...")
        tichobj.computeAllGene(args.weightType,fixedFunctionType=args.fixedFunctionType,halfDistance=args.halfDistance,
                               setpromoter1=args.setpromoter1,threads=1,ifUseHiCRef=args.ifUseHiCRef,)
        
        tichobj.RgxDf.to_csv(args.outdir + "/RgX.tsv",header=None,sep="\t",index=None)
        tichobj.RgDF.to_csv(args.outdir + "/Rg.tsv",header=None,sep="\t",index=None)
        print("Finish Computing...")


        if args.tpmfile and args.structureTypeList and args.structureFileList and args.structureWeightList:
            print("Start adjust RgX and Rg...")
            args.tpmfile = os.path.abspath(args.tpmfile)
            args.structureTypeList = args.structureTypeList.split(",")
            args.structureFileList = args.structureFileList.split(",")
            args.structureWeightList = [float(i) for i in args.structureWeightList.split(",")]
            if len(args.structureTypeList) != len(args.structureFileList) or len(args.structureTypeList) != len(args.structureWeightList):
                raise ValueError("The length of structureTypeList, structureFileList, and structureWeightList must be the same.")
            
            args.tmpcolrep= [int(i) for i in args.tmpcolrep.split(",")]
            adjvalue(args.outdir + "/RgX.tsv",args.outdir + "/Rg.tsv",args.outdir,args.tpmfile,
                    args.structureTypeList,args.structureFileList,args.structureWeightList,
                    tmpcolrep=args.tmpcolrep,ignorehead=args.ignorehead,tmpgeneID=args.tmpgeneID,
                    ranktype=args.ranktype,)
            print("Finish adjust RgX and Rg...")
        

    #input file
    parser_calcu = subparsers.add_parser("calcu", help="Calculate Rg and RgX based on multiomics data")
    input_group = parser_calcu.add_argument_group("Input and output argument")
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
    input_group.add_argument("--refgene_file",help="Reference gene file in the same format as candidateGeneFile. This file is used to define gene promoter regions. Typically, it can be the same as candidateGeneFile",type=str)
    input_group.add_argument("--TSSrange",help="Defines the promoter region as transcription start site (TSS) ± this range.",type=int,default=500)
    input_group.add_argument("--peakToGeneMaxDistance",help="Maximum distance (in base pairs) allowed between a peak and a gene for linking",type=int,default=100000)
    input_group.add_argument("--hicfilepath",help="Path to hic files in Juicer .hic format.",type=str)
    input_group.add_argument("--readFileList2",help="A second set of epigenome data files, in the same format as readFileList. \
                             For example, DNase signals can be provided in readFileList and H3K27ac signals in readFileList2. \
                             These two signals are combined using the geometric mean",type=str,default=None)
    input_group.add_argument("--outdir",help="Output directory",default="outdir",type=str)

    #process epigenome command
    processEpi_group = parser_calcu.add_argument_group("Process epigenome data arguments")
    processEpi_group.add_argument("--macs2species",help='''Used only when candidateSite is set to "denovo_peak". Specifies the effective genome size for MACS2 peak calling. It can be a numeric value (e.g., 1000000000) or a shortcut string (‘hs’ for human, ‘mm’ for mouse, ‘ce’ for C. elegans, ‘dm’ for Drosophila)''',type=str,default="hs")
    processEpi_group.add_argument("--binResolution",help='''Used only when candidateSite is set to “surronding_bin”. Defines the bin size (in base pairs) for creating windowed candidate sites.''',type=int,default=100)
    processEpi_group.add_argument("--blackregion",help="Regions to exclude, such as ENCODE blacklist sites, provided in Bed3 format",type=str,default=None)
    processEpi_group.add_argument("--fixPeakWidth",help=": Applicable only for given BED3 candidate sites. If set to True, each peak’s width is fixed to 500 bp by centering and extending ±250 b",type=int,default=None)
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
    processHiC_group.add_argument("--further_normalize_type",type=str,default=None,help="default: default normalize; abc: similar normalization to the ABC model; \
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

    #Adjustment arguments
    
    # adjust grounp
    adjust_group = parser_calcu.add_argument_group("Adjustment arguments")
    adjust_group.add_argument("--tpmfile",type=str, default=None, help="The TPM file contains gene expression values in a single-column format, \
                              where each row corresponds to the TPM value of a gene.")
    adjust_group.add_argument("--tmpcolrep",type=str,default="3,4",help="The column number in the TPM file that contains the TPM, could be multiple columns. like 1,2")
    adjust_group.add_argument("--ignorehead",action='store_true',default=False,help="If set, the first row of the TPM file is ignored. \
                              This is useful when the first row contains column headers.")
    adjust_group.add_argument("--tmpgeneID",type=int,default=2,help="The column name in the TPM file that contains gene IDs." )
    adjust_group.add_argument("--structureTypeList",type=str,default=None,help="A comma-separated list of structure types to be used for adjustment. \
                              must be supplied in this way 'boundary','tad','loop','stripe','compartmentSame'")
    adjust_group.add_argument("--structureFileList",type=str,default=None,help="A comma-separated list of files containing structure information. \
                              must be supplied in this way 'boundary.bed','tad.bed','loop.bed','stripe.bed','compartmentSame.bed'")
    adjust_group.add_argument("--structureWeightList",type=str,default=None,help="A comma-separated list of weights corresponding to each structure type. \
                              must be supplied in this way 0.5,1.2,5,2,2 ")
    adjust_group.add_argument("--ranktype",type=int,default=0,help="ranktype could be sumrank or diffrak")

    parser_calcu.set_defaults(func=func_calcu)


#------------------------------------------------------------------
    #Function2 DEG analysis for Rg and RgX
    parser_deg = subparsers.add_parser("deg", help="Differential analysis based on Rg and RgX")

#------------------------------------------------------------------
    #Function3 Predict candidate enhancers or target genes
    parser_ep = subparsers.add_parser("ep", help="Predict candidate enhancers or target genes")


#------------------------------------------------------------------
    #Function4 identification of context-specific functions
    parser_diff = subparsers.add_parser("context", help="Identification of context-specific functions")
    context_input = parser_diff.add_argument_group("Input file argument for context-specific analysis")
    context_input.add_argument("type",help="choose a mode from [test,extract]. Use 'test' to exam if a factor has context-specific function. \
                               Use 'extract' to extract CRM pairs with negative function",default="test",type=str)
    context_input.add_argument("mergedRgFile",help="Files providing geneRg_ctrl, geneRg_treat, geneTPM, geneLogfc, geneID and geneFDR")
    context_input.add_argument("--rg_ctrl_col",type=int,)
    context_input.add_argument("--rg_treat_col",type=int)
    context_input.add_argument("--tpm_col",type=int)
    context_input.add_argument("--logfc_col",type=int)
    
    context_test = parser_diff.add_argument_group("Argument for 'test' analysis")
    context_test.add_argument("--basedon",default="rg")
    
    context_extract = parser_diff.add_argument_group("Argument for 'extract' analysis")
    context_extract.add_argument("--mergedRgxFile",help="Files providing sites(chr,start,end),rgx_geneID_col,rgx_ctrl_col,rgx_treat_col")
    context_extract.add_argument("--geneid_col",type=int)
    context_extract.add_argument("--geneFDR_col",type=int)
    context_extract.add_argument("--rgx_geneID_col",type=int)
    context_extract.add_argument("--rgx_ctrl_col",type=int)
    context_extract.add_argument("--rgx_treat_col",type=int)

    context_output = parser_diff.add_argument_group("Output argument for context-specific analysis")
    context_output.add_argument("--outname",default="TF")

    def contextfunc(args):
        if args.type ==  "test":
            prepare_select_by_rank(args.mergedRgFile, args.rg_ctrl_col, args.rg_treat_col, args.tpm_col, args.logfc_col, 
                               basedon = args.basedon,label=args.outname)
        elif args.type ==  "extract":
            extractNeg(args.mergedRgFile, args.mergedRgxFile, args.rg_ctrl_col, args.rg_treat_col, args.tpm_col, args.logfc_col, 
                       args.geneid_col,args.geneFDR_col, args.rgx_geneID_col,args.rgx_ctrl_col,args.rgx_treat_col,
                       negboolRgX=None,iteration=True,showInteration=True,iteration_count=0,outdir="identify_context")

    parser_diff.set_defaults(func=contextfunc)


#------------------------------------------------------------------
    #Function5 large-scale analysis of Rg and RgX
    parser_large = subparsers.add_parser("large", help="Large-scale analysis of Rg and RgX")

#------------------------------------------------------------------
    #Function6 time series analysis of Rg and RgX
    parser_time = subparsers.add_parser("time", help="Time series analysis of Rg and RgX")

#------------------------------------------------------------------
    parser.add_argument("-V","--version",help="Show tichr version",action='store_true',default=False)
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        print('\nerror: No command specified')
        sys.exit(0)
        
    if args.version:
        print("tichr version 0.0.2")
        exit(0)
    try:
        func = args.func
    except AttributeError:
        parser.error("Too few arguments, please specify more parameters")
    func(args)

if __name__ == '__main__':
    main()  