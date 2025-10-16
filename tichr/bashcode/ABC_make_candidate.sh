macs2gs=$1
outdir=$2
bampath=$3
gt=$4
blackregion=$5
refgene_file=$6

peak_extend=250
topn=150000
minPeakWidth=500

tssrange=250

awk -v OFS="\t" \
    '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'"} \
     $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'"}' \
     $refgene_file |sortBed > $outdir/includeTSS.region


macs2 callpeak -f AUTO -g $macs2gs -p 0.1 -n macs2 --shift -75 --extsize 150 \
    --nomodel --keep-dup all --call-summits --outdir $outdir \
    -t $bampath

awk 'BEGIN {{OFS="\t"}} {{if (NF > 0) print $1,"0",$2 ; else print $0}}' $gt > gt.bed.tmp

bedtools intersect -u -a $outdir/macs2_peaks.narrowPeak -b gt.bed.tmp |\
    bedtools sort -faidx $gt -i stdin > $outdir/macs2_peaks.narrowPeak.sorted

cut -f 1-3 $outdir/macs2_peaks.narrowPeak.sorted |\
    coverageBed -a stdin -b $bampath -counts |\
    bedtools sort -i stdin -faidx $gt |\
    bedtools merge -i stdin -c 4 -o max |\
    sort -nr -k 4 | head -n $topn |\
    bedtools intersect -b stdin -a $outdir/macs2_peaks.narrowPeak.sorted -wa |\
    awk '{print $1 "\t" $2 + $10 "\t" $2 + $10}' |\
    bedtools slop -i stdin -b $peak_extend -g $gt |\
    bedtools sort -i stdin -faidx $gt |\
    bedtools merge -i stdin |\
    bedtools intersect -v -wa -a stdin -b $blackregion |\
    cut -f 1-3 |\
    (bedtools intersect -a $outdir/includeTSS.region -b gt.bed.tmp -wa | cut -f 1-3 && cat) |\
    bedtools sort -i stdin -faidx $gt |\
    bedtools merge -i stdin > $outdir/candidatebed.enhancer.bed




#bash /home/wang/github/Tichr/bashcode/ABC_make_candidate.sh hs test \
#   ~/Tichr/2024April/ABC-Enhancer-OnlyPowerLaw/example_chr/chr22/ENCFF860XAE.chr22.sorted.se.bam \
#   ~/Tichr/2024April/ABC-Enhancer-OnlyPowerLaw/reference/hg38/GRCh38_EBV.no_alt.chrom.sizes.tsv \
#   ~/Tichr/2024April/ABC-Enhancer-OnlyPowerLaw/reference/hg38/GRCh38_unified_blacklist.bed \
#   ~/Tichr/2024April/ABC-Enhancer-OnlyPowerLaw/reference/hg38/CollapsedGeneBounds.hg38.bed