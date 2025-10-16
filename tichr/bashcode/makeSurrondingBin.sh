candidateGeneFile=$1
peakToGeneMaxDistance=$2
gtfile=$3
binResolution=$4
outdir=$5

awk -v OFS="\t" '$6 == "+" {print $1,$2,$2} $6 == "-" {print $1,$3,$3}' $candidateGeneFile |\
    sortBed | slopBed -i stdin -b $peakToGeneMaxDistance -g $gtfile |\
    awk -v OFS="\t" -v binResolution=$binResolution '{print $1,int($2 / binResolution) * binResolution + binResolution,int($3 / binResolution) * binResolution + binResolution }' |\
    bedtools makewindows -b stdin -w $binResolution | sortBed |uniq > $outdir/surrondingbin.bed
