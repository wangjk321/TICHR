bamfile=$1
gt=$2
species=$3
outdir=$4
bedfile=$5
coverageMethod=$6
spmr_flag=$7
refgene_file=$8
tssrange=$9

#标记每个位点是否是promoter
awk -v OFS="\t" \
    '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'",$4,$7} \
     $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'",$4,$7}' \
     $refgene_file |sortBed > $outdir/refgene_tss.temp

# bam to bedGraph
if [ "$coverageMethod" == "macs2RP" ]
then
    macs2 callpeak -t $bamfile -f BAM \
        --SPMR -B -q 0.01 --keep-dup 1 --extsize 146 \
        --nomodel -g $species --outdir $outdir -n macs2RP  &> $outdir/macs2RP.log

    sortBed -i $bedfile |\
        bedtools map -a stdin -b $outdir/macs2RP_treat_pileup.bdg -c 4 -o mean -null 0 |\
        intersectBed -c -a stdin -b $outdir/refgene_tss.temp > $outdir/candidateSiteCoverage.bdg
        
elif [ "$coverageMethod" == "coverageBed" ]
then
    sortBed -i $bedfile | coverageBed -a stdin -b $bamfile -counts |\
        intersectBed -c -a stdin -b $outdir/refgene_tss.temp >  $outdir/candidateSiteCoverage.bdg

    if [ "$spmr_flag" == "yes" ]
    then
        echo spmr_flag yes
        mv $outdir/candidateSiteCoverage.bdg $outdir/candidateSiteCoverage.bdg.tmp
        scaling_factor=$(echo "scale=8; 1000000/$(samtools view -c $bamfile)" | bc)
        awk -v OFS="\t" -v scale=$scaling_factor '{ $4 = $4 * scale; print }' $outdir/candidateSiteCoverage.bdg.tmp > $outdir/candidateSiteCoverage.bdg
    fi
fi



