bedGraphfile=$1
gt=$2
species=$3
outdir=$4
bedfile=$5
coverageMethod=$6
spmr_flag=$7
refgene_file=$8
tssrange=$9
bashdir=$10

#标记每个位点是否是promoter
awk -v OFS="\t" \
    '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'",$4,$7} \
     $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'",$4,$7}' \
     $refgene_file |sortBed > $outdir/refgene_tss.temp

echo "bedfile $bedfile, bedGraphfile $bedGraphfile"

bedtools intersect -a $bedfile -b $bedGraphfile -wo | awk '
{
    key = $1 "\t" $2 "\t" $3
    result[key] += $7 * $8
}
END {
    for (key in result) {
        print key "\t" result[key]
    }
}
' | sort -k1,1 -k2,2n |\
intersectBed -c -a stdin -b $outdir/refgene_tss.temp > $outdir/candidateSiteCoverage.bdg

if [ "$spmr_flag" == "yes" ]
then
    echo spmr_flag yes
    mv $outdir/candidateSiteCoverage.bdg $outdir/candidateSiteCoverage.bdg.tmp
    scaling_factor=$(echo "scale=8; 1000000/$(awk '{sum += ($3 - $2) * $4} END {print sum}' $bedGraphfile)" | bc)
    awk -v OFS="\t" -v scale=$scaling_factor '{ $4 = $4 * scale; print }' $outdir/candidateSiteCoverage.bdg.tmp > $outdir/candidateSiteCoverage.bdg
fi
