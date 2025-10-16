input=$1
outdir=$2
refgene_file=$3
tssrange=$4

awk -v OFS="\t" \
    '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'",$4,$7} \
     $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'",$4,$7}' \
     $refgene_file |sortBed > $outdir/refgene_tss.temp

intersectBed -c -a $input -b $outdir/refgene_tss.temp > $outdir/candidateSiteCoverage.bdg