candidatesite_file=$1
readcoverage_file=$2
refgene_file=$3
tssrange=$4
outname=$5

awk -v OFS="\t" \
    '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'",$4,$7} \
     $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'",$4,$7}' \
     $refgene_file |sortBed > refgene_tss.temp

#标记每个位点是否是promoter
sortBed -i $candidatesite_file |\
    coverageBed -a stdin -b $readcoverage_file -counts |\
    intersectBed -c -a stdin -b refgene_tss.temp > $outname

#rm refgene_tss.temp