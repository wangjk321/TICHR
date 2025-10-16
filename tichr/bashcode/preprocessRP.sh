bamfile=$1
gt=$2
species=$3

# bam to bedGraph
macs2 callpeak -t $bamfile -f BAM \
    --SPMR -B -q 0.01 --keep-dup 1 --extsize 146 \
    --nomodel -g $species --outdir midfiletmp -n inputbam  &> info.txt

shpath=$(dirname "$(readlink -f "$0")")

#clip
#$shpath/bedClip midfiletmp/inputbam_treat_pileup.bdg $gt midfiletmp/inputbam_treat_pileup.clip.bdg

#bedgraph to biwgwig
$shpath/bedGraphToBigWig midfiletmp/inputbam_treat_pileup.bdg $gt midfiletmp/inputbam_treat_pileup.bw

#for RPEP
enhancersite=$4
refgenebed=$5
if [ -n "$enhancersite" ]; then
    awk -v OFS="\t" \
        '$6 == "+" {print $1,$2-"'$tssrange'",$2+"'$tssrange'",$4,$7} \
        $6 == "-" {print $1,$3-"'$tssrange'",$3+"'$tssrange'",$4,$7}' \
        $refgenebed |sortBed > midfiletmp/refgene_tss.temp

    awk -v OFS='\t' '{print $1,$2,$3,"peak"NR}' $enhancersite |sort -k1,1 -k2,2n -k3,3n > midfiletmp/enhancersite.bed4
    $shpath/bigWigAverageOverBed -bedOut="midfiletmp/candidatepeak.read.bed" midfiletmp/inputbam_treat_pileup.bw midfiletmp/enhancersite.bed4 midfiletmp/read.cov     
    cut -f 1,2,3,5 midfiletmp/candidatepeak.read.bed | intersectBed -c -a stdin -b midfiletmp/refgene_tss.temp > midfiletmp/candidatepeak.read.final
fi
