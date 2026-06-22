golddf=$1
predictdf=$2
golddf_gene_col=$3
predictdf_gene_col=$4
predictdf_ABC_col=$5
outname=$6
golddf_has_colname=${7:-True}
#prefix=$(openssl rand -base64 6)
prefix=$(head /dev/urandom | tr -dc A-Za-z0-9 | head -c 6)

golddf_colnum=$(awk -F'\t' '{print NF; exit}' $golddf)
awk -v col=$golddf_gene_col '{print $col}' $golddf | sort | uniq > ${prefix}_unique_genes.tmp 
# Column 4 contains the gene IDs in the ground-truth file

awk -v col2=$predictdf_gene_col 'NR==FNR{a[$1];next}($col2 in a)' ${prefix}_unique_genes.tmp $predictdf |sortBed > ${prefix}_predictdf_samegene.tmp 
# Column 18 contains the gene IDs in the prediction file

golddf_genecol_plus=$(expr $golddf_colnum + $predictdf_gene_col)
# Column index of the prediction gene ID after merging

#sed '1d' $golddf |sortBed|\
#	intersectBed -wa -wb -a stdin -b ${prefix}_predictdf_samegene.tmp |\
#        awk -v n=$golddf_gene_col -v m=$golddf_genecol_plus '$n==$m'	> ${prefix}_golddf_predictdf.tmp

if [ "$golddf_has_colname" = True ]; then
	sed '1d' $golddf | sortBed |\
		intersectBed -wa -wb -a stdin -b ${prefix}_predictdf_samegene.tmp |\
		awk -v n=$golddf_gene_col -v m=$golddf_genecol_plus '$n==$m' > ${prefix}_golddf_predictdf.tmp
else
	cat $golddf | sortBed |\
		intersectBed -wa -wb -a stdin -b ${prefix}_predictdf_samegene.tmp |\
		awk -v n=$golddf_gene_col -v m=$golddf_genecol_plus '$n==$m' > ${prefix}_golddf_predictdf.tmp
fi


golddf_abccol_plus=$(expr $golddf_colnum + $predictdf_ABC_col) #Column index of the score after merging

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

python3 ${SCRIPT_DIR}/overlap_predict_true.py ${prefix}_golddf_predictdf.tmp $golddf_colnum $golddf_abccol_plus $outname.tmp
if [ "$golddf_has_colname" = True ]; then
	#echo "golddf_has_colname"
	{ head -1 "$golddf"; echo "matchedScore";} | paste - - |cat - $outname.tmp > $outname
	rm $outname.tmp
else
	#echo "golddf_no_colname"
	mv $outname.tmp $outname
fi

rm ${prefix}*tmp 

