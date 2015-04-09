#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to use annotations in a reference species to estimate those for other
# species, via sequence alignments and Wise2.
#

# these location values may be changed by "make install"
SCRIPTS=..
JARS=..
BIN=../bin
RESOURCES=../resources

GENEWISE=genewise                      # location of the installed GeneWise executable
USER_SEQ=seq.d                         # original user-provided sequences
TEMP=temp.d                            # assorted temporary files
TMP1=$TEMP/per_species                 # used in loops; removed with each pass
TMP2=$TEMP/per_gene
SCRIPT=`echo $0 | sed -e 's;.*/;;'`    # script name from command line; path removed for msgs

#--------------------------
# $TEMP/m_temp
# $TEMP/temp_loc
# $TEMP/temp_gene
# $TEMP/temp_dna
# $TEMP/temp_aa
#
# $TMP1/$sp_name.maf
# $TMP1/$sp_name.loc
# $TMP1/$sp_name.loc_bound
# $TMP1/temp.exons
#
# $TMP2/temp_maf
# $TMP2/m_genes
# $TMP2/cur_prot
# $TMP2/p_seq_temp
# $TMP2/p_seq
# $TMP2/gene_loc
# $TMP2/loc_file
# $TMP2/final_loc_file
# $TMP2/other_temp_dna
# $TMP2/other_temp_aa
#--------------------------

if [ $# -ne 2 ]
then
	echo "Usage:  $SCRIPT refseq_name annot_dir"
	exit 1
fi
ref_sp=$1
ANNOT_DIR=$2

if [ ! -f $USER_SEQ/$ref_sp ]
then
	echo "Sequence matching reference name \"$ref_sp\" is not found in $USER_SEQ."
	exit 1
fi

if [ ! -f $ANNOT_DIR/$ref_sp.codex ]
then
	echo "Reference annotation file \"$ANNOT_DIR/$ref_sp.codex\" is not found."
	exit 1
fi

num_lines=`cat $ANNOT_DIR/$ref_sp.codex | wc -l`
if [ $num_lines -eq 0 ]
then
	echo "Reference annotation file \"$ANNOT_DIR/$ref_sp.codex\" is empty."
	exit 2
fi

rm -rf $TEMP; mkdir -p $TEMP

$BIN/make_ugname $ANNOT_DIR/$ref_sp.codex > $TEMP/m_temp
$BIN/get_gene_bound $TEMP/m_temp > $TEMP/temp_loc
$BIN/pull_c $USER_SEQ/$ref_sp $TEMP/temp_loc > $TEMP/temp_gene    # first nt of seq counts as "1"
$BIN/pull_c $USER_SEQ/$ref_sp $TEMP/m_temp > $TEMP/temp_dna
$BIN/dna2aa -v $TEMP/temp_dna 1 > $TEMP/temp_aa

for sp_name in `ls $USER_SEQ`
do
	if [ ! -f $ANNOT_DIR/$sp_name.codex ]
	then
		echo "creating $ANNOT_DIR/$sp_name.codex"
		#echo "#### Resetting $TMP1 for $sp_name"
		rm -f $TMP1/*    # kluge for weird bug in "rm -rf" (e.g. with AFS on Macs)
		rm -rf $TMP1; mkdir -p $TMP1
		$BIN/lastz $USER_SEQ/$sp_name $TEMP/temp_gene T=2 Y=3400 --format=maf > $TMP1/$sp_name.maf
		$BIN/extract_gene_cluster_whole $TMP1/$sp_name.maf > $TMP1/$sp_name.loc
		$BIN/gene_boundaries $TMP1/$sp_name.loc 1 > $TMP1/$sp_name.loc_bound

		num=0
		exec 0<$TMP1/$sp_name.loc_bound
		while read line
		do
			num=`expr $num + 1`
			#TMP2=$TEMP/per_gene.$num
			#echo "#### Resetting $TMP2 for $line"
			rm -f $TMP2/*    # kluge for weird bug in "rm -rf" (e.g. with AFS on Macs)
			rm -rf $TMP2; mkdir -p $TMP2
			b=`echo $line | cut -d: -f2 | cut -d " " -f2`
			e=`echo $line | cut -d: -f2 | cut -d " " -f3`
			#echo $b $e
			$BIN/lastz "$USER_SEQ/$sp_name[$b,$e]" $TEMP/temp_dna --format=maf > $TMP2/temp_maf
			$BIN/find_match $TMP2/temp_maf > $TMP2/m_genes
			#cat $TMP2/m_genes

			is_done=f
			while read aline
			do
				cur_name=`echo $aline | awk '{print $1}'`
				direction=`echo $aline | awk '{print $2}'`
				#echo "$cur_name $direction"
				if [ $direction = '+' ] || [ $direction = '-' ]
				then
					$BIN/pull_one_prot $TEMP/temp_aa $cur_name > $TMP2/cur_prot
					#cat $TMP2/cur_prot
					num_lines=`cat $TMP2/cur_prot | wc -l`
					num_lines=`expr $num_lines - 1`
					len=`tail -$num_lines $TMP2/cur_prot | wc -c`
					len=`expr $len - $num_lines`
					#echo "len = $len"
				fi

				if [ $direction = '-' ]
				then
					$BIN/dna $b,$e $USER_SEQ/$sp_name > $TMP2/p_seq_temp
					$BIN/dna -c $TMP2/p_seq_temp > $TMP2/p_seq
				elif [ $direction = '+' ]
				then
					$BIN/dna $b,$e $USER_SEQ/$sp_name > $TMP2/p_seq
				fi

				num_lines=0
				cur_lines=0

				if [ "$direction" = "+" ] || [ "$direction" = "-" ]
				then
					#echo "direction: $direction"
					$GENEWISE -genes $TMP2/cur_prot $TMP2/p_seq > $TMP2/gene_loc
					$BIN/wise2sim4 $TMP2/gene_loc $b $e $direction > $TMP2/loc_file

					$BIN/ch_gname $TMP2/loc_file $sp_name$num > $TMP2/final_loc_file
					num_lines=`cat $TMP2/final_loc_file | wc -l`
					#echo "$num_lines"
					if [ "$num_lines" -ne 0 ]
					then
						$BIN/pull_c $USER_SEQ/$sp_name $TMP2/final_loc_file > $TMP2/other_temp_dna
						$BIN/dna2aa -v $TMP2/other_temp_dna 1 > $TMP2/other_temp_aa
						#echo "dna2aa completed successfully"
						#cat $TMP2/other_temp_aa
						cur_lines=`cat $TMP2/other_temp_aa | wc -l`
						cur_lines=`expr $cur_lines - 1`
						cur_len=`tail -$cur_lines $TMP2/other_temp_aa | wc -c`
						#echo "num of lines: $cur_lines; num of chars: $cur_len"
						cur_len=`expr $cur_len - $cur_lines`
						comp=`expr $len \* 7 / 10`
						if [ "$comp" -gt 50 ]
						then
							comp=50    # the average length of exons - about 170 bp in human
						fi

						if [ "$cur_len" -lt "$comp" ]
						then
							#echo "short aa seq $cur_len < $comp"
							tf=f
						else
							tf=`$BIN/filter_out $TMP2/other_temp_aa $USER_SEQ/$sp_name $len`
						fi
					fi
				fi

				if [ "$cur_lines" -eq 0 ] || [ "$num_lines" -eq 0 ]
				then
					#echo "no lines"
					tf=f
				fi

				#echo "tf = $tf"
				if [ $is_done = 't' ]
				then
					:    # do nothing
				elif [ $tf = 't' ]
				then
					cat $TMP2/final_loc_file >> $TMP1/temp.exons
					is_done=t
				elif [ $tf = 'b' ] || [ $tf = 'M' ] || [ $tf = 'P' ]
				then
					$BIN/ext_loc_info $TMP2/final_loc_file "$tf" >> $TMP1/temp.exons
					is_done=t
				fi
			done < $TMP2/m_genes
		done

		if [ ! -f $TMP1/temp.exons ]
		then
			echo "" > $ANNOT_DIR/$sp_name.codex
		else
			$BIN/sort_genes $TMP1/temp.exons > $ANNOT_DIR/$sp_name.codex
		fi
	fi
done

#echo "#### Cleaning up $TEMP"    # rm's '-v' option may help with debugging
rm -f $TMP2/*    # kluge for weird bug in "rm -rf" (e.g. with AFS on Macs)
rm -rf $TMP2
rm -f $TMP1/*
rm -rf $TMP1
rm -f $TEMP/*
rm -rf $TEMP

