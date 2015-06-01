#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to run the CHAP conversion-detection pipeline.
#

# these location values may be changed by "make install"
SCRIPTS=..
JARS=..
BIN=../bin
RESOURCES=../resources

#REPEATMASKER=/srv/gs1/projects/cherry/apps/RepeatMasker/RepeatMasker/RepeatMasker
REPEATMASKER=RepeatMasker
USER_SEQ=seq.d                          # original user-provided sequences
USER_ANNOT=annot.d                      # original user-provided gene annotations
ALL_HDR=$RESOURCES/all_gc.headers.txt   # column headers for conversion output
#OPTION1="one-to-one"                   # not used
OPTION2="one-to-many"                   # parameters for cage_ortho program
OPTION3="many-to-many"
CRIT_BOUND=0.8                          # coverage threshold for choosing conversion criterion
MASK_SEQ=masked_seq.d                   # masked sequences
RM_OUT=rm_out.d                         # RepeatMasker output
CAGE_SELF=ortho.d/self.d                # temporary output from LASTZ
CAGE_PAIR=ortho.d/pair.d                # temporary output from LASTZ
CAGE_CHAIN=ortho.d/ch_pair.d            # full chained pairwise alignments from CAGE
CAGE_DATA=ortho.d/data.d                # consolidated alignments from CAGE
CAGE_OTM=ortho.d/one_to_many_ortho.d    # preliminary "one-to-many" ortholog calls from CAGE
CAGE_MTM=ortho.d/many_to_many_ortho.d   # preliminary "many-to-many" ortholog calls from CAGE
# $sp1.d                                # temp alignments and p-values from conversion detector
CONV_SELF=self.d                        # self alignments from conversion detector
TEMP=temp.d                             # assorted temporary files
CONTIGS=contigs.d                             # assorted temporary contig-info files
SCRIPT=`echo $0 | sed -e 's;.*/;;'`     # script name from command line; path removed for msgs


# ---- get arguments and check tree ----

if [ $# -ne 1 ] && [ $# -ne 2 ] && [ $# -ne 3 ]
then
	echo "Usage:  $SCRIPT species_tree_file [--no_rm]"
	exit 1
fi

if [ ! -d $USER_SEQ ]
then
	echo "Sequence directory not found.  Please create $USER_SEQ and copy sequence files into it."
	exit 1
fi

rm -rf $TEMP; mkdir -p $TEMP
rm -rf $CONTIGS; mkdir -p $CONTIGS
ls $USER_SEQ > $TEMP/all_species.txt
$BIN/check_species_tree $1 $TEMP/all_species.txt

# ---- make working copies of sequence files, and run RepeatMasker ----
# sequence names in each alignment block are formed as $(species name).$(contig name) for multiple contigs

if [ \( "$#" -eq 1 \) -o \( \( "$#" -eq 2 -a "$2" = "--distant" \) \) ]
then
	rm -rf $MASK_SEQ; mkdir -p $MASK_SEQ
	rm -rf $RM_OUT;   mkdir -p $RM_OUT

	for sp in `ls $USER_SEQ`
	do
		num_contigs=`cat $USER_SEQ/$sp | grep ">" | wc -l`
		if [ "$num_contigs" -le 1 ] 
		then
			echo ">"$sp > $TEMP/name
			sed -e '1d' $USER_SEQ/$sp > $TEMP/$sp
			cat $TEMP/name $TEMP/$sp > $MASK_SEQ/$sp
		else 
			$BIN/conv_head $USER_SEQ/$sp $sp > $MASK_SEQ/$sp
		fi	
	done

	for sp in `ls $MASK_SEQ`
	do
		rm -rf $TEMP; mkdir -p $TEMP
		mv $MASK_SEQ/$sp $TEMP/$sp
		$REPEATMASKER -xsmall $TEMP/$sp
		#cp $TEMP/$sp.masked $MASK_SEQ/$sp

		echo ">"$sp > "$TEMP/$sp"_1
		$BIN/remove_headers $USER_SEQ/$sp >> "$TEMP/$sp"_1

    if [ -f $TEMP/$sp.masked ]
    then
      echo ">"$sp > "$TEMP/$sp"_2
      $BIN/remove_headers $TEMP/$sp.masked >> "$TEMP/$sp"_2

      comp=`$BIN/check_mask "$TEMP/$sp"_1 "$TEMP/$sp"_2 | awk '{print $1}'`
      if [ "$comp" == "first" ]
      then
        cat "$TEMP/$sp"_1 > $MASK_SEQ/$sp
      else
        cat $TEMP/$sp.masked > $MASK_SEQ/$sp
      fi
    	rm -rf "$TEMP/$sp"_1 "$TEMP/$sp"_2
			cp $TEMP/$sp.out $RM_OUT/$sp.out
    else
      cat "$TEMP/$sp"_1 > $MASK_SEQ/$sp
    	rm -rf "$TEMP/$sp"_1
			cat ../rm_out.header > $RM_OUT/$sp.out
    fi

	done
elif [ \( $2 = "--no_rm" \) -o \( \( "$#" -eq 3 -a "$3" = "--no_rm" \) \) ]
then
	for sp in `ls $USER_SEQ`
	do
		if [ ! -f $MASK_SEQ/$sp ] || [ ! -f $RM_OUT/$sp.out ]
		then
			echo "Sequence $sp is not masked.  Please remove the --no_rm option."
			exit 1
		fi
	done
else
	echo "Usage:  $SCRIPT species_tree_file [--no_rm]"
	exit 1
fi

# ---- find orthologs ----

rm -rf ortho.d; mkdir -p ortho.d
mkdir $CAGE_SELF
mkdir $CAGE_PAIR
#mkdir ortho.d/ch_self.d
mkdir $CAGE_CHAIN
mkdir $CAGE_DATA
#mkdir ortho.d/$OPTION1-ortho.d
mkdir $CAGE_OTM
mkdir $CAGE_MTM

for sp in `ls $MASK_SEQ`
do
	num_contigs=`cat $MASK_SEQ/$sp | grep ">" | wc -l`
	if [ "$num_contigs" -le 1 ] 
	then
		$BIN/lastz T=2 Y=3400 $MASK_SEQ/$sp --self --nomirror --ambiguous=iupac --format=maf > $TEMP/$sp.maf
	else 
		$BIN/lastz T=2 Y=3400 $MASK_SEQ/$sp[multiple] $MASK_SEQ/$sp --nomirror --ambiguous=iupac --format=maf > $TEMP/$sp.maf
	fi
	$BIN/remove_lower_tri_and_diag $TEMP/$sp.maf > $CAGE_SELF/$sp.maf
done

for sp1 in `ls $MASK_SEQ`
do
	for sp2 in `ls $MASK_SEQ`
	do
		if [ $sp1 != $sp2 ]
		then
			rm -rf $TEMP; mkdir -p $TEMP
			num_contigs=`cat $MASK_SEQ/$sp1 | grep ">" | wc -l`
			if [ $# -eq 1 ] 
			then
				y_drop=3400
			elif [ \( "$2" = "--distant" \) -o \( \( "$#" -eq 3 -a "$3" = "--distant" \) \) ]
			then
				y_drop=5300
			else
				y_drop=3400
			fi

			if [ "$num_contigs" -le 1 ] 
			then
				#BIN/lastz T=2 Y=3400 $MASK_SEQ/$sp1 $MASK_SEQ/$sp2 --ambiguous=iupac --format=maf > $CAGE_PAIR/$sp1.$sp2.maf
				$BIN/lastz T=2 Y=$y_drop $MASK_SEQ/$sp1 $MASK_SEQ/$sp2 --ambiguous=iupac --format=maf > $CAGE_PAIR/$sp1.$sp2.maf
			else
				$BIN/lastz T=2 Y=$y_drop $MASK_SEQ/$sp1[multiple] $MASK_SEQ/$sp2 --ambiguous=iupac --format=maf > $CAGE_PAIR/$sp1.$sp2.maf
			fi

			$BIN/cage_mask $CAGE_SELF/$sp1.maf $CAGE_PAIR/$sp1.$sp2.maf $TEMP/a.maf $TEMP/ab_new.maf 1
			$BIN/cage_mask $CAGE_SELF/$sp2.maf $TEMP/ab_new.maf $TEMP/b.maf $TEMP/new.maf 2
			if [ -f $USER_ANNOT/$sp1.codex ]
			then
				$BIN/cage_chain $TEMP/a.maf $TEMP/new_a.maf $RM_OUT/$sp1.out $RM_OUT/$sp1.out $USER_ANNOT/$sp1.codex $USER_ANNOT/$sp1.codex $sp1 $sp1
			else
				$BIN/cage_chain $TEMP/a.maf $TEMP/new_a.maf $RM_OUT/$sp1.out $RM_OUT/$sp1.out $sp1 $sp1
			fi
			if [ -f $USER_ANNOT/$sp2.codex ]
			then
				$BIN/cage_chain $TEMP/b.maf $TEMP/new_b.maf $RM_OUT/$sp2.out $RM_OUT/$sp2.out $USER_ANNOT/$sp2.codex $USER_ANNOT/$sp2.codex $sp2 $sp2
			else
				$BIN/cage_chain $TEMP/b.maf $TEMP/new_b.maf $RM_OUT/$sp2.out $RM_OUT/$sp2.out $sp2 $sp2
			fi
			if [ -f $USER_ANNOT/$sp1.codex ] && [ -f $USER_ANNOT/$sp2.codex ]
			then
				$BIN/cage_chain $TEMP/new.maf $TEMP/new_ab.maf $RM_OUT/$sp1.out $RM_OUT/$sp2.out $USER_ANNOT/$sp1.codex $USER_ANNOT/$sp2.codex $sp1 $sp2
			else
				$BIN/cage_chain $TEMP/new.maf $TEMP/new_ab.maf $RM_OUT/$sp1.out $RM_OUT/$sp2.out $sp1 $sp2
			fi
			#cp $TEMP/new_a.maf ortho.d/ch_self.d/$sp1.maf
			#cp $TEMP/new_b.maf ortho.d/ch_self.d/$sp2.maf
			cp $TEMP/new_ab.maf $CAGE_CHAIN/$sp1.$sp2.maf
			cat $TEMP/new_a.maf > $CAGE_DATA/$sp1.$sp2.maf
			cat $TEMP/new_b.maf >> $CAGE_DATA/$sp1.$sp2.maf
			cat $TEMP/new_ab.maf >> $CAGE_DATA/$sp1.$sp2.maf
			echo "preliminary orthologous mappings between $sp1 and $sp2"
			#$BIN/cage_ortho $CAGE_DATA/$sp1.$sp2.maf $OPTION1 > ortho.d/$OPTION1-ortho.d/$sp1.$sp2.maf
			$BIN/cage_ortho $CAGE_DATA/$sp1.$sp2.maf $OPTION2 $sp1 $sp2 > $CAGE_OTM/$sp1.$sp2.maf
			$BIN/cage_ortho $CAGE_DATA/$sp1.$sp2.maf $OPTION3 $sp1 $sp2 > $CAGE_MTM/$sp1.$sp2.maf
		fi
	done
done
rm -rf $CAGE_SELF
rm -rf $CAGE_PAIR
#rm -rf ortho.d/ch_self.d           # no longer created
#rm -rf $CAGE_DATA                  # keep for orthology pipeline
#rm -rf ortho.d/$OPTION1-ortho.d    # no longer created

# ensure that all orthologous files are non-empty
for file in `ls $CAGE_OTM`
do
	if [ ! -s $CAGE_OTM/$file ]
	then
		echo "$CAGE_OTM/$file is empty"
		exit 1
	fi
done
for file in `ls $CAGE_MTM`
do
	if [ ! -s $CAGE_MTM/$file ]
	then
		echo "$CAGE_MTM/$file is empty"
		exit 1
	fi
done

# ---- get maximum descent info ----

echo "get maximum descent"
rm -rf $TEMP; mkdir -p $TEMP
rm -rf $CONV_SELF; mkdir -p $CONV_SELF

for sp1 in `ls $MASK_SEQ`
do
	echo $sp1
	rm -rf $sp1.d; mkdir -p $sp1.d
	num_contigs=`cat $MASK_SEQ/$sp1 | grep ">" | wc -l`
	if [ "$num_contigs" -le 1 ] 
	then
		$BIN/lastz $MASK_SEQ/$sp1 --self Q=$RESOURCES/human_dog.q --ambiguous=iupac --maf > $sp1.d/$sp1.maf
	else 
		$BIN/lastz $MASK_SEQ/$sp1[multiple] $MASK_SEQ/$sp1 Q=$RESOURCES/human_dog.q --ambiguous=iupac --maf > $sp1.d/$sp1.maf
	fi
	$BIN/remove_repeats $sp1.d/$sp1.maf > $sp1.d/$sp1.remove_repeats.maf
	$BIN/remove_lower_tri $sp1.d/$sp1.remove_repeats.maf > $sp1.d/$sp1.remove_lower_tri.maf
	$BIN/chain $sp1.d/$sp1.remove_lower_tri.maf > $sp1.d/$sp1.chain
	for sp2 in `ls $MASK_SEQ`
	do
		if [ $sp1 != $sp2 ]
		then
			echo "#$sp1.$sp2" >> $TEMP/all.md
			$BIN/md_quadruplet_entire_region $sp1.d/$sp1.chain $CAGE_OTM/$sp1.$sp2.maf $CAGE_MTM/$sp1.$sp2.maf $CONTIGS/$sp1.contigs.list $CONTIGS/$sp2.contigs.list $CRIT_BOUND >> $TEMP/all.md
			mkdir $sp1.d/$sp2.d    # needed later, for p-values
		fi
	done
	# don't clean up yet; GC_history needs $sp1.d/$sp1.remove_repeats.maf
done

# ---- calculate p-values and call conversions ----

echo "calculate p-values and find conversions"
if [ -f $RESOURCES/p_value_table.txt ]
then
	$BIN/3seq_2D_cluster_fast $TEMP/all.md $RESOURCES/p_value_table.txt > $TEMP/all.pvalue
else
	$BIN/3seq_2D_cluster $TEMP/all.md > $TEMP/all.pvalue
fi
$BIN/separate_pvalue $TEMP/all.pvalue    # output goes to $sp1.d/$sp2.d/all.pvalue

#echo
#echo "======================Results of Conversion=========================="
#echo "Species                 Paralogous pairs          Conversion"
for sp1 in `ls $MASK_SEQ`
do
	parameters=""
	for sp2 in `ls $MASK_SEQ`
	do
		if [ $sp1 != $sp2 ]
		then
			parameters="$parameters $sp1.d/$sp2.d/all.pvalue"
		fi
	done
	$BIN/FDR_cluster $parameters > $TEMP/FDR.txt
	pvalue=`$BIN/FDR_test $TEMP/FDR.txt`
	$BIN/combine_all_pvalue_cluster $parameters $pvalue $TEMP/all.corrected.pvalue
	$BIN/show_GC_cluster $TEMP/all.corrected.pvalue >> $TEMP/all.gc
	rm -f $TEMP/FDR.txt $TEMP/all.corrected.pvalue    # reused
done

# ---- remove redundancies and finalize output ----

rm -f all.gc non-redundant.gc
rm -f all.remove_redundancy.gc all_1.gc species_tree_with_index.txt conversion_in_tree.txt all.gc.tmp
$BIN/GC_history $1 $TEMP/all.gc > all.remove_redundancy.gc
	# also creates all_1.gc, species_tree_with_index.txt, conversion_in_tree.txt
cat species_tree_with_index.txt $ALL_HDR all_1.gc > all.gc.tmp
$BIN/expand_branches all.gc.tmp all.remove_redundancy.gc > all.gc
$BIN/reformat_nonredundant all.gc all.remove_redundancy.gc $CONTIGS > non-redundant.gc

# ---- clean up ----

rm -f all.remove_redundancy.gc all_1.gc all.gc.tmp
rm -f species_tree_with_index.txt    # keep this, as suggested by Federico
rm -f conversion_in_tree.txt          # misleading (based only on upper-bound edges)
for sp1 in `ls $MASK_SEQ`
do
	mv $sp1.d/$sp1.remove_repeats.maf $CONV_SELF    # keep this for Gmaj
	rm -rf $sp1.d
done
rm -rf $TEMP

## ---- print some instructions for viewing results ----
#
#cat << EOF
#
#----------------------------------------------------------------------------------
#
#To view the conversion results with Gmaj (requires Java), run a command like this:
#  $SCRIPTS/gmaj-conv.sh refseq_name
#Example:
#  $SCRIPTS/gmaj-conv.sh human
#
#Or, to view the ortholog calls (also with Gmaj), use:
#  $SCRIPTS/gmaj-ortho.sh refseq_name seq2_name orth_type
#where orth_type is "cage", "context", or "content".  Note that the data files
#for the latter two types will not be available until the ortho.sh script
#has been run.
#Example:
#  $SCRIPTS/gmaj-ortho.sh human dusky_titi cage
#
#These commands have additional options for changing the annotations, sequence
#numbering, etc.; please see the README file for details.
#
#----------------------------------------------------------------------------------
#
#EOF
