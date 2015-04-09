#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to run the CHAP pipeline for detecting conversions and orthology.
#

# these location values may be changed by "make install"
SCRIPTS=..
JARS=..
BIN=../bin
RESOURCES=../resources

CONV_SCRIPT=$SCRIPTS/conversion.sh
FIGURE_SCRIPT=$SCRIPTS/ortho-fig.sh
ANNOT_SCRIPT=infer-annot.sh               # for info message
USER_SEQ=seq.d                            # original user-provided sequences
USER_ANNOT=annot.d                        # original user-provided gene annotations
CAGE_CHAIN=ortho.d/ch_pair.d              # full chained pairwise alignments from CAGE
CAGE_DATA=ortho.d/data.d                  # consolidated alignments from CAGE
CONV_SELF=self.d                          # self alignments from conversion detector
ORTHO_EVENTS=ortho.d/events.d             # event lists from orthology mapper
X_ORTHO=ortho.d/x-ortho.d                 # X-orthology results
N_ORTHO=ortho.d/n-ortho.d                 # N-orthology results
FIG_ANNOT=fig_annot.d                     # gene annotations, including inferred pseudo-genes
TEMP=temp.d                               # assorted temporary files
CONTIGS=contigs.d                             # assorted temporary contig-info files
ORTHO_X1=$TEMP/x-pre1.d
ORTHO_X2=$TEMP/x-pre2.d
ORTHO_N1=$TEMP/n-pre.d
SCRIPT=`echo $0 | sed -e 's;.*/;;'`       # script name from command line; path removed for msgs

if [ $# -ne 2 ] && [ $# -ne 3 ]
then
	echo "Usage:  $SCRIPT species_tree_file refseq_name [--no_rm]"
	exit 1
else
	ref_sp=$2
fi

if [ ! -f $USER_SEQ/$ref_sp ]
then
	echo "Sequence matching reference name \"$ref_sp\" is not found in $USER_SEQ."
	exit 1
fi

num_sp=0       # number of species
num_annot=0    # number of annotation files provided
skip_vis=      # boolean: whether to skip visualization (due to lack of annotations)
if [ ! -d $USER_ANNOT ]
then
	echo "Annotation directory \"$USER_ANNOT\" is not found;"
	echo "creating an empty one to avoid errors."
	mkdir -p $USER_ANNOT
	skip_vis=1
else
	for sp in `ls $USER_SEQ`
	do
		num_sp=`expr $num_sp + 1`
		if [ ! -f $USER_ANNOT/$sp.codex ]
		then
			if [ $sp = $ref_sp ]    # '=' is more portable than '=='
			then
				echo "Reference annotation file \"$USER_ANNOT/$sp.codex\" is not found."
				skip_vis=1
			else
				echo "Annotation file \"$USER_ANNOT/$sp.codex\" is not found."
			fi
		else
			num_annot=`expr $num_annot + 1`
		fi
	done
fi
if [ $num_annot -lt 2 ]
then
	skip_vis=1
fi
if [ $skip_vis ] || [ $num_annot -ne $num_sp ]
then
	if [ $skip_vis ]
	then
		echo
		echo "Generation of summary figures will be skipped, due to insufficient"
		echo "annotations.  You must supply gene annotation files for the reference"
		echo "and at least one other sequence to create gene orthology figures."
	fi
	echo "However, orthologous sequence alignments can still be determined,"
	echo "and should only be slightly affected by the lack of annotations."
	echo
	echo "If you don't know the gene locations in a non-reference species,"
	echo "you may be able to estimate them via the \"$ANNOT_SCRIPT\" script"
	echo "included in the package, which uses GeneWise to infer annotations"
	echo "for other species from the reference (see package documentation)."
	echo
fi

if [ $# -eq 3 ] && [ $3 = "--no_rm" ]
then
	$CONV_SCRIPT $1 --no_rm
else
	$CONV_SCRIPT $1
fi

if [ ! -f non-redundant.gc ]
then
	echo "Conversion output file \"non-redundant.gc\" is not found."
	echo "Please try running the \"$CONV_SCRIPT\" script separately to find the problem."
	exit 1
fi

rm -rf $TEMP; mkdir -p $TEMP

rm -f $TEMP/all.codex; touch $TEMP/all.codex
for fname in `ls $USER_ANNOT`
do
	sp=`echo $fname | cut -d '.' -f1`
	ext=`echo $fname | cut -d '.' -f2`
	if [ "$ext" = "codex" ]
	then
		echo "@ $sp" >> $TEMP/all.codex
		cat $USER_ANNOT/$fname >> $TEMP/all.codex
		echo "" >> $TEMP/all.codex
	fi
done

rm -rf $ORTHO_X1; mkdir -p $ORTHO_X1
rm -rf $ORTHO_N1; mkdir -p $ORTHO_N1

for sp1 in `ls $USER_SEQ`
do
	for sp2 in `ls $USER_SEQ`
	do
		if [ $sp2 != $sp1 ]
		then
			echo "orthologous mappings by context and content between $sp1 and $sp2"
			$BIN/ortho_map $CAGE_DATA/$sp1.$sp2.maf position-ortho non-redundant.gc $TEMP/all.codex > $ORTHO_X1/$sp1.$sp2.maf

			if [ $sp1 = $ref_sp ]
			then
				$BIN/ortho_map $CAGE_DATA/$sp1.$sp2.maf content-ortho non-redundant.gc $TEMP/all.codex $TEMP/$sp1.$sp2.ops $CONTIGS/$sp1.contigs.list $CONTIGS/$sp2.contigs.list > $ORTHO_N1/$sp1.$sp2.maf
			fi
		fi
	done
done

rm -rf $ORTHO_X2; mkdir -p $ORTHO_X2

i=0
for sp1 in `ls $USER_SEQ`
do
	echo "adjusting orthology for $sp1"
	j=0
	for sp2 in `ls $USER_SEQ`
	do
		if [ $j -gt $i ]
		then
			$BIN/adjust_symmetry $ORTHO_X1/$sp1.$sp2.maf $ORTHO_X1/$sp2.$sp1.maf $CAGE_CHAIN/$sp1.$sp2.maf $CAGE_CHAIN/$sp2.$sp1.maf $ORTHO_X2/$sp1.$sp2.maf $ORTHO_X2/$sp2.$sp1.maf $CONTIGS/$sp1.contigs.list $CONTIGS/$sp2.contigs.list
		fi
		j=`expr $j + 1`
	done
	i=`expr $i + 1`
done

mkdir -p $X_ORTHO; rm -f $X_ORTHO/$ref_sp.*.maf
mkdir -p $N_ORTHO; rm -f $N_ORTHO/$ref_sp.*.maf

rm -f $TEMP/combined.ops; touch $TEMP/combined.ops
for sp in `ls $USER_SEQ`
do
	echo "# $sp" >> $CONTIGS/all.contigs.list
	cat $CONTIGS/$sp.contigs.list >> $CONTIGS/all.contigs.list
	if [ $sp != $ref_sp ]
	then
		cat $TEMP/$ref_sp.$sp.ops >> $TEMP/combined.ops
	fi
done
$BIN/adjust_events non-redundant.gc $ref_sp $TEMP/combined.ops $ORTHO_X2 $CONTIGS/all.contigs.list > $TEMP/$ref_sp.events

for sp in `ls $USER_SEQ`
do
	if [ $sp != $ref_sp ]
	then
		echo "refining orthologous mappings between $ref_sp and $sp"
		$BIN/refine_maf $ORTHO_X2/$ref_sp.$sp.maf $CAGE_CHAIN/$ref_sp.$sp.maf $TEMP/$ref_sp.events $TEMP/$ref_sp.$sp.ops $X_ORTHO/$ref_sp.$sp.maf $N_ORTHO/$ref_sp.$sp.maf $CONTIGS/all.contigs.list
	fi
done

mkdir -p $ORTHO_EVENTS

rm -f $ORTHO_EVENTS/$ref_sp.events; touch $ORTHO_EVENTS/$ref_sp.events
cat $TEMP/$ref_sp.events | while read line
do
	event=`echo $line | cut -d ' ' -f1`
	if [ $event = '#' ]
	then
		echo $line >> $ORTHO_EVENTS/$ref_sp.events
	else
		echo $line | cut -d ' ' -f1-5,8 >> $ORTHO_EVENTS/$ref_sp.events
	fi
done

if [ $skip_vis ]
then
#	rm -rf $TEMP
	exit 0
fi
## otherwise:

rm -rf $FIG_ANNOT; mkdir -p $FIG_ANNOT

#------ Giltae will rework this section to make reference-independent annotations ------

cp $USER_ANNOT/$ref_sp.codex $TEMP/$ref_sp.codex
for sp in `ls $USER_SEQ`
do
	if [ $sp != $ref_sp ] && [ -f $USER_ANNOT/$sp.codex ]
	then
		$BIN/find_ps $ORTHO_X2/$sp.$ref_sp.maf $TEMP/$ref_sp.codex $USER_ANNOT/$sp.codex $ref_sp $CONTIGS/all.contigs.list > $TEMP/temp.codex
		mv $TEMP/temp.codex $TEMP/$ref_sp.codex
	fi
done

$BIN/find_ps $CONV_SELF/$ref_sp.remove_repeats.maf $TEMP/$ref_sp.codex $TEMP/$ref_sp.codex $ref_sp $CONTIGS/all.contigs.list > $FIG_ANNOT/$ref_sp.codex

for sp in `ls $USER_SEQ`
do
	if [ $sp != $ref_sp ] && [ -f $USER_ANNOT/$sp.codex ]
	then
		$BIN/find_ps $ORTHO_X2/$ref_sp.$sp.maf $USER_ANNOT/$sp.codex $FIG_ANNOT/$ref_sp.codex $sp $CONTIGS/all.contigs.list > $FIG_ANNOT/$sp.codex
	fi
done

#---------------------------------------------------------------------------------------

# rm -rf $TEMP

$FIGURE_SCRIPT $ref_sp

