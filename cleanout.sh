#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to clean out a cluster directory, removing files and directories
# added by the CHAP pipeline.  Any material added by users to pipeline-created
# directories will be wiped out when the directories are removed, but other
# user files will be left alone.
#
# The level parameter controls which files and directories are removed, with
# higher values specifying increasingly drastic cleanup.  They are cumulative,
# with each level including all lower ones.
#
#  0: temporary scratch files normally deleted automatically by the pipeline;
#     useful if the script did not finish due to an error
#
#  1: additional intermediate output from pipeline programs; however
#     final results and files needed by Gmaj and figure generator are kept
#
#  2: result files for all reference sequences other than the specified one
#
#  3: all output except RepeatMasker results and Gmaj user preferences; useful
#     with the pipeline's "--no_rm" option to avoid the delay of re-masking
#
#  4: all output; only the original user input should remain
#

USER_SEQ=seq.d                            # original user-provided sequences (always keep)
USER_ANNOT=annot.d                        # original user-provided gene annotations (always keep)
MASK_SEQ=masked_seq.d                     # 4: masked sequences
RM_OUT=rm_out.d                           # 4: RepeatMasker output
ORTHO_ALL=ortho.d                         # 3: container for orthology files
CAGE_SELF=ortho.d/self.d                  # 0: temporary output from LASTZ
CAGE_PAIR=ortho.d/pair.d                  # 0: temporary output from LASTZ
CAGE_CHAIN=ortho.d/ch_pair.d              # 2: full chained pairwise alignments from CAGE
CAGE_DATA=ortho.d/data.d                  # 1: consolidated alignments from CAGE
CAGE_OTM=ortho.d/one_to_many_ortho.d      # 1: preliminary "one-to-many" ortholog calls from CAGE
CAGE_MTM=ortho.d/many_to_many_ortho.d     # 2: preliminary "many-to-many" ortholog calls from CAGE
# $sp.d                                   # 0: temporary alignments, etc. from conversion detector
CONV_SELF=self.d                          # 2: self alignments from conversion detector
GC_TEMP="all.remove_redundancy.gc all_1.gc conversion_in_tree.txt all.gc.tmp"    # 0: temp files
GC_OUT="all.gc non-redundant.gc species_tree_with_index.txt"                     # 3: conv output
ORTHO_EVENTS=ortho.d/events.d             # 2: event lists from orthology mapper
X_ORTHO=ortho.d/x-ortho.d                 # 2: X-orthology results
N_ORTHO=ortho.d/n-ortho.d                 # 2: N-orthology results
FIG_ANNOT=fig_annot.d                     # 3: gene annotations, including inferred pseudo-genes
FIGURES=figures.d                         # 3: files for automatic PostScript figures
FIG_IN="x-ortho.fig n-ortho.fig"          # 1: suffixes for orthology figure input
UNDERLAYS=temp_underlays.d                # 1: staging area for underlays in gmaj* scripts
GMAJ_SYM="chained.maf ortho.maf"          # 1: symlinks for alignfiles in gmaj-ortho.sh
GMAJ_PARAM=gmaj.param                     # 1: parameters file for Gmaj
GMAJ_PREFS=gmaj.prefs                     # 4: user preferences file for Gmaj
TEMP=temp.d                               # 0: assorted temporary files
SCRIPT=`echo $0 | sed -e 's;.*/;;'`       # script name from command line; path removed for msgs


# ---- get arguments ----

if [ $# -lt 2 ] || [ $# -gt 3 ]; then
	echo "Usage:  $SCRIPT cluster_dir clean_level [refseq_name]"
	exit 1
fi
if [ ! -d $1 ]; then
	echo "Directory \"$1\" is not found."
	exit 1
else
	cd $1
fi
if [ $2 -lt 0 ] || [ $2 -gt 4 ]; then
	echo "Level must be an integer in the range 0-4."
	exit 1
else
	LEVEL=$2
fi
if [ $LEVEL -eq 2 ]; then
	if [ $# -lt 3 ]; then
		echo "Refseq_name is required for level 2."
		exit 1
	else
		REFNAME=$3
	fi
else
	if [ $# -gt 2 ]; then
		echo "Refseq_name only applies for level 2."
		exit 1
	else
		REFNAME=
	fi
fi

# ---- clean up by level ----

if [ $LEVEL -ge 0 ]; then
	rm -rf $TEMP
	rm -rf $CAGE_SELF
	rm -rf $CAGE_PAIR
	if [ -d $USER_SEQ ]; then
		for sp in `ls $USER_SEQ`; do
			rm -rf $sp.d
		done
	fi
	rm -f $GC_TEMP
fi

if [ $LEVEL -ge 1 ]; then
	rm -rf $CAGE_DATA
	rm -rf $CAGE_OTM
	for suff in $FIG_IN; do
		rm -f $FIGURES/*.$suff
	done
	rm -rf $UNDERLAYS
	for sym in $GMAJ_SYM; do
		if [ -h $sym ]; then    # symlink
			rm $sym
		fi
	done
	rm -f $GMAJ_PARAM
fi

if [ $LEVEL -eq 2 ]; then    # we have $REFNAME
	if [ -d $USER_SEQ ]; then
		for sp in `ls $USER_SEQ`; do
			if [ $sp != $REFNAME ]; then
				rm -f $CAGE_CHAIN/$sp.*
				rm -f $CAGE_MTM/$sp.*
				rm -f $CONV_SELF/$sp.*
				rm -f $ORTHO_EVENTS/$sp.*
				rm -f $X_ORTHO/$sp.*
				rm -f $N_ORTHO/$sp.*
				rm -f $FIGURES/$sp.*
			fi
		done
	fi
fi

if [ $LEVEL -ge 3 ]; then
	rm -rf $ORTHO_ALL
#	rm -rf $CAGE_CHAIN      # included in $ORTHO_ALL
#	rm -rf $CAGE_MTM        # included in $ORTHO_ALL
	rm -rf $CONV_SELF
	rm -f $GC_OUT
#	rm -rf $ORTHO_EVENTS    # included in $ORTHO_ALL
#	rm -rf $X_ORTHO         # included in $ORTHO_ALL
#	rm -rf $N_ORTHO         # included in $ORTHO_ALL
	rm -rf $FIG_ANNOT
	rm -rf $FIGURES
fi

if [ $LEVEL -ge 4 ]; then
	rm -rf $MASK_SEQ
	rm -rf $RM_OUT
	rm -f $GMAJ_PREFS
fi

