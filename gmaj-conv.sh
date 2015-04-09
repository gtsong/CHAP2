#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to set up and run Gmaj for the gene-conversion pipeline.
# It uses the exons2underlays tool in the package's bin directory.
#

# these location values may be changed by "make install"
SCRIPTS=..
JARS=..
BIN=../bin
RESOURCES=../resources

JAVA=java
JAR=$JARS/gmaj.jar
E2U="$BIN/exons2underlays -laj -intron Clear -fexon"
CAGE_CHAIN=ortho.d/ch_pair.d               # full chained pairwise alignments from CAGE
CONV_SELF=self.d                           # self alignments from conversion detector
GC_FILE=./all.gc                           # output from conversion detector
USER_ANNOT=annot.d                         # original user-provided gene annotations
FIG_ANNOT=fig_annot.d                      # gene annotations, including inferred pseudo-genes
UNDERLAYS=temp_underlays.d                 # staging area for underlays
GMAJ_PARAM=./gmaj.param                    # parameters file for Gmaj
GMAJ_PREFS=./gmaj.prefs                    # user preferences file for Gmaj
HEAP_MB=200                                # memory allowance for Gmaj
SCRIPT=`echo $0 | sed -e 's;.*/;;'`        # script name from command line; path removed for msgs

# ---- get arguments ----

if [ $# -lt 1 ] || [ $# -gt 5 ]; then
	echo "Usage:  $SCRIPT refseq_name [annot_dir] [genomic_offset] [\"title\"] [exon_color]"
	exit 1
fi

REFNAME=$1
if [ -d $FIG_ANNOT ]; then
	ANNOT_DIR=$FIG_ANNOT
else
	ANNOT_DIR=$USER_ANNOT
fi
OFFSET=0
TITLE="Gene conversions in $REFNAME"
COLOR=LightGray

if [ $# -gt 1 ] && [ $2 ]; then
	ANNOT_DIR=$2
fi
if [ $# -gt 2 ]; then
	OFFSET=$3
fi
if [ $# -gt 3 ]; then
	TITLE=$4
fi
if [ $# -gt 4 ]; then
	COLOR=$5
fi
if [ $COLOR ]; then
	if [ $COLOR = "None" ] || [ $COLOR = "none" ] || [ $COLOR = "NONE" ]; then
		COLOR=
	fi
fi

# ---- check availability of Java ----

set +e
which $JAVA > /dev/null 2>&1
if [ $? -ne 0 ]
then
	echo "Program \"$JAVA\" is not found."
	echo "Please install Java and make it available from your command path."
	exit 1
fi
set -e

# ---- set up data files ----

rm -rf $UNDERLAYS; mkdir -p $UNDERLAYS
SEQS=`cat $GC_FILE | grep '^[ \t]*(' -m 1 | sed -e 's;[[][^]]*[]]; ;g' | tr -s '(),;' ' '`

SELFMAF=$CONV_SELF/$REFNAME.remove_repeats.maf
#ORTHMAFS=`echo $CAGE_CHAIN/$REFNAME.*.maf`    # we want tree order instead
ORTHMAFS=""
for seq in $SEQS; do
	if [ $seq != $REFNAME ]; then
		ORTHMAFS="$ORTHMAFS $CAGE_CHAIN/$REFNAME.$seq.maf"
	fi
done
for seq in $SEQS; do
	if [ ! $COLOR ]; then
		:    # no underlays
	elif [ -f $ANNOT_DIR/$seq.underlays ]; then
		echo "Using existing underlay file $ANNOT_DIR/$seq.underlays"
		cp $ANNOT_DIR/$seq.underlays $UNDERLAYS
	elif [ -f $USER_ANNOT/$seq.underlays ]; then
		echo "Using existing underlay file $USER_ANNOT/$seq.underlays"
		cp $USER_ANNOT/$seq.underlays $UNDERLAYS
	elif [ -f $ANNOT_DIR/$seq.codex ]; then
		echo "Creating underlay file $UNDERLAYS/$seq.underlays from $ANNOT_DIR/$seq.codex"
		$E2U $COLOR $ANNOT_DIR/$seq.codex | sed -e 's;\<intron\>;gene;g' \
			> $UNDERLAYS/$seq.underlays
	else
		echo "No annotations found for $seq"
	fi
done

# ---- write parameters file ----

rm -f $GMAJ_PARAM
seqno=0
fn=

echo "#:gmaj"                              > $GMAJ_PARAM
echo                                       >> $GMAJ_PARAM
echo "title     = \"$TITLE\""              >> $GMAJ_PARAM
echo "alignfile = $SELFMAF $ORTHMAFS"      >> $GMAJ_PARAM
echo "refseq    = $REFNAME"                >> $GMAJ_PARAM
echo "geneconv  = $GC_FILE"                >> $GMAJ_PARAM
echo "nowarn    = geneconv_nomaf"          >> $GMAJ_PARAM
echo                                       >> $GMAJ_PARAM
for seq in $REFNAME $REFNAME; do
	echo "seq $seqno:"                    >> $GMAJ_PARAM
	seqno=`expr $seqno + 1`
	echo "seqname   = $seq"               >> $GMAJ_PARAM
	fn=$ANNOT_DIR/$seq.codex; [ -f $fn ] || fn=
	echo "exons     = $fn"                >> $GMAJ_PARAM
	fn=$UNDERLAYS/$seq.underlays; ([ -f $fn ] && [ $COLOR ]) || fn=
	echo "underlays = $fn"                >> $GMAJ_PARAM
	echo "offset    = $OFFSET"            >> $GMAJ_PARAM
	echo                                  >> $GMAJ_PARAM
done
for seq in $SEQS; do
	if [ $seq != $REFNAME ]; then
		echo "seq $seqno:"               >> $GMAJ_PARAM
		seqno=`expr $seqno + 1`
		echo "seqname   = $seq"          >> $GMAJ_PARAM
		fn=$ANNOT_DIR/$seq.codex; [ -f $fn ] || fn=
		echo "exons     = $fn"           >> $GMAJ_PARAM
		fn=$UNDERLAYS/$seq.underlays; ([ -f $fn ] && [ $COLOR ]) || fn=
		echo "underlays = $fn"           >> $GMAJ_PARAM
		echo "offset    = "              >> $GMAJ_PARAM
		echo                             >> $GMAJ_PARAM
	fi
done

# ---- run Gmaj with extra heap space, and user prefs if present ----

if [ -f $GMAJ_PREFS ]; then
	exec $JAVA -Xmx${HEAP_MB}m -jar $JAR -prefs $GMAJ_PREFS $GMAJ_PARAM
else
	exec $JAVA -Xmx${HEAP_MB}m -jar $JAR $GMAJ_PARAM
fi

