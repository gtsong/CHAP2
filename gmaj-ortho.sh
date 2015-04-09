#!/bin/sh
set -e    # exit on error (except commands bracketed with "set +e ... set -e")
#
# Script to set up and run Gmaj for the CHAP pipeline, showing
# orthology calls instead of conversion evidence.  It uses the
# exons2underlays tool in the package's bin directory.
#

# these location values may be changed by "make install"
SCRIPTS=..
JARS=..
BIN=../bin
RESOURCES=../resources

JAVA=java
JAR=$JARS/gmaj.jar
E2U="$BIN/exons2underlays -laj -intron PaleGray -fexon"
CAGE_CHAIN=ortho.d/ch_pair.d              # full chained pairwise alignments from CAGE
CAGE_MTM=ortho.d/many_to_many_ortho.d     # preliminary "many-to-many" ortholog calls from CAGE
X_ORTHO=ortho.d/x-ortho.d                 # X-orthology results
N_ORTHO=ortho.d/n-ortho.d                 # N-orthology results
PRE1_ORTHO=temp.d/x-pre1.d                 # N-orthology results
PRE2_ORTHO=temp.d/x-pre2.d                 # N-orthology results
USER_ANNOT=annot.d                        # original user-provided gene annotations
FIG_ANNOT=fig_annot.d                     # gene annotations, including inferred pseudo-genes
UNDERLAYS=temp_underlays.d                # staging area for underlays
GMAJ_PARAM=./gmaj.param                   # parameters file for Gmaj
GMAJ_PREFS=./gmaj.prefs                   # user preferences file for Gmaj
TITLE="Orthology calls vs. all chained alignments"
COLOR=LightGray
HEAP_MB=200                               # memory allowance for Gmaj
SCRIPT=`echo $0 | sed -e 's;.*/;;'`       # script name from command line; path removed for msgs

# ---- get arguments ----

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
	echo "Usage:  $SCRIPT refseq_name seq2_name orth_type [annot_dir]"
	echo "where orth_type is 'cage', 'context', or 'content'"
	exit 1
fi

REFNAME=$1
SEQ2NAME=$2
if [ $3 = 'cage' ]; then
	ORTHO_DIR=$CAGE_MTM
elif [ $3 = 'context' ]; then
	ORTHO_DIR=$X_ORTHO
elif [ $3 = 'content' ]; then
	ORTHO_DIR=$N_ORTHOp
elif [ $3 = 'x-pre1' ]; then
	ORTHO_DIR=$PRE1_ORTHO
elif [ $3 = 'x-pre2' ]; then
	ORTHO_DIR=$PRE2_ORTHO
else
	echo "Usage:  $SCRIPT refseq_name seq2_name orth_type [annot_dir]"
	echo "  where orth_type is 'cage', 'context', or 'content'"
	exit 1
fi
if [ $# -gt 3 ]; then
	ANNOT_DIR=$4
elif [ -d $FIG_ANNOT ]; then
	ANNOT_DIR=$FIG_ANNOT
else
	ANNOT_DIR=$USER_ANNOT
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

if [ -h chained.maf ]; then    # symlink
	rm chained.maf
fi
if [ ! -e chained.maf ]; then
	ln -s $CAGE_CHAIN/$REFNAME.$SEQ2NAME.maf chained.maf
else
	echo "File \"chained.maf\" already exists; please remove or rename it"
	exit 1
fi

if [ -h ortho.maf ]; then    # symlink
	rm ortho.maf
fi
if [ ! -e ortho.maf ]; then
	ln -s $ORTHO_DIR/$REFNAME.$SEQ2NAME.maf ortho.maf
else
	echo "File \"ortho.maf\" already exists; please remove or rename it"
	exit 1
fi

for seq in $REFNAME $SEQ2NAME; do
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

echo "#:gmaj"                                    > $GMAJ_PARAM
echo                                             >> $GMAJ_PARAM
echo "title = \"$TITLE\""                        >> $GMAJ_PARAM
echo "alignfile = ortho.maf chained.maf"         >> $GMAJ_PARAM
echo "refseq = $REFNAME"                         >> $GMAJ_PARAM
echo                                             >> $GMAJ_PARAM
for seq in $REFNAME $SEQ2NAME; do
	echo "seq $seqno:"                          >> $GMAJ_PARAM
	seqno=`expr $seqno + 1`
	echo "seqname   = $seq"                     >> $GMAJ_PARAM
	fn=$ANNOT_DIR/$seq.codex; [ -f $fn ] || fn=
	echo "exons     = $fn"                      >> $GMAJ_PARAM
	fn=$UNDERLAYS/$seq.underlays; ([ -f $fn ] && [ $COLOR ]) || fn=
	echo "underlays = $fn"                      >> $GMAJ_PARAM
	echo                                        >> $GMAJ_PARAM
done

# ---- run Gmaj with extra heap space, and user prefs if present ----

if [ -f $GMAJ_PREFS ]; then
	exec $JAVA -Xmx${HEAP_MB}m -jar $JAR -prefs $GMAJ_PREFS $GMAJ_PARAM
else
	exec $JAVA -Xmx${HEAP_MB}m -jar $JAR $GMAJ_PARAM
fi

