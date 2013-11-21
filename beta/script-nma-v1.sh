#!/bin/bash

basepath="/home/rpb/worker-big/membrane-repository/membrane-v614-analyze/"
gro="md.part0002.resnr.gro"
xtc="md.part0002.skip10.xtc"
tpr="md.part0002.tpr"
echo -e "keep 0\na 1-358\nkeep 1\nq\n" |  make_ndx \
	-f $basepath$gro \
	-o $basepath"index-protein1.ndx"
g_covar \
	-f $basepath$xtc \
	-s $basepath$tpr \
	-n $basepath"index-protein1.ndx" \
	-av $basepath"average.pdb" \
	-ascii $basepath"covar.dat" \
	-o $basepath"eigenval.xvf" \
	-v $basepath"eigenvec.trr" \
	-l $basepath"covar.log"

