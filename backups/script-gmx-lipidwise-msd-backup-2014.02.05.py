#!/bin/bash

<<'deprec'
#---NOTE. THIS CODE ONLY WORKS ON ONE STEP. YOU MUST GENERALIZE IT AFTER SORTING OUT THE START/EQUIL/PRODUCTION TIMES
#---settings

basedir=../s8-kraken/
sysinputgro=../s3-start-to-lonestar/system-input.gro
systoploc=../s3-start-to-lonestar/
inputfile=input-md-in.mdp
declare -a systopfiles=(charmm36.ff lipids-tops system.top)
nlipids=80
subjstartres=761
nats=150

#---copy topology files to the current directory for tpr generation
echo "CALCULATING MSD, LIPID-WISE"
for file in ${systopfiles[@]}; do
        cp -r $systoploc/$file ./
done

#---generate a lipid-only trajectory
echo -e "keep 0\nr PI2P\nkeep 1\nq\n" | make_ndx \
	-f $sysinputgro \
	-o index-pi2p.ndx \
&> log-make-ndx-pi2p

#---generate a single lipid-only structure
trjconv \
        -f $sysinputgro \
        -s $(echo `ls -1t $basedir/md.partndex
????.tpr | head -n 1`) \
        -n index-pi2p.ndx \
        -o system-input.pi2p.gro \
&> log-trjconv

#---convert xtc files to lipid-only xtc
for file in $(echo `ls $basedir/md.part????.xtc`); do
	NRUN=$(echo ${file:$((${#basedir}+8)):4} | sed 's/^0*//')
	echo $NRUN
	trjconv \
		-f $(printf $basedir/md.part%04d.xtc $NRUN) \
		-o $(printf md.part%04d.pi2p.xtc $NRUN) \
		-pbc nojump \
		-s $(printf $basedir/md.part%04d.tpr $NRUN) \
		-n index-pi2p.ndx \
	&> $(printf log-trjconv-md.part%04d.log $NRUN)
done
trjcat \
	-f $(echo `ls md.part????.pi2p.xtc`) \
	-o md.parts.msd.pi2p.xtc \
&> log-trjcat
deprec

basetraj=s6-kraken-md.part0018.40000-90000-100.pi2p.xtc
basedir=../s8-kraken #---Nb only necessary to get a single tpr
sig=pi2p
resmolname=PI2P
sysinputgro=../s3-start-to-lonestar/system-input.gro
systoploc=../s3-start-to-lonestar/
inputfile=input-md-in.mdp
declare -a systopfiles=(charmm36.ff lipids-tops system.top)
nlipids=80
nats=150

#---copy topology files to the current directory for tpr generation
echo "CALCULATING MSD, LIPID-WISE"
for file in ${systopfiles[@]}; do
        cp -r $systoploc/$file ./
done

#---generate a lipid-only trajectory
echo -e "keep 0\nr PI2P\nkeep 1\nq\n" | make_ndx \
    -f $sysinputgro \
    -o index-$sig.ndx \
&> log-make-ndx-$sig

#---calculate the starting residue for the selected block of lipids
startatom=$(awk 'NR == 2 {print $1}' index-pi2p.ndx)
resnrline=$(awk -v startat="$startatom" '$0 ~ startat' $sysinputgro | awk -v name="$resmolname" '$0 ~ name {print $0}')
subjstartres=$(echo "$resnrline" | cut -c 1-5)
echo "start lipid is "$subjstartres

#---re-write the topology with only a single lipid
sed '/molecules/q' system.top > system.$sig.top
echo $resmolname" "$nlipids >> system.$sig.top

#---generate a single lipid-only structure
trjconv \
        -f $sysinputgro \
        -s $(echo `ls -1t $basedir/md.part????.tpr | head -n 1`) \
        -n index-$sig.ndx \
        -o system-input.$sig.gro \
&> log-trjconv

#---generate tpr for a single lipid
echo "writing master tpr"
grompp \
	-f $inputfile \
	-c system-input.$sig.gro \
	-p system.$sig.top \
	-o md.parts.msd.$sig.tpr \
	-po mdp-master.mdp \
&> log-grompp-$sig

#---re-write the topology with only a single lipid
sed '/molecules/q' system.top > system.single.top
echo $resmolname" 1" >> system.single.top

#---for each lipid, calculate msd
for nlip in $(eval echo {0..$(($nlipids-1))}); do
	echo "starting lipid "$nlip

	#---identify one lipid and write index file
	echo "writing index"
	echo -e "keep 0\nr "$(($nlip+$subjstartres))"\nkeep 1\nq\n" | make_ndx \
		-f system-input.$sig.gro \
		-o $(printf index-%04d.ndx $nlip) \
	&> $(printf log-make-ndx-%04d.log $nlip)

<<'thiswasstupid'
	#---write structure for single lipid
	echo "writing structure"
	trjconv \
		-f system-input.$sig.gro \
		-s md.parts.msd.$sig.tpr \
		-n $(printf index-%04d.ndx $nlip) \
		-o $(printf conf-%04d.gro $nlip) \
	&> $(printf log-trjconv1-%04d.log $nlip)

	#---generate tpr for a single lipid
	echo "writing tpr"
	grompp \
		-f $inputfile \
		-c $(printf conf-%04d.gro $nlip) \
		-p system.single.top \
		-po $(printf mdpout-%04d.mdp $nlip) \
		-o $(printf conf-%04d.tpr $nlip) \
	&> $(printf log-grompp-%04d.log $nlip)

	#---isolate the single lipid from the lipid-only trajectory
	echo "writing xtc"
	trjconv \
		-f $basetraj \
		-o $(printf conf-%04d.xtc $nlip) \
		-s md.parts.msd.$sig.tpr \
		-n $(printf index-%04d.ndx $nlip) \
		-pbc nojump \
	&> $(printf log-trjconv2-%04d.log $nlip)
	
	#---calculate MSD
	echo "calculating msd"
	echo -e "0\n" | g_msd \
		-f $(printf conf-%04d.xtc $nlip) \
		-o $(printf msd-%04d.xvg $nlip) \
		-s $(printf conf-%04d.tpr $nlip) \
		-lateral z \
	&> $(printf log-msd-%04d.log $nlip)
thiswasstupid

	#---calculate MSD
	echo "calculating msd"
	echo -e "0\n" | g_msd \
		-f $basetraj \
		-n $(printf index-%04d.ndx $nlip) \
		-o $(printf msd-%04d.xvg $nlip) \
		-s md.parts.msd.$sig.tpr \
		-lateral z \
	&> $(printf log-msd-%04d.log $nlip)


done



