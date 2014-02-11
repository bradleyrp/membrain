#!/bin/bash

#---settings
#---Nb the following locations are required, as is the mdp file and topologies.
basetraj=s6-kraken-md.part0018.40000-90000-100.dopc.xtc
basedir=../s8-kraken #---Nb only necessary to get a single tpr
sig=dopc
resmolname=DOPC
subsetstring="keep 0\nr DOPC\nkeep 1\nq\n"
sublipidsel=" & a P"
subsubsetstring="keep 0\nr DOPC"$sublipidsel"\nkeep 1\nq\n"
sysinputgro=../s3-start-to-lonestar/system-input.gro
systoploc=../s3-start-to-lonestar/
syssubsetgro=./system.dopc.gro #---from time-spider
inputfile=input-md-in.mdp
declare -a systopfiles=(charmm36.ff lipids-tops system.top)

#---copy topology files to the current directory for tpr generation
echo "CALCULATING MSD, LIPID-WISE"
for file in ${systopfiles[@]}; do
        cp -r $systoploc/$file ./
done
nlipids=$(awk -v lname="$resmolname " '$0 ~ lname {print $2}' system.top)
echo "nlipids = "$nlipids

#---generate index relative to full system
echo "generating index for the full system"
echo -e $subsetstring | make_ndx \
    -f $sysinputgro \
    -o index-$sig.ndx \
&> log-make-ndx-$sig

#---calculate the starting residue for the selected block of lipids
startatom=$(awk 'NR == 2 {print $1}' index-$sig.ndx)
#---replaced previous method with general one
#resnrline=$(awk -v startat="$startatom" '$0 ~ startat' $sysinputgro | awk -v name="$resmolname" '$0 ~ name {print $0}')
resnrline=$(sed -n -e "/^.\{5\}$(printf %-5s $resmolname).\{5\}$(printf %5s $startatom)/p" ../s3-start-to-lonestar/system-input.gro)
subjstartres=$(echo "$resnrline" | cut -c 1-5)
echo "start lipid is "$subjstartres

#---re-write the topology with only a single lipid
sed '/molecules/q' system.top > system.$sig.top
echo $resmolname" "$nlipids >> system.$sig.top

#---take a subset of the full system
trjconv \
        -f $sysinputgro \
        -s $(echo `ls -1t $basedir/md.part????.tpr | head -n 1`) \
        -n index-$sig.ndx \
        -o system.$sig.gro  \
&> log-trjconv-gro

#---generate tpr for the subset of the full system
echo "writing master tpr"
grompp \
    -f $inputfile \
    -c system.$sig.gro \
    -p system.$sig.top \
    -o md.parts.msd.$sig.tpr \
    -po mdp-master.mdp \
&> log-grompp-$sig

#---generate indices for the subset relative to the subset trajectory
echo -e $subsetstring | make_ndx \
    -f system.$sig.gro \
    -o index-$sig-subset.ndx \
&> log-make-ndx-$sig-subset

#---generate a subset trajectory
trjconv \
        -f $basetraj \
        -s md.parts.msd.$sig.tpr \
        -n index-$sig-subset.ndx \
        -o system.$sig.xtc \
		-pbc nojump \
&> log-trjconv-xtc

#---generate indices for the sub-subset
echo -e $subsubsetstring | make_ndx \
	-f system.$sig.gro \
	-o index-$sig-subsubset.ndx \
&> log-make-ndx-$sig-subset

#---re-write the topology with only a single lipid
sed '/molecules/q' system.top > system.single.top
echo $resmolname" 1" >> system.single.top

#---for each lipid, calculate msd
for nlip in $(eval echo {0..$(($nlipids-1))}); do
	echo "starting lipid "$nlip

	#---generate index for the sub-subset single lipid
	echo "writing index"
	echo -e "keep 0\nr "$(($nlip+$subjstartres))$sublipidsel "\nkeep 1\nq\n" | make_ndx \
		-f system.$sig.gro \
		-o $(printf index-%04d.ndx $nlip) \
	&> $(printf log-make-ndx-%04d.log $nlip)

	#---calculate MSD
	echo "calculating msd"
	g_msd \
		-f system.$sig.xtc \
		-n $(printf index-%04d.ndx $nlip) \
		-o $(printf msd-"$sig"-%04d.xvg $nlip) \
		-s md.parts.msd.$sig.tpr \
		-lateral z \
	&> $(printf log-msd-%04d.log $nlip)
done

#---final calculation with all molecules selected above, sub-subset
echo "calculating MSD for the full system"
echo -e "0\n" | g_msd \
	-f system.$sig.xtc \
	-n index-$sig-subsubset.ndx \
	-o msd-$sig-all.xvg \
	-s md.parts.msd.$sig.tpr \
	-lateral z \
&> log-gmsd-all

