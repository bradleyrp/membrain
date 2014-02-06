#!/bin/bash

#---method
masterdir=/home/ryb/membrane-v5xx
loc=$(pwd)
tmpname=$(echo ${loc:$((${#masterdir}))})
sysname=$(echo $tmpname | awk 'BEGIN { FS = "/" } ; { print $2}')
timestepps=2
slicestart=40000
sliceend=90000
sliceskip=100

#---flags for subset of system
subset=1
subsetstring="keep 0\nr PI2P\nkeep 1\nq\n"
subsetid=".pi2p"
sysinputgro=../s3-start-to-lonestar/system-input.gro

#-------------------------------------------------------------------------------------------------------------

float_scale=3
function float_eval()
{
    local stat=0
    local result=0.0
    if [[ $# -gt 0 ]]; then
        result=$(echo "scale=$float_scale; $*" | bc -q 2>/dev/null)
        stat=$?
        if [[ $stat -eq 0  &&  -z "$result" ]]; then stat=1; fi
    fi
    echo $result
    return $stat
}
function float_cond()
{
    local cond=0
    if [[ $# -gt 0 ]]; then
        cond=$(echo "$*" | bc -q 2>/dev/null)
        if [[ -z "$cond" ]]; then cond=0; fi
        if [[ "$cond" != 0  &&  "$cond" != 1 ]]; then cond=0; fi
    fi
    local stat=$((cond == 0))
    return $stat
}

#-------------------------------------------------------------------------------------------------------------

echo "TIMESLICE: Finding trajectories in the timeslice "$slicestart","$sliceend","$sliceskip
declare -a filelist=( )
for step in $(find $masterdir/$sysname/ -maxdepth 1 -type d -regex ".*/[s-z][0-9]\-.+"| sort); do
	edrlist=(`ls ${step}/md.part????.edr 2>/dev/null`)
	if [[ $((${#edrlist[@]})) -ne "0" ]]; then
		for file in $(find $step/ -name md.part????.edr | sort); do
			gmxcheck -e $file &> tmp.log
			start=$(sed -e 's/\r/\n/g' tmp.log | awk '/Reading energy frame      0 time/ {print $6}')
			endo=$(sed -e 's/\r/\n/g' tmp.log | awk '/Last energy frame read/ {print $7}')
			if ( $(float_cond $slicestart"<="$endo) && $(float_cond $endo"<="$sliceend) ) \
			|| ( $(float_cond $slicestart"<="$start) && $(float_cond $start"<="$sliceend) ) ; then
				echo $file
				filelist+=(${file:0:-4}.xtc)
			fi
		done
	fi
done
rm tmp.log

#---SUBSET

if [[ $subset -eq 1 ]]; then

echo "TIMESLICE: Generating group file."
echo -e $subsetstring | make_ndx \
	-f $sysinputgro \
	-o index-subset.ndx \
	&> log-make-ndx

declare -a filelist2=( )
echo "TIMESLICE: Creating temporary trajectory slices."
for part in ${filelist[@]}; do
	echo $part
	tmpname=$(echo ${part:$((${#masterdir}+${#sysname}+2)):-4} | sed -ne 's/\//-/p')
	if [[ -f $part ]]; then
		echo "writing "$tmpname".xtc"
		trjconv \
			-f $part \
			-s ${part:0:-4}".tpr" \
			-o $tmpname".xtc" \
			-b $slicestart \
			-n index-subset.ndx \
			-e $sliceend \
			-skip $(float_eval $sliceskip"/"$timestepps) \
			&> log-trjconv-$tmpname
		filelist2+=($tmpname".xtc")
	else
		echo "missing "$part
	fi
done

#-------------------------------------------------------------------------------------------------------------

#---FULL SYSTEM

else

declare -a filelist2=( )
echo "TIMESLICE: Creating temporary trajectory slices."
for part in ${filelist[@]}; do
	echo $part
	tmpname=$(echo ${part:$((${#masterdir}+${#sysname}+2)):-4} | sed -ne 's/\//-/p')
	if [[ -f $part ]]; then
		echo "writing "$tmpname".xtc"
		echo -e "0\n" | trjconv \
			-f $part \
			-s ${part:0:-4}".tpr" \
			-o $tmpname".xtc" \
			-b $slicestart \
			-e $sliceend \
			-skip $(float_eval $sliceskip"/"$timestepps) \
			&> log-trjconv-$tmpname
		filelist2+=($tmpname".xtc")
	else
		echo "missing "$part
	fi
done

fi

#-------------------------------------------------------------------------------------------------------------

echo "TIMESLICE: Concatenating."
basename=${filelist2[0]}
trjcat \
	-f $(echo ${filelist2[@]}) \
	-o ${basename:0:-4}"."$slicestart"-"$sliceend"-"$sliceskip$subsetid".xtc" \
	&> log-trjcat
echo "TIMESLICE: Cleaning up."

for part in ${filelist2[@]}; do
	rm $part
done

#-------------------------------------------------------------------------------------------------------------

if 
