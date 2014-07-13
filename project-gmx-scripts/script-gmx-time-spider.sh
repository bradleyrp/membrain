#!/bin/bash

<<'originalcode'
rm log-gmx-time-spider.log
for file in $(find ./ -name md.part????.edr); do
	gmxcheck -e $file &> tmp.log
	start=$(sed -e 's/\r/\n/g' tmp.log | awk '/Reading energy frame      0 time/ {print $6}')
	endo=$(sed -e 's/\r/\n/g' tmp.log | awk '/Last energy frame read/ {print $7}')
	echo $file" : "$start","$endo >> log-gmx-time-spider.log
done
rm tmp.log
originalcode

masterdir=../
logfilename="log-script-gmx-time-spider-$(echo `date +%Y-%m-%d`).log"
echo "TIME-SPIDER-REPORT" > $logfilename
for sysname in $(find $masterdir -name "membrane-v*" -type d -maxdepth 1); do
	echo "SYSTEM: "$sysname >> $logfilename
	for step in $(find $sysname/ -maxdepth 1 -type d -regex ".*/[s-z][0-9]\-.+"| sort); do
		edrlist=(`ls ${step}/md.part????.edr 2>/dev/null`)
		if [[ $((${#edrlist[@]})) -ne "0" ]]; then
			echo -e "\nPART: "$step >> $logfilename
			if [ -e $step/input-md-in.mdp ]; then
				echo -e "---------------------------------------------" >> $logfilename
				echo "input comparison: " >> $logfilename
				diff input-md-in.mdp $step/input-md-in.mdp >> $logfilename
				echo -e "---------------------------------------------" >> $logfilename
			fi
			for file in $(find $step/ -name md.part????.edr | sort); do
				gmxcheck -e $file &> tmp.log
				start=$(sed -e 's/\r/\n/g' tmp.log | awk '/Reading energy frame      0 time/ {print $6}')
				endo=$(sed -e 's/\r/\n/g' tmp.log | awk '/Last energy frame read/ {print $7}')
				echo "$file	: "$start","$endo >> $logfilename
			done
		fi
	done
	echo -e "" >> $logfilename
done
rm tmp.log


