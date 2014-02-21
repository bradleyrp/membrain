#!/bin/bash

#---INSTRUCTIONS
#-------------------------------------------------------------------------------------------------------------

<<'notes'
notes

#---DEFINITIONS
#-------------------------------------------------------------------------------------------------------------

#---Always call interactive python. Modules load /etc/pythonstart. Always import directly to the namespace.
python="python -i"

#---COMMON FUNCTIONS
#-------------------------------------------------------------------------------------------------------------

function eg_aamd_gr_voronoi_single {
$python ./script-aamd-gr-lipid-ion-voronoi-v1.py \
	-o gr \
	-mode voronoi_bin \
	-r aamd \
	-c /home/rpb/worker-big/membrane-test-set/look-v509/system-input.gro \
	-t /home/rpb/worker-big/membrane-test-set/look-v509/md.part0033.skip10.pbcmol.xtc
	-pickle /home/rpb/worker-big/membrane-test-set/look-v509/
}

function eg_aamd_gr_standard_internal {
./membrain-tool.sh \
	-o gr \
	-p v509.part0033 \
	-c $BASEPATH/look-v509/system-input.gro \
	-t $BASEPATH/look-v509/md.part0033.skip10.pbcmol.xtc \
	-r aamd \
	-skip 10 \
	-gr_group1 resname PI2P and name P5 \
	-gr_group2 name CL
}

function eg_aamd_gr_standard {
$python script-gr-v1.py \
	-c $BASEPATH/look-v509/system-input.gro \
	-t $BASEPATH/look-v509/md.part0033.skip10.pbcmol.xtc \
	-r aamd \
	-skip 10
}

function eg_aamd_gr_1d {
$python script-gr-lipid-ion-1d.py \
	-c ${trajpaths[0]} \
	-t ${trajpaths[1]} \
	-r aamd \
	-skip 20
}

function eg_meso_precomputed {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform regular \
	-prefix "membrane_XYZ-" \
	-suffix "\-ip.xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/meso-square-ncmembrane-precomputed \
	-nbase 52 \
	-length 1 \
	-rounder 10
}

function eg_meso_nc {
./membrain-tool.sh \
	-o undulations \
	-r meso \
	-xyzform square \
	-prefix "membrane_" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/DUMP \
	-nbase 52 \
	-start 90 \
	-end 99 \
	-viz \
	-length 1 \
	-rounder 10
}

function eg_meso_rect_k20_s40 {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform rect \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/meso-rect-controls/kappa-20/rep-0/equilibrate \
	-nbase 40 \
	-start 40 \
	-skip 1
}

function eg_meso_rect_k40_s40 {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform rect \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/meso-rect-controls/kappa-40/rep-0/equilibrate \
	-nbase 40 \
	-start 40 \
	-skip 1
}

function eg_meso_rect_k10_s40 {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform rect \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/meso-rect-controls/kappa-10/rep-0/equilibrate \
	-nbase 40 \
	-start 40 \
	-skip 1
}

function eg_cgmd_v600 {
./membrain-tool.sh \
	-viz \
	-o undulations \
	-mode best \
	-r cgmd \
	-c /home/rpb/worker-big/look-v600/system-input.gro \
	-t /home/rpb/worker-big/look-v600/md.part0006.correct.xtc \
	-skip 10
}

function eg_meso_rect_v2002_bare {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform rect \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/look-v2002/bare-rep-0/equilibrate \
	-nbase 22 \
	-start 1500 \
	-skip 1
}

function eg_meso_rect_v2002_anis {
./membrain-tool.sh \
	-o undulations \
	-viz \
	-r meso \
	-xyzform rect \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/look-v2002/anis-rep-0/equilibrate \
	-nbase 22 \
	-skip 1
}

function eg_meso_nc2 {
./membrain-tool.sh \
	-o undulations \
	-r meso \
	-xyzform square2 \
	-prefix "conf-" \
	-suffix ".xyz" \
	-dir /home/rpb/worker-big/membrane-test-set/blength-1p1 \
	-nbase 44 \
	-start 90 \
	-end 99 \
	-viz \
	-length 1.1 \
	-rounder 1
}


#---MAIN
#-------------------------------------------------------------------------------------------------------------

#eg_cgmd_v600
#eg_meso_precomputed
#eg_meso_nc
#eg_meso_nc2
#eg_meso_rect_k40_s40
#eg_meso_rect_k20_s40
#eg_meso_rect_k10_s40
#eg_meso_rect_v2002_bare
#eg_meso_rect_v2002_anis
