#!/usr/bin/python -i

from membrainrunner import *

location = ''
execfile('locations.py')

pickleset = ['pkl.apl-cells.v533.md.part0003.2000-5000-4.pkl',
	'pkl.apl-cells.v534.md.part0003.2000-5000-4.pkl']
names = ['PI(3,5)P2, Mg','PI(3,5)P2, Ca']
figname = 'v533.v534'

do_P35P_area_compare = 1
do_addon_visualize = 0

#---MAIN
#-------------------------------------------------------------------------------------------------------------

#---Compare the areas per lipid for P35P residues and plot.
if do_P35P_area_compare:
	msets = []
	target_data = []
	target_residues = []
	for i in range(len(pickleset)):
		mset = unpickle(pickles+pickleset[i])
		msets.append(mset)
		#---Identify the P35P residues, first by absolute number, then by index in the monolayer
		#---Beware, hack: only works if you are studying one lipid type in a single monolayer.
		target_residues_abs = msets[i].monolayer_by_resid[0][msets[i].resnames.index('P35P')]
		target_residues.append([msets[0].monolayer_residues[0].index(i) for i in target_residues_abs])
		rawdat = mset.getdata('cells').get(['type','areas'])
		#---Only save the areas for the selected residues
		target_data.append([rawdat[j][k][i] for i in target_residues[-1] 
			for j in range(len(rawdat)) for k in range(2)])

	clrs = brewer2mpl.get_map('Set1', 'qualitative', 5).mpl_colors
	hists = []
	fig = plt.figure()
	ax = plt.subplot(111)
	for i in range(len(target_data)):
		hist0 = numpy.histogram(target_data[i],bins=30,normed=True)
		hists.append(hist0)
		c = clrs[i%len(clrs)]
		ax.plot(hist0[1][1:],hist0[0],'-',label=names[i],lw=2,color=c)
	plt.xlabel('lipid area (\AA)')
	plt.ylabel('frequency')
	plt.legend(loc=2)
	plt.savefig(pickles+'fig-apl-cells-'+figname+'.png', dpi=300)
	plt.show()
	
#---Temporary addition to script-aamd-apl.py
if do_addon_visualize:
	msets = [mset]
	tmp = mset.getdata('cells').get(['frame',0,'monolayer',0,'type','points'])
	tmp2 = mset.wrappbc(tmp,vecs=mset.vecs[0],mode='grow')
	meshplot(tmp2)
	target_residues_abs = msets[0].monolayer_by_resid[0][msets[0].resnames.index('CHL1')]
	target_residues = [msets[0].monolayer_residues[0].index(i) for i in target_residues_abs]
	meshpoints([tmp[i] for i in target_residues],scale_factor=[10. for i in range(len(target_residues))])
	target_residues_abs = msets[0].monolayer_by_resid[0][msets[0].resnames.index('P35P')]
	target_residues = [msets[0].monolayer_residues[0].index(i) for i in target_residues_abs]
	meshpoints([tmp[i] for i in target_residues],scale_factor=[10. for i in range(len(target_residues))],color=(1,1,1))
