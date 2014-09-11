#!/usr/bin/python

from ModeCouple import *

#---? COMMON FILES
execfile('header.py')
execfile('script_header.py')

#---connect
if 'conn' not in globals(): 
	conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	#try: conn = psycopg2.connect("dbname='membrain_simbank' user='rpb' host='localhost' password=''")
	#except: raise Exception('except: cannot connect to database')
	cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
	
#---retrieve a set of experiments for comparison
cur.execute('SELECT * from mesosims')
select_mesosims = [dict(i) for i in cur.fetchall()]
cur.execute('SELECT * from dataref_structure')
select = [dict(i) for i in cur.fetchall()]
#---note that I am filtering the table in python and not postgresql
select = [i for i in select if all([
	i[j] == select_criteria_meso[j] for j in select_criteria_meso.keys()])]
#---populate analysis descriptors from the database
for params in select:
	ind = where([i['id']==params['parent_mesosims'] for i in select_mesosims])[0][0]
	combo = dict(params.items() + select_mesosims[ind].items())
	rundirnum = int(combo['rundir'])
	key = params['callsign']+'-rundir-'+str(rundirnum)
	#---always fix uppercase naming when importing to python
	if 'c_0' in combo.keys(): combo['C_0'] = combo['c_0']
	combo['detail_name'] = combo['shortname']
	analysis_descriptors[key] = combo

if 'compare_parameters_grid' in routine or 1:

	#---settings
	master_details_listing = []
	do_couple_plot = False
	whichtest = ['fullscan','oneplot'][1]
	if whichtest == 'oneplot':
		cgmd_avail_list = ['v616-210000-310000-200','v550-300000-400000-200']
		meso_avail_list = ['v2016-rundir-0','v2016-rundir-4','v2016-rundir-8','v2016-rundir-12']
	elif whichtest == 'fullscan':
		cgmd_avail_list = cgmd_avail
		meso_avail_list = [i for i in analysis_descriptors.keys() if i not in cgmd_avail]

	#---master loop
	master_mscs = []
	for batch_cgmd in cgmd_avail_list:
		for batch_meso in meso_avail_list:

			match_scales = analysis_names = plot_reord = [batch_cgmd,batch_meso]
			status('status: comparing '+batch_cgmd+' '+batch_meso)
			bigname = '-'.join(analysis_names)

			#---BEGIN MATCHING SECTION

			#---allocate if empty
			if 'msets' not in globals(): msets = []
			if 'mscs' not in globals(): mscs = []
			if 'collect_c0s' not in globals(): collect_c0s = []

			#---load and interpolate
			lenscale = 1.0
			for a in analysis_names:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				if 'mset' in globals(): del mset
				mset = MembraneSet()
				if simtype == 'meso':
					params = analysis_descriptors[a]
					#---? type on rundir is int but returns as a float
					pklname = 'pkl.structures.meso.'+\
						params['callsign']+'-'+\
						'C_0-'+str(params['C_0'])+'-'+\
						'len-'+str(params['lenscale'])+'-'+\
						'rundir-'+str(int(params['rundir']))+\
						'.pkl'
					print 'unpickle '+str(pklname)
					mset = unpickle(pickles+pklname)
					print array(mset.surf[0]).max()
					c0s = mset.getdata('c0map').data
					collect_c0s.append(c0s)
					msets.append(mset)
				elif simtype == 'md':
					status('status: loading from MD')
					if 'mset' in globals(): del mset
					mset = unpickle(pickles+locate)
					msets.append(mset)
					#---here we set the hypothetical curvature equal to the 
					#---...induced curvature at the mesoscale
					c0ask = (analysis_descriptors[analysis_names[1]])['C_0']
					hypo[0] = c0ask
					#---compute hypothetical curvature field for the MD simulation
					if hypo != None:
						vecs = mean(mset.vecs,axis=0)
						m,n = mset.griddims
						#---getgrid is xyz points in nm
						getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] 
							for i in linspace(0,vecs[0]/mset.lenscale,m)])
						#---convert everything to nm
						#---recall z0,c0,x0,y0,sx,sy,th = params
						#---params sets C0 in native units, x0,y0 in proportional units, and sx,sy in nm
						params = [0,hypo[0],
							vecs[0]*hypo[1]/mset.lenscale,
							vecs[1]*hypo[2]/mset.lenscale,
							hypo[3],hypo[4],
							hypo[5]]
						c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)]
							for i in range(m)])
						collect_c0s.append([c0hypo for i in range(len(mset.surf))])
					else:
						collect_c0s.append([])
				elif simtype == 'meso_precomp':
					if 'mset' in globals(): del mset
					mset = unpickle(pickles+locate)
					#---for precomp simulations these units are relative to a0 and already in a reglar grid
					collect_c0s.append(mset.getdata('c0map').data)
					msets.append(mset)

			#---match mesoscale length scales to an MD simulation
			if match_scales != None:
				ref_ind = analysis_names.index(match_scales[0])
				move_ind = analysis_names.index(match_scales[1])
				#---matching the average box vectors here by average and not maximum
				lenscale = mean(mean(msets[move_ind].vecs,axis=0)[:2])/\
					(mean(mean(msets[ref_ind].vecs,axis=0)[:2])/msets[ref_ind].lenscale)
				for a in analysis_names:
					if (analysis_descriptors[a])['simtype'] == 'meso' or \
						(analysis_descriptors[a])['simtype'] == 'meso_precomp':
						#analysis_descriptors[a]['lenscale'] = lenscale
						for i in analysis_descriptors[a]: 
							if i != 'lenscale': vars()[i] = (analysis_descriptors[a])[i]
						msets[analysis_names.index(a)].lenscale = lenscale
				'''
				#---here we set the hypothetical curvature equal to the induced curvature at the mesoscale
				hypo[0] = c0ask
				#---reset the hypothetical C0 field according to the new scaling 
				#---... this replaces the calculation above
				for a in analysis_names:
					for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
					if analysis_descriptors[a]['simtype'] == 'md':
						anum = analysis_names.index(a)
						mset = msets[anum]
						vecs = mean(mset.vecs,axis=0)
						m,n = mset.griddims
						#---getgrid is xyz points in nm
						getgrid = array([[[i,j] for j in linspace(0,3*vecs[1]/mset.lenscale,3*n)] 
							for i in linspace(0,3*vecs[0]/mset.lenscale,3*m)])
						#---needs checked, "key step used to have a 0.5*hypo[0] here 
						#---...possibly due to convention"
						params = [0,
							hypo[0]*msets[1].lenscale/mset.lenscale,
							vecs[0]*(1+hypo[1])/mset.lenscale,
							vecs[1]*(1+hypo[2])/mset.lenscale,
							sqrt(r_2)/msets[1].lenscale/sqrt(2),
							sqrt(r_2)/msets[1].lenscale/sqrt(2),
							hypo[5]]
						#---handle curvature fields at the mesoscale which cross the PBC boundary
						c0hypo_nopbc = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) 
							for j in range(3*n)] for i in range(3*m)])
						c0hypo = [[max([c0hypo_nopbc[i+sh[0]*m,j+sh[1]*n] for sh in [[k,l] 
							for k in range(3) for l in range(3)]]) for j in range(n)] for i in range(m)]
						collect_c0s[anum] = [c0hypo for i in range(len(mset.surf))]
				'''
				construct_hypofield(c0ask)

			#---calculate coupled modes
			for a in analysis_names:
				for i in analysis_descriptors[a]: vars()[i] = (analysis_descriptors[a])[i]
				m = analysis_names.index(a)
				if 'msc' in globals(): del msc
				msc = ModeCouple()
				msc.calculate_mode_coupling(msets[m],collect_c0s[m])
				mscs.append(msc)

			test_kappas = [18,20,22,24,26]
			test_sigmas = [-0.5,0.,0.5,1.0,2.0]
			colors = ['g','b']
			simlabel = ['cgmd','meso']

			#---separate residuals
			if 1:
				#---optimize
				popts = []
				for m in range(2):
					p_opt = leastsq(equipart_resid,array([20.0,0.0]))
					print 'optimized, '+simlabel[m]+' = '+str(list(p_opt[0]))
					test_kappas.append(p_opt[0][0])
					test_sigmas.append(p_opt[0][1])
					popts.append(p_opt)
				master_details_listing.append({
					'batch_meso':batch_meso,
					'batch_cgmd':batch_cgmd,
					'kappa_cgmd':popts[0][0][0],
					'sigma_cgmd':popts[0][0][1],
					'kappa_meso':popts[1][0][0],
					'sigma_meso':popts[1][0][1],
					'resid_meso':mean([i**2 for i in equipart_resid([popts[1][0][0],popts[1][0][1]])]),
					'resid_cgmd':mean([i**2 for i in equipart_resid([popts[0][0][0],popts[0][0][1]])]),
					})
			#---deprecated previously used with equipart_resid(both=True) to do combination residual
			if 0:
				p_opt,cov,infodict,mesg,ie = leastsq(equipart_resid,array([20.0,0.0]),
					full_output=True)
				master_details_listing.append({
					'batch_meso':batch_meso,
					'batch_cgmd':batch_cgmd,
					'kappa':p_opt[0],
					'sigma':p_opt[1],
					'resid':mean([i**2 for i in equipart_resid([p_opt[0],p_opt[1]])]),
					})
			master_mscs.append(mscs)

			#---grid of comparisons
			if whichtest == 'oneplot': 
				gs = gridspec.GridSpec(len(test_kappas),len(test_sigmas),wspace=0.0,hspace=0.0)
				fig = plt.figure(figsize=(12,12))
				for ki in range(len(test_kappas)):
					for si in range(len(test_sigmas)):
						ax = fig.add_subplot(gs[ki,si])
						resid = [0,0]
						for m in range(2):
							newkappa = test_kappas[ki]
							newsigma = test_sigmas[si]
							resid[m] = mean([i**2 for i in equipart_resid([newkappa,newsigma])])
							mset = msets[m]
							qmagst = mscs[m].qmagst
							area = double(mean([mset.vec(i)[0]*mset.vec(i)[1] 
								for i in mset.surf_index])/mset.lenscale**2)
							scalefac = newkappa*area
							tsum2d = scalefac*(mscs[m].t2d[0]*qmagst**4-\
								mscs[m].t2d[1]*qmagst**2-mscs[m].t2d[2]*qmagst**2+mscs[m].t2d[3])+\
								1*qmagst**2*mscs[m].t2d[0]*newsigma*area*msets[1].lenscale**2
							xdat = collapse_spectrum(mscs[m].qmagst,mscs[m].qmagst)
							ydat = collapse_spectrum(mscs[m].qmagst,tsum2d)
							ax.plot(xdat,ydat,'-',lw=2,label=shortname,c=colors[m],alpha=1.0)
						ax.text(0.05,0.95,
							'$\kappa=$'+'{:.1f}'.format(newkappa)+'\n'+\
							'$\gamma=$'+'{:.1f}'.format(newsigma)+'\n'+\
							'{:.2f}'.format(resid[0])+','+'{:.2f}'.format(resid[1]),
							transform=ax.transAxes,va='top')
						if ki == len(test_kappas)-2 and si == len(test_kappas)-2:
							ax.text(0.05,0.15,'CGMD,opt',transform=ax.transAxes,va='top')
						if ki == len(test_kappas)-1 and si == len(test_kappas)-1:
							ax.text(0.05,0.15,'MESO,opt',transform=ax.transAxes,va='top')
						ax.set_xscale('log')
						ax.set_yscale('log')
						ax.set_ylim((0.1,10))
						ax.set_xlim((0.01,10))			
						ax.axhline(y=1,xmin=0.,xmax=1,lw=2,color='k')
						ax.axvline(x=1.0,ymin=0.,ymax=1,lw=2,color='k')
						ax.grid(True)
						ax.get_xaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower',nbins=3))
						ax.get_yaxis().set_major_locator(mpl.ticker.MaxNLocator(prune='lower',nbins=3))
						if ki < len(test_kappas)-1: ax.set_xticklabels([])
						if si > 0: ax.set_yticklabels([])
				plt.suptitle(analysis_descriptors[batch_cgmd]['label']+' vs '+'meso $C_0='+\
					'{:.3f}'.format(analysis_descriptors[batch_meso]['c_0']*msets[1].lenscale)+\
					'\:{(nm)}^{-1}$',fontsize=fsaxlabel+2)
				plt.savefig(pickles+'fig-coupling-signatures-'+batch_cgmd+'-'+batch_meso+'.png',dpi=300)
				plt.clf()		
			#---clean up the loop
			del mscs,msets
