#!/usr/bin/python

#---calculations
mset = msets[0]
vecs = mean(mset.vecs,axis=0)
m,n = mset.griddims
getgrid = array([[[i,j] for j in linspace(0,vecs[1]/mset.lenscale,n)] for i in linspace(0,vecs[0]/mset.lenscale,m)])
if 'tmp' not in globals(): tmp = array(msets[0].surf)
msets[0].surf = []
for i in range(len(tmp)):
	msets[index_md].surf.append(1*tmp[i])
maxpos = unravel_index((mean(msets[0].surf,axis=0)/msets[0].lenscale).argmax(),msets[0].griddims)
maxposgrid = [maxpos[i]/vecs[i]*msets[0].griddims[i] for i in range(2)]
hypo = [0.04,0.4,0.5,20,8,0]
params = [0,hypo[0],vecs[0]*hypo[1]/mset.lenscale,vecs[1]*hypo[2]/mset.lenscale,hypo[3],hypo[4],hypo[5]]
c0hypo = array([[gauss2d(params,getgrid[i,j,0],getgrid[i,j,1]) for j in range(n)] for i in range(m)])
mscs[index_md].calculate_mode_coupling(msets[index_md],[c0hypo for i in range(len(mset.surf))])

#---figure
fig = plt.figure()
gs = gridspec.GridSpec(3,4,wspace=0.0,hspace=0.0)

#---plot curvature field
extrem = max([mscs[i].t1d[1][:,1].max() for i in range(2)])
extrem = 0.1
ax = plt.subplot(gs[0,0])
ax.imshow(mean(mscs[0].c0s,axis=0),vmax=extrem,vmin=0.,cmap=mpl.cm.binary)
ax.set_title('CGMD')
ax = plt.subplot(gs[0,1])
ax.set_title('MESO')
im = ax.imshow(mean(mscs[1].c0s,axis=0),vmax=extrem,vmin=0.,cmap=mpl.cm.binary)
axins = inset_axes(ax,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,
	borderpad=0)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
plt.setp(axins.get_yticklabels())
plt.setp(axins.get_xticklabels())
axins.set_ylabel(r'$\left\langle C_0 \right\rangle (\mathrm{{nm}^{-1}})$',
	rotation=270)

#---plot average structure
vmax = max([mean(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
vmin = min([mean(msets[i].surf,axis=0).min()/msets[i].lenscale for i in range(2)])
extrem = max(abs(vmax),abs(vmin))
vmax,vmin = extrem,-extrem
ax = plt.subplot(gs[1,0])
ax.imshow(mean(msets[0].surf,axis=0)/msets[0].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r)
ax = plt.subplot(gs[1,1])
im = ax.imshow(mean(msets[1].surf,axis=0)/msets[1].lenscale,vmin=vmin,vmax=vmax,cmap=mpl.cm.RdBu_r)
axins = inset_axes(ax,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,
	borderpad=0)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
plt.setp(axins.get_yticklabels())
plt.setp(axins.get_xticklabels())
axins.set_ylabel(r'$\left\langle z(x,y)\right\rangle (\mathrm{nm})$',rotation=270)

#---plot standard deviations
extrem = max([std(msets[i].surf,axis=0).max()/msets[i].lenscale for i in range(2)])
ax = plt.subplot(gs[2,0])
ax.imshow(std(msets[0].surf,axis=0)/msets[0].lenscale,vmin=0.,vmax=extrem,cmap=mpl.cm.jet)
ax = plt.subplot(gs[2,1])
im = ax.imshow(std(msets[1].surf,axis=0)/msets[1].lenscale,vmin=0.,vmax=extrem,cmap=mpl.cm.jet)
axins = inset_axes(ax,width="5%",height="100%",loc=3,
	bbox_to_anchor=(1.,0.,1.,1.),
	bbox_transform=ax.transAxes,
	borderpad=0)
cbar = plt.colorbar(im,cax=axins,orientation="vertical")
plt.setp(axins.get_yticklabels())
plt.setp(axins.get_xticklabels())
axins.set_ylabel(r'$\left\langle \left(z-\overline{z}\right)^{2} \right\rangle (\mathrm{{nm}^2})$',
	rotation=270)

#---spectrum plot
axl = plt.subplot(gs[0,3])
m = 1
axl = plt.subplot(133)
axl.set_title(r'$\left\langle C_{0,q} h_{-q} \right\rangle $')
axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1]/msets[m].lenscale,marker='.',color='r',s=40,label='CGMD')
m = 0
axl.scatter(mscs[m].t1d[1][:,0],mscs[m].t1d[1][:,1]/msets[m].lenscale,marker='.',color='b',s=40,label='MESO')
axl.set_ylim((10**-11,10**-1))
axl.set_xlim((0.06,1.5))
axl.set_yscale('log')
axl.set_xscale('log')
axl.grid(True)
axl.yaxis.set_ticks_position("right")
axl.legend(loc='lower left')
gs.tight_layout(fig,h_pad=0.,w_pad=1.0)
plt.savefig(pickles+'fig-bilayer-couple-c0qhq-compare.png',bbox_inches='tight')
plt.show()

