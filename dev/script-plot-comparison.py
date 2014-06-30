#!/usr/bin/python -i

from membrainrunner import *

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import brewer2mpl

#---Table of results from dimple fitting
#---Headers are: system,cutoff,meanmaxh,sigmax,sigmay,fitted,total_frames
curve_fit_data_header = ['system','cutoff (nm)','$H_{max}$','$\sigma_a$','$\sigma_b$','$N_{fit}$','$N_{tot}$']
curve_fit_data = [
	['4-ENTH',15,0.010,3.131,3.171,184,352],
	['1-ENTH',15,0.010,3.226,3.951,218,446],
	['control',15,0.007,3.640,3.454,381,1152],
	['control',10,0.009,3.550,3.098,335,1152],
	['4-ENTH, 1',15,0.012,2.911,2.876,184,352],
	['4-ENTH, 2',15,0.008,3.977,3.395,126,352],		
	['4-ENTH, 3',15,0.010,3.197,3.455,162,352],
	['4-ENTH, 4',15,0.009,3.680,3.820,133,352],
	['4-ENTH',10,0.010,3.268,3.415,145,352],
	['4-ENTH',20,0.011,3.097,3.127,221,352],
	['4-ENTH',25,0.009,3.287,3.330,222,352],
	['4-ENTH',30,0.009,3.005,3.264,216,352]]

#---Updated data
curve_fit_data = [
	['4-ENTH',15,0.014,34.132,28.120,184,352],
	['1-ENTH',15,0.014,35.071,47.514,218,446],
	['control',15,0.013,63.826,48.456,381,1152],
	['4-ENTH, 1',15,0.018,22.960,20.319,184,352],
	['4-ENTH, 2',15,0.012,60.562,46.386,126,352],		
	['4-ENTH, 3',15,0.014,42.037,37.008,162,352],
	['4-ENTH, 4',15,0.014,81.631,60.302,133,352]]

	
#---Print table to terminal. Check out collections.namedtuple to improve this
if 0:
	for i in curve_fit_data_header:
		print i+'\t',
	print '\n',
	for i in curve_fit_data:
		for j in i:
			print str(j)+'\t',
		print '\n',

data = curve_fit_data
data_t = [list(i) for i in (np.array(curve_fit_data).T)]
header = curve_fit_data_header
		
columns = tuple([i[0] for i in data])
rows = tuple(header[1:])

#---Define colors
colorcodes_all = brewer2mpl.get_map('paired','qualitative',10).mpl_colors
color_bar_dict = {'4-ENTH':1,'1-ENTH':9,'control':3,'4-ENTH, 1':0,'4-ENTH, 2':0,'4-ENTH, 3':0,'4-ENTH, 4':0}
bar_colors = [colorcodes_all[color_bar_dict[i]] for i in [j[0] for j in data]]
reorder_rows = [0,4,5,6,7,1,2,3,8,9,10,11]
reorder_rows = range(len(curve_fit_data))

#---Bar plot
bar_width=0.6
index = np.arange(len(data)) + 0.2
mpl.rcParams.update({'font.size': 14})
fig = plt.figure(figsize=(15,6))
ax = plt.subplot2grid((1,1),(0,0))
for r in range(len(reorder_rows)):
	row = reorder_rows[r]
	print r
	print row
	plt.bar(index[r], data[row][2],bar_width,color=bar_colors[row],alpha=0.8)
#---Table
datatable = ax.table(
	cellText=(data_t[1:]),
	loc='bottom',
	colLoc='center',
	rowLoc='center',
	rowLabels=header[1:],
	colLabels=data_t[0])

table_props = datatable.properties()
table_cells = table_props['child_artists']
for cell in table_cells: 
	cell.set_fontsize(12)
	cell.set_fontsize(12)
	cell.set_lw(0.8)
	cell.set_height(0.06)
	tmp = cell.get_text()
	tmp.set_horizontalalignment('center')

plt.subplots_adjust(left=0.2, bottom=0.4)
plt.xticks([])
maxy = round(max([i[2] for i in data])*1.2,3)
ax.set_ylim((0,maxy))
plt.title('Curvature Comparison')
ax.grid()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_position(('outward',20))
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
plt.show()
