#---file was formerly called compute_membrane_prop_vtu.py

def parse_membrane_data(filename):
	"parse the vtu file in xml format for membrane coordinates"
	tree=ET.parse(filename)
	root=tree.getroot()
	for d in root.iter('UnstructuredGrid'):
		coord=d.find('Piece').find('Points').find('DataArray').text.split()
		coord=map(float,coord)                           # convert string to float
		n=coord.__len__()/3
		vcoord=(numpy.array(coord).reshape(n,3))

		triangle=d.find('Piece').find('Cells').find('DataArray').text.split()
		triangle=map(int,triangle)
		n1=triangle.__len__()/3
		trilist=(numpy.array(triangle).reshape(n1,3))
		print 'File has ',n, 'vertices and ',n1,' triangles'
	data=[vcoord,trilist]
	return data

def compute_triangle_area(coord,trilist):
	area,proj_area=0.0,0.0
	for i in trilist:
		v12x,v12y,v12z=coord[i[1],0]-coord[i[0],0],coord[i[1],1]-coord[i[0],1],coord[i[1],2]-coord[i[0],2]
		v13x,v13y,v13z=coord[i[2],0]-coord[i[0],0],coord[i[2],1]-coord[i[0],1],coord[i[2],2]-coord[i[0],2]
		tparea=(v12x*v13y-v13x*v12y)*0.5
		ax,ay,az=(v12y*v13z-v13y*v12z),(v13x*v12z-v12x*v13z),(v12x*v13y-v13x*v12y)
		tarea=0.5*numpy.sqrt(ax**2+ay**2+az**2)
		area += tarea
		proj_area += tparea
	print area,proj_area
		
#!/usr/bin/env python
import numpy,glob,sys
import xml.etree.ElementTree as ET
import scipy,os

fname=sys.argv[1]
#os.chdir(sys.argv[1])
#filestring='*.vtu'
#filelist=glob.glob(filestring)
#print filelist
#for fname in filelist:
print 'reading file ',fname
print '------------------------------------>>>>>>>>>>>>>>>>'
data=parse_membrane_data(fname)
compute_triangle_area(data[0],data[1])

	

