#!/usr/bin/python

if 0:
	dat = hqs[0:1]
	#cm,cn = [int(round(i/2.-1)) for i in shape(dat[0])]
	m,n = shape(dat[0])
	print m,n
	v1 = array(dat)[:,slice((1 if m%2==0 else None),None),slice((1 if n%2==0 else None),None)]
	v2 = array(dat)[:,slice((-1 if m%2==1 else None),(0 if m%2==0 else None),-1),slice((-1 if n%2==1 else None),(0 if n%2==0 else None),-1)]
	print shape(v1)
	print shape(v2)
	qmags = array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) for j in range(0,n)] for i in range(0,m)])
	ax = plt.subplot(311)
	plt.imshow(array(abs(v1[0])).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	ax = plt.subplot(312)
	plt.imshow(array(abs(v2[0])).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	ax = plt.subplot(313)

	qmags = array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) for j in range(0,n)] for i in range(0,m)])
	qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]

	plt.imshow(array(abs(qmagst)).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	plt.show()
if 1:
	ax = plt.subplot(311)
	plt.imshow(array(mean(abs(hqsa*hqsb),axis=0)**2).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	ax = plt.subplot(312)
	plt.imshow(array(abs(mean(hqsb,axis=0))).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	ax = plt.subplot(313)

	qmags = array([[sqrt(((i-cm)/((Lx)/1.)*2*pi)**2+((j-cn)/((Ly)/1.)*2*pi)**2) for j in range(0,n)] for i in range(0,m)])
	qmagst = qmags[0:(-1 if m%2==0 else None),0:(-1 if n%2==0 else None)]

	plt.imshow(array(abs(qmagst)).T,interpolation='nearest',origin='lower',norm=mpl.colors.LogNorm());
	plt.show()
