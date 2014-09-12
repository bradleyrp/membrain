#!/usr/bin/python
if 0:
	tmp = [[ms.kqs[(int(u/m2)+int(v/n2))%m2,(u%m2+v%n2)%n2] for v in range(m2*n2)] for u in range(m2*n2)]
	tmpspecplotter(tmp,callsign+'-term-'+data[1]+'-homog',logplot=logplot,show=True)
	qis = array([[[i,j] for j in range(n2)] for i in range(m2)])
	qqp = (qis+qis)%[64,64]
	tmp = abs(fft.ifftshift(ms.kqs))
	tmp2 = array([[tmp[(i+k)%m2,(j+l)%n2] for j in range(n2) for l in range(n2)] for i in range(m2) for k in range(m2)])
	tmp3 = fft.fftshift(tmp2)
	tmp4 = outer(ms.kqs,ms.kqs)
	tmpspecplotter(tmp4,callsign+'-term-'+data[1]+'-homog',logplot=logplot,show=True)
if 0:
	cm,cn = [int(round(i/2.-1)) for i in shape(ms.msc.t2d[0])]
	st = time.time()
	kqqp = array([[(ms.kqs[i+k-cm,j+l-cn] if (i+k-cm<n2 and i+k-cm>0 and j+l-cn<m2 and j+l-cn>0) else 0) 
		for j in range(-n2/2,n2/2) for l in range(-n2/2,n2/2)] for i in range(-m2/2,m2/2) for k in range(-m2/2,m2/2)])
	print 1./60*(time.time()-st)
	tmpspecplotter(ms.kqs,callsign+'-term-'+data[1]+'-homog',logplot=logplot,show=True)
if 0:
	inds = array([[[int(u/n2)-int(v/n2),u%n2-v%n2] for v in range(m2*n2)] for u in range(m2*n2)])
	inds = array([[[int(u/n2)-int(v/n2),u%n2-v%n2] for v in range(m2*n2)] for u in range(m2*n2)])

kqs = zeros((m2,n2));kqs[cm,cn]=43;
inds = array([[[int(u/n2)-int(v/n2),u%n2-v%n2] for v in range(m2*n2)] for u in range(m2*n2)])
((1*(inds[...,0]>0)+1*(inds[...,1]>0)+1*(inds[...,0]<m2)+1*(inds[...,1]<n2))==4)
plt.imshow(inds[...,1],interpolation='nearest',origin='lower');plt.show()
plt.imshow(inds[...,0],interpolation='nearest',origin='lower');plt.show()
