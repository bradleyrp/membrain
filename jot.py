#!/usr/bin/python

zone = zones[z]
print 1
inzone = array([2==(1*(ionstraj[ionsel][:,i,2]>zone[0])+1*(ionstraj[ionsel][:,i,2]<zone[1])) for i in range(nframes)])
print 2
inzonesliced = [[mean(inzone[i:i+d+1],axis=0) for i in range(nframes-d)] for d in range(nframes)]
print 3
starttime = time.time()
bigmask = [numpy.ma.masked_greater_equal(inzonesliced[i],occupancy) for i in range(nframes)]
print time.time()-starttime
print 4
mastermsd = [array((1*numpy.ma.mask_rows(bigmask[i]).mask+1*numpy.ma.mask_cols(bigmask[i]).mask)) for i in range(nframes)]
