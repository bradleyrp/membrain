#from MDAnalysis import *
#mset_get_charges = Universe('/home/rpb/worker-big/membrane-repository/look-v510/md.part0033.tpr','/home/rpb/worker-big/membrane-repository/look-v510/system-input.gro')
#fp = open('/home/rpb/pkl.v510.charges.pkl','w')
#sel = mset_get_charges.selectAtoms('all')
#pickle.dump(tmp.charges(),fp)
#fp.close()
import re
fp = open('/home/rpb/worker-big/membrane-repository/look-v510/log-gmxdump','r')
rawfile = []
for line in fp:
	rawfile.append(line)
fp.close()
for l in rawfile:
	if re.findall('q=',l):
		print l
