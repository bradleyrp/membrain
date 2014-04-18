index = [i for i in range(len(msdats[anum])) if all([
	msdats[anum][i].getnote('cutoff') == cutoff,
	msdats[anum][i].getnote('zfiltdir') == zfiltdir,
	msdats[anum][i].getnote('decay_z0_min') == decay_z0_min,
	msdats[anum][i].getnote('decay_z0') == decay_z0,
	all(msdats[anum][i].getnote('tl')[msdats[anum][i].getnote('this_test')] == test)])]

tmp = array([[msdats[anum][i].getnote('cutoff') == cutoff,
	msdats[anum][i].getnote('zfiltdir') == zfiltdir,
	msdats[anum][i].getnote('decay_z0_min') == decay_z0_min,
	msdats[anum][i].getnote('decay_z0') == decay_z0,
	all(msdats[anum][i].getnote('tl')[msdats[anum][i].getnote('this_test')] == test)] for i in range(len(msdats[anum]))])

tmp2 =[msdats[anum][i].getnote('decay_z0_min') for i in range(len(msdats[anum]))]
