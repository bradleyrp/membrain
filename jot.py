#!/usr/bin/python

if 0:
	rewrapz = [ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)]
	t = array(rewrapz)[0]

	rewrapz = array([ionstraj[i,:,2]-1*(array(ionstraj[i,:,2])>vecs[2])*vecs[2]+(array(ionstraj[i,:,2])<0.)*vecs[2] for i in range(nions)])
	binws = [linspace(0,i[2],14)[1] for i in mset_surf.vecs]
	t = array([[ionstraj[i,j,2]-1*(array(ionstraj[i,j,2])>mset_surf.vecs[j][2])*mset_surf.vecs[j][2]+(array(ionstraj[i,j,2])<0.)*mset_surf.vecs[j][2] for j in range(len(ionstraj[i]))] for i in range(len(ionstraj))])
	t1 = array([[int(round(t[i][j]/binws[j])) for j in range(len(t[i]))] for i in range(len(t))])
	

