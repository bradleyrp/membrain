#!/usr/bin/python

for i in range(len(hists)):
	plt.scatter(range(len(hists[i])), hists[i])
plt.show()

