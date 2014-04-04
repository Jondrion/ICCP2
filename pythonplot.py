from __future__ import division, print_function
from matplotlib import pyplot
import numpy as np

filedata=np.genfromtxt("polymerdata.txt")

sorteddata=filedata[np.argsort(filedata[:,0])]

lengths=np.unique(sorteddata[:,0])

errEtoE=np.zeros(len(lengths))
avrEtoE=np.zeros(len(lengths))
errGyradius=np.zeros(len(lengths))
avrGyradius=np.zeros(len(lengths))
datasize=np.zeros(len(lengths))

i=0
for l in lengths:
	lengthdata=sorteddata[sorteddata[:,0]==l]

	datasize[i]=len(lengthdata[:,1])
	
	avrEtoE[i]=np.average(lengthdata[:,2], weights=lengthdata[:,1])
	errEtoE[i]=(np.average((lengthdata[:,2] - avrEtoE[i])**2, weights=lengthdata[:,1]))/(datasize[i]-1)
	
	avrGyradius[i]=np.average(lengthdata[:,3], weights=lengthdata[:,1])
	errGyradius[i]=np.average((lengthdata[:,3] - avrGyradius[i])**2, weights=lengthdata[:,1])/(datasize[i]-1)
	datasize[i]=len(lengthdata[:,1])
	del lengthdata
	i=i+1
print(datasize)
theo=(lengths-1)**(1.5)

fig = pyplot.figure("Polymerplot")
ax = fig.add_subplot(111)
ax.set_xlim(1,np.max(lengths))
ax.set_xlabel('Length of the polymer $N$')
ax.set_ylabel('End to End distance $R^2$')
ax.set_xscale('log')
ax.set_yscale('log')

ax.errorbar(lengths, avrEtoE,yerr=errEtoE**(1/2),fmt='o', label="End to End distance")
ax.plot(lengths, theo, 'r--', label="Theoretical End to End")
ax.plot(lengths, datasize, 'go', label="Datasize")
ax.legend(bbox_to_anchor=(0.8, 1.25))


fig2 = pyplot.figure("Polymerplot")
ax2 = fig2.add_subplot(111)
ax2.set_xlim(1,np.max(lengths))

ax2.errorbar(lengths, avrGyradius,yerr=errGyradius**(1/2),fmt='ro', label="Gyradius")

ax2.legend(bbox_to_anchor=(0.8, 1.05))
pyplot.show()