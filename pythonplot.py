from __future__ import division, print_function
from matplotlib import pyplot
import numpy as np

filedata=np.genfromtxt("polymerdata.txt")

sorteddata=filedata[np.argsort(filedata[:,0])]

lengths=np.unique(sorteddata[:,0])

varEtoE=np.zeros(len(lengths))
avrEtoE=np.zeros(len(lengths))
varGyradius=np.zeros(len(lengths))
avrGyradius=np.zeros(len(lengths))
datasize=np.zeros(len(lengths))

i=0
for l in lengths:
	lengthdata=sorteddata[sorteddata[:,0]==l]
	
	avrEtoE[i]=np.average(lengthdata[:,2], weights=lengthdata[:,1])
	varEtoE[i]=np.average((lengthdata[:,2] - avrEtoE[i])**2, weights=lengthdata[:,1])
	
	avrGyradius[i]=np.average(lengthdata[:,3], weights=lengthdata[:,1])
	varGyradius[i]=np.average((lengthdata[:,3] - avrGyradius[i])**2, weights=lengthdata[:,1])
	datasize[i]=len(lengthdata[:,1])
	del lengthdata
	i=i+1

theo=(lengths-1)**(1.5)

fig = pyplot.figure("Polymerplot")
ax = fig.add_subplot(111)
ax.set_xlim(1,np.max(lengths))
ax.set_xlabel('Length of the polymer $N$')
ax.set_ylabel('End to End distance $R^2$')
ax.set_xscale('log')
ax.set_yscale('log')

ax.errorbar(lengths, avrEtoE,yerr=varEtoE**(1/2),fmt='o', label="End to End distance")
ax.plot(lengths, theo, 'r--', label="Theoretical End to End")
ax.plot(lengths, datasize, 'go', label="Datasize")
ax.legend(bbox_to_anchor=(0.8, 1.05))


# fig2 = pyplot.figure("Polymerplot")
# ax2 = fig2.add_subplot(111)
# ax2.set_xlim(1,np.max(lengths))
# ax2.set_xlabel('Length of the polymer $N$')
# ax2.set_ylabel('Gyradius $G$')
# ax2.set_xscale('log')
# ax2.set_yscale('log')

# ax2.errorbar(lengths, avrGyradius,yerr=varGyradius**(1/2),fmt='o', label="Gyradius")

# ax2.legend(bbox_to_anchor=(0.8, 1.05))
pyplot.show()