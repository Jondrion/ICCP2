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


i=0
for l in lengths:
	lengthdata=sorteddata[sorteddata[:,0]==l]
	
	avrEtoE[i]=np.average(lengthdata[:,2], weights=lengthdata[:,1])
	varEtoE[i]=np.average((lengthdata[:,2] - avrEtoE[i])**2, weights=lengthdata[:,1])
	
	avrGyradius[i]=np.average(lengthdata[:,3], weights=lengthdata[:,1])
	varGyradius[i]=np.average((lengthdata[:,3] - avrGyradius[i])**2, weights=lengthdata[:,1])
	i=i+1
	


theo=(lengths-1)**(1.5)

fig = pyplot.figure("Polymerplot")
ax = fig.add_subplot(111)
ax.set_xlim(1,np.max(lengths))
ax.set_xlabel('Length of the polymer $N$')
ax.set_ylabel('End to End distance $R^2$')
ax.set_xscale('log')
ax.set_yscale('log')
print(varEtoE**(1/2))
ax.errorbar(lengths, avrEtoE**2,yerr=varEtoE,fmt='o', label="End to End distance")
ax.plot(lengths, theo, 'r--', label="Theoretical End to End")
ax.legend(bbox_to_anchor=(0.8, 1.05))
pyplot.show()