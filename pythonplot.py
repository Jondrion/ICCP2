from __future__ import division, print_function
from matplotlib import pyplot
import numpy as np

filedata=np.genfromtxt("polymerdata.txt")

print(filedata)
sorteddata=filedata[np.argsort(filedata[:,0])]

lengths=np.unique(sorteddata[:,0])

stdEtoE=np.zeros(len(lengths))
avrEtoE=np.zeros(len(lengths))
stdGyradius=np.zeros(len(lengths))
avrGyradius=np.zeros(len(lengths))

i=0
for l in lengths:
	lengthdata=sorteddata[sorteddata[:,0]==l]
	
	stdEtoE[i]=np.std(lengthdata[:,1])
	avrEtoE[i]=np.average(lengthdata[:,1])
	stdGyradius[i]=np.std(lengthdata[:,2])
	avrGyradius[i]=np.average(lengthdata[:,2])
	i=i+1


theo=lengths**(1.5)

fig = pyplot.figure("Polymerplot")
ax = fig.add_subplot(111)
ax.set_xlim(0,len(lengths))
ax.set_xlabel('Length of the polymer $N$')
ax.set_ylabel('End to End distance $R^2$')
ax.plot(lengths, avrEtoE**2,'r-',label="End to End distance")
ax.plot(lengths, theo, 'b--', label="Theoretical End to End")
ax.legend(bbox_to_anchor=(1.1, 1.05))
pyplot.show()