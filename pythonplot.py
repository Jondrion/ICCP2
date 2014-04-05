from __future__ import division, print_function
from matplotlib import pyplot
from scipy.optimize import curve_fit
import numpy as np

filedata=np.genfromtxt("polymerdata.txt")



lengths=np.unique(filedata[:,0])

errEtoE=np.zeros(len(lengths))
avrEtoE=np.zeros(len(lengths))
errGyradius=np.zeros(len(lengths))
avrGyradius=np.zeros(len(lengths))
datasize=np.zeros(len(lengths))

i=0
for l in lengths:
	lengthdata=filedata[filedata[:,0]==l]

	datasize[i]=len(lengthdata[:,1])
	
	avrEtoE[i]=np.average(lengthdata[:,2], weights=lengthdata[:,1])
	errEtoE[i]=(np.average((lengthdata[:,2] - avrEtoE[i])**2, weights=lengthdata[:,1]))/(datasize[i]-1)
	
	avrGyradius[i]=np.average(lengthdata[:,3], weights=lengthdata[:,1])
	errGyradius[i]=np.average((lengthdata[:,3] - avrGyradius[i])**2, weights=lengthdata[:,1])/(datasize[i]-1)
	datasize[i]=len(lengthdata[:,1])
	del lengthdata
	i=i+1


# define function to fit the end to end distance
def fitfunc(x,a):
	return a*(x-1)**(1.5)

# define function to fit the gyradius
def fitfunc2(x,a,b):
	return a*(x-1)**(b)

# fit the data to the theoretical function for end to end distance
a1,b1=curve_fit(fitfunc,lengths,avrEtoE)
theo=fitfunc(lengths,a1)


# Determinde the fitted function for gyradius
a2,b2=curve_fit(fitfunc2,lengths,avrGyradius)

theo2=fitfunc2(lengths,a2[0],a2[1])

print("Fitconstant EtoE:",a1)
print("Fitconstants Gyradius:",a2)


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
ax.legend(bbox_to_anchor=(0.8, 1.05))

pyplot.draw()

fig2 = pyplot.figure("Gyradius")
ax2 = fig2.add_subplot(111)
ax2.set_xlim(1,np.max(lengths))

ax2.errorbar(lengths, avrGyradius,yerr=errGyradius**(1/2),fmt='ro', label="Gyradius")
ax2.plot(lengths, theo2, 'g--', label="Fitted Gyradius")
ax2.legend(bbox_to_anchor=(0.8, 1.05))

pyplot.draw()

pyplot.ion()
pyplot.show()
raw_input('Press <ENTER> to continue')
