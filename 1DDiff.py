#!/usr/local/bin/python
from math import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
import os
import errno

## 1DDiff.py
## Remi Carmigniani
## 

tend=0.5
##################################################################################################################
############################			Useful functions		      ############################
##################################################################################################################
#Trap Integration
def trapInt(f,dx,nx):
	error = 0.0
	for i in range(0,nx-1):
		error = error+f[i]+f[i+1]	
	return error*.5*dx


##################################################################################################################
############################				End			      ############################
##################################################################################################################
#IMPORTANT COMMENT : To copy a list do not simply use the command unew = u this link the two togethers
L=2*pi
error=[]
## Series of N to try
N = [25,50,100]
#first calculate h = dt/dx =cte
dx = L/float(100-1)
dt = .1*dx*dx
h=dt/dx

for k in range(0,3):
        t=0
	## Discretization parameter 
	Nx=N[k]
	dx=L/float(Nx-1)	
	dt = 0.1*dx*dx
	h=dt/dx
	lx = 0.1
	
	u = [0 for i in xrange(Nx)]
	## build u
	for i in range(0,Nx):	
		u[i]=cos(i*dx)
	plt.show()
	while t<tend :
		t=t+dt
		uold= u[:]
		u[0]    = uold[Nx-2]*lx + (1.-2.*lx)*uold[0]    + uold[1]*lx
		u[Nx-1] = uold[Nx-2]*lx + (1.-2.*lx)*uold[Nx-1] + uold[1]*lx
		for i in range(1,Nx-1):	
			u[i] = uold[i-1]*lx+(1.-2.*lx)*uold[i]+uold[i+1]*lx
        #Error
        	
	f = [0 for i in xrange(Nx)] 
	for i in range(0,Nx):
			f[i] = u[i]-cos(i*dx)*exp(-t)
			f[i]=f[i]**2
	
        err=sqrt(trapInt(f,dx,Nx))
	error.append(err)
	print 'Calculate the error  for ' + repr(Nx) + ' : ' + repr(err)

plt.loglog(N, error)
plt.title('Error')
plt.savefig('Error.png')


plt.close()
plt.plot(u)
plt.show()
print 'Simulation Completed without error'

		
        	

       
	
	


