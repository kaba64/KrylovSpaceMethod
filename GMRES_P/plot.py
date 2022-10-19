#!/usr/loca/bin/python
from matplotlib import rc
from matplotlib import rc
import matplotlib.pyplot as plt
from numpy import *
from pylab import *

G_data = loadtxt("G_solution.txt");

#data = loadtxt("rh.txt");
m = max(G_data[:,0])+1;
n = max(G_data[:,1])+1

m = int(m)
n = int(n)

X = zeros([m, n]);
Y = zeros([m, n]);

G_field = zeros([m, n]);

for q in arange(0,m*n,1):
        i = int(G_data[q,0])
        j = int(G_data[q,1])
        G_field[i,j] = G_data[q,4]
        x = G_data[q,2]
        y = G_data[q,3]
        X[i, j] = x
        Y[i, j] = y

plt.figure(1)
plt.contourf(X,Y,G_field,20,cmap=cm.jet)
plt.colorbar()
plt.show()
#plt.title("GMRES solution field",fontsize=20)
