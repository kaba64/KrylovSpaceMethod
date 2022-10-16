#!/usr/loca/bin/python
from matplotlib import rc
from matplotlib import rc
import matplotlib.pyplot as plt
from numpy import *
from pylab import *

#C_data = loadtxt("C_solution.txt");
#G_data = loadtxt("G_solution.txt");
B_data = loadtxt("B_solution.txt");
#data = loadtxt("rh.txt");
m = max(B_data[:,0])+1;
n = max(B_data[:,1])+1

m = int(m)
n = int(n)

X = zeros([m, n]);
Y = zeros([m, n]);
#field = zeros([m, n]);
#G_field = zeros([m, n]);
B_field = zeros([m, n]);
for q in arange(0,m*n,1):
        i = int(B_data[q,0])
        j = int(B_data[q,1])
        #field[i,j] = C_data[q,4]
        #G_field[i,j] = G_data[q,4]
        B_field[i,j] = B_data[q,4]
        x = B_data[q,2]
        y = B_data[q,3]
        X[i, j] = x
        Y[i, j] = y

plt.figure(1)
plt.contourf(X,Y,B_field,20,cmap=cm.jet)
plt.colorbar()
plt.title("BIGCSTAB solution field",fontsize=20)
plt.show()
