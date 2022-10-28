#!/usr/loca/bin/python
from matplotlib import rc
from matplotlib import rc
import matplotlib.pyplot as plt
from numpy import *
from pylab import *
 
a_data = loadtxt("a_solution.txt");
b_data = loadtxt("b_solution.txt");
c_data = loadtxt("c_solution.txt");
 
m = max(a_data[:,0])+1;
n = max(a_data[:,1])+1

m = int(m)
n = int(n)

X = zeros([m, n]);
Y = zeros([m, n]);

a_field = zeros([m, n]);
b_field = zeros([m, n]);
c_field = zeros([m, n]);

for q in arange(0,m*n,1):
        i = int(a_data[q,0])
        j = int(a_data[q,1])
        a_field[i,j] = a_data[q,4]
        b_field[i,j] = b_data[q,4]
        c_field[i,j] = c_data[q,4]
        x = a_data[q,2]
        y = a_data[q,3]
        X[i, j] = x
        Y[i, j] = y

plt.figure(1)
plt.contourf(X,Y,a_field,levels=20, vmin=-1.0, vmax=1.0,cmap=cm.jet)
plt.colorbar()
plt.title("a",fontsize=20)
#plt.savefig('0.25(s).png')
plt.figure(2)
plt.contourf(X,Y,b_field,levels=20, vmin=-1.0, vmax=1.0,cmap=cm.jet)
plt.colorbar()
plt.title("b",fontsize=20)
#plt.savefig('2(s).png')
plt.figure(3)
plt.contourf(X,Y,c_field,levels=20, vmin=-1.0, vmax=1.0,cmap=cm.jet)
plt.colorbar()
plt.title("c",fontsize=20)
plt.show()
